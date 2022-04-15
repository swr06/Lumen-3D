#include "Voxelizer.h"

#include "ModelRenderer.h"

namespace Lumen {

	constexpr int CascadeCount = 6;

	static bool GenerateMips = true;


	static float VoxelRanges[CascadeCount] = { 128.0f, 256.0f, 384.0f, 512.0f, 1024.0f, 2048.0f };
	static GLuint VoxelVolumes[CascadeCount]{ 0, 0, 0, 0, 0, 0 };
	static GLuint VoxelVolumesNormals[CascadeCount]{ 0, 0, 0, 0, 0, 0 };

	static glm::vec3 VoxelCenters[6];

	static GLClasses::ComputeShader* ClearShader;

	static GLClasses::Shader* VoxelizeShader;

	static float Align(float value, float size)
	{
		return std::floor(value / size) * size;
	}

	static glm::vec3 SnapPosition(glm::vec3 p, int Cascade) {

		float amts[6] = { 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f };

		p.x = Align(p.x, amts[Cascade]);
		p.y = Align(p.y, amts[Cascade]);
		p.z = Align(p.z, amts[Cascade]);

		return p;
	}


	void Voxelizer::CreateVolumes()
	{
		for (int Volume = 0; Volume < CascadeCount; Volume++) {

			glGenTextures(1, &VoxelVolumes[Volume]);
			glBindTexture(GL_TEXTURE_3D, VoxelVolumes[Volume]);

			if (GenerateMips)
			{
				glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
			}

			else {
				glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

			}
			glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
			glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
			glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
			glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
			glTexImage3D(GL_TEXTURE_3D, 0, GL_RGBA16F, 128, 128, 128, 0, GL_RGBA, GL_FLOAT, nullptr);

			// Generate normal texture (ui16 with cube normal encoding)
			glGenTextures(1, &VoxelVolumesNormals[Volume]);
			glBindTexture(GL_TEXTURE_3D, VoxelVolumesNormals[Volume]);
			glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
			glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
			glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
			glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
			glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
			glTexImage3D(GL_TEXTURE_3D, 0, GL_R16UI, 128, 128, 128, 0, GL_RED_INTEGER, GL_UNSIGNED_SHORT, nullptr);
		}


		ClearShader = new GLClasses::ComputeShader();
		ClearShader->CreateComputeShader("Core/Shaders/ClearVolume.comp");
		ClearShader->Compile();

		VoxelizeShader = new GLClasses::Shader();
		VoxelizeShader->CreateShaderProgramFromFile("Core/Shaders/Voxelizer/VoxelizationVertex.glsl",
													"Core/Shaders/Voxelizer/VoxelizationRadiance.glsl",
													"Core/Shaders/Voxelizer/VoxelizationGeometry.geom");
		VoxelizeShader->CompileShaders();

		for (int x = 0; x < 6; x++) {
			VoxelCenters[x] = glm::vec3(0.0f);
		}
	}

	void Voxelizer::VoxelizeCascade(int Cascade, glm::vec3 Position, const glm::mat4& Projection, const glm::mat4& View, GLuint Shadowmap, const glm::mat4& LVP, const glm::vec3& SunDirection, const std::vector<Entity*>& EntityList)
	{
		if (Cascade >= CascadeCount) {
			throw "wtf.";
		}

		Position = SnapPosition(Position,Cascade);

		VoxelCenters[Cascade] = Position;

		glBindTexture(GL_TEXTURE_3D, VoxelVolumes[Cascade]);
		glDisable(GL_CULL_FACE);
		glDisable(GL_DEPTH_TEST);
		glDisable(GL_BLEND);

		glBindFramebuffer(GL_FRAMEBUFFER, 0);
		glUseProgram(0);

		int GROUP_SIZE = 8;

		ClearShader->Use();
		glBindImageTexture(0, VoxelVolumes[Cascade], 0, GL_TRUE, 0, GL_READ_WRITE, GL_RGBA16F);
		glBindImageTexture(1, VoxelVolumesNormals[Cascade], 0, GL_TRUE, 0, GL_READ_WRITE, GL_R16UI);
		glDispatchCompute(128 / GROUP_SIZE, 128 / GROUP_SIZE, 128 / GROUP_SIZE);
		glMemoryBarrier(GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);

		// Voxelize ->

		VoxelizeShader->Use();

		glBindTexture(GL_TEXTURE_3D, VoxelVolumes[Cascade]);
		//glBindImageTexture(0, VoxelVolumes[Cascade], 0, GL_TRUE, 0, GL_READ_WRITE, GL_RGBA8);
		glBindImageTexture(0, VoxelVolumes[Cascade], 0, GL_TRUE, 0, GL_READ_WRITE, GL_RGBA16F);
		glBindImageTexture(1, VoxelVolumesNormals[Cascade], 0, GL_TRUE, 0, GL_READ_WRITE, GL_R16UI);

		glm::mat4 VP = Projection * View;

		glViewport(0, 0, 128, 128);

		VoxelizeShader->SetMatrix4("u_VP", VP);
		VoxelizeShader->SetVector3f("u_CoverageSize", glm::vec3(VoxelRanges[Cascade]));
		VoxelizeShader->SetFloat("u_CoverageSizeF",(VoxelRanges[Cascade]));
		VoxelizeShader->SetVector3f("u_VoxelGridCenter", Position);
		VoxelizeShader->SetVector3f("u_VoxelGridCenterF", Position);
		VoxelizeShader->SetVector3f("u_SunDir", SunDirection);
		VoxelizeShader->SetInteger("u_VolumeSize", 128);
		VoxelizeShader->SetInteger("u_Shadowmap", 10);
		VoxelizeShader->SetInteger("u_CascadeNumber", Cascade);
		VoxelizeShader->SetMatrix4("u_LightVP", LVP);

		for (int i = 0; i < 6; i++) {
			std::string s = "u_VoxelRanges[" + std::to_string(i) + ("]");
			VoxelizeShader->SetFloat(s.c_str(), (float)VoxelRanges[i]);
		}

		glActiveTexture(GL_TEXTURE10);
		glBindTexture(GL_TEXTURE_2D, Shadowmap);

		for (auto& e : EntityList) {

			RenderEntity(*e, *VoxelizeShader);

		}


		if (GenerateMips)
		{
			// Generate mipmap levels 

			glBindTexture(GL_TEXTURE_3D, VoxelVolumes[Cascade]);
			glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
			glGenerateMipmap(GL_TEXTURE_3D);
			glBindTexture(GL_TEXTURE_3D, 0);

		}
	}

	GLuint* Voxelizer::GetVolumes()
	{
		return VoxelVolumes;
	}

	GLuint* Voxelizer::GetVolumeNormals()
	{
		return VoxelVolumesNormals;
	}

	float* Voxelizer::GetVolumeRanges()
	{
		return VoxelRanges;
	}

	glm::vec3* Voxelizer::GetVolumeCenters()
	{
		return VoxelCenters;
	}

	void Voxelizer::RecompileShaders()
	{
		delete ClearShader;
		delete VoxelizeShader;

		ClearShader = new GLClasses::ComputeShader();
		ClearShader->CreateComputeShader("Core/Shaders/ClearVolume.comp");
		ClearShader->Compile();

		VoxelizeShader = new GLClasses::Shader();
		VoxelizeShader->CreateShaderProgramFromFile("Core/Shaders/Voxelizer/VoxelizationVertex.glsl",
			"Core/Shaders/Voxelizer/VoxelizationRadiance.glsl",
			"Core/Shaders/Voxelizer/VoxelizationGeometry.geom");
		VoxelizeShader->CompileShaders();
	}


}