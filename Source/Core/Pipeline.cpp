#define DUPLICATE(x) (x, x)

#include "Pipeline.h"

#include "FpsCamera.h"
#include "GLClasses/Shader.h"
#include "Object.h"
#include "Entity.h"
#include "ModelFileLoader.h"
#include "ModelRenderer.h"
#include "GLClasses/Fps.h"
#include "GLClasses/Framebuffer.h"
#include "ShaderManager.h"
#include "GLClasses/DepthBuffer.h"
#include "ShadowRenderer.h"
#include "GLClasses/CubeTextureMap.h"

#include "ProbeMap.h"

#include "TAAJitter.h"

#include <string>

Lumen::FPSCamera Camera(90.0f, 800.0f / 600.0f);

static bool vsync = false;
static float SunTick = 50.0f;
static glm::vec3 SunDirection = glm::vec3(0.1f, -1.0f, 0.1f);

// Flags 
static bool TAA = true;
static float SpecularIndirectRes = 0.25f;
static float SpecularIndirectUpsampleRes = 0.5f;


class RayTracerApp : public Lumen::Application
{
public:

	RayTracerApp()
	{
		m_Width = 800;
		m_Height = 600;
	}

	void OnUserCreate(double ts) override
	{
	
	}

	void OnUserUpdate(double ts) override
	{
		glfwSwapInterval((int)vsync);

		GLFWwindow* window = GetWindow();
		float camera_speed = 0.525f;

		if (GetCursorLocked()) {
			if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
				Camera.ChangePosition(Camera.GetFront() * camera_speed);

			if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
				Camera.ChangePosition(-(Camera.GetFront() * camera_speed));

			if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
				Camera.ChangePosition(-(Camera.GetRight() * camera_speed));

			if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
				Camera.ChangePosition(Camera.GetRight() * camera_speed);

			if (glfwGetKey(window, GLFW_KEY_SPACE) == GLFW_PRESS)
				Camera.ChangePosition(Camera.GetUp() * camera_speed);

			if (glfwGetKey(window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS)
				Camera.ChangePosition(-(Camera.GetUp() * camera_speed));

		}
	}

	void OnImguiRender(double ts) override
	{
		ImGui::Text("Position : %f,  %f,  %f", Camera.GetPosition().x, Camera.GetPosition().y, Camera.GetPosition().z);
		ImGui::Text("Front : %f,  %f,  %f", Camera.GetFront().x, Camera.GetFront().y, Camera.GetFront().z);
		ImGui::SliderFloat("Sun Time ", &SunTick, 0.1f, 256.0f);
		ImGui::SliderFloat3("Sun Dir : ", &SunDirection[0], -1.0f, 1.0f);

		ImGui::NewLine();
		ImGui::NewLine();

		ImGui::SliderFloat("Specular Resolution", &SpecularIndirectRes, 0.1f, 1.0f);
		ImGui::SliderFloat("Specular Upsample Resolution", &SpecularIndirectUpsampleRes, 0.1f, 1.0f);

		ImGui::Checkbox("TAA", &TAA);
	}

	void OnEvent(Lumen::Event e) override
	{
		if (e.type == Lumen::EventTypes::MouseMove && GetCursorLocked())
		{
			Camera.UpdateOnMouseMovement(e.mx, e.my);
		}

		if (e.type == Lumen::EventTypes::WindowResize)
		{
			Camera.SetAspect((float)e.wx / (float)e.wy);
		}

		if (e.type == Lumen::EventTypes::KeyPress && e.key == GLFW_KEY_ESCAPE) {
			exit(0);
		}

		if (e.type == Lumen::EventTypes::KeyPress && e.key == GLFW_KEY_F1)
		{
			this->SetCursorLocked(!this->GetCursorLocked());
		}

		if (e.type == Lumen::EventTypes::KeyPress && e.key == GLFW_KEY_F2 && this->GetCurrentFrame() > 5)
		{
			Lumen::ShaderManager::RecompileShaders();
		}

		if (e.type == Lumen::EventTypes::KeyPress && e.key == GLFW_KEY_V && this->GetCurrentFrame() > 5)
		{
			vsync = !vsync;
		}
	}


};

struct CommonUniforms {
	glm::mat4 View, Projection, InvView, InvProjection, PrevProj, PrevView, InvPrevProj, InvPrevView;
	int Frame;
};

void SetCommonUniforms(GLClasses::Shader& shader, CommonUniforms& uniforms) {
	shader.SetFloat("u_Time", glfwGetTime());
	shader.SetInteger("u_Frame", uniforms.Frame);
	shader.SetInteger("u_CurrentFrame", uniforms.Frame);
	shader.SetMatrix4("u_ViewProjection", Camera.GetViewProjection());
	shader.SetMatrix4("u_Projection", uniforms.Projection);
	shader.SetMatrix4("u_View", uniforms.View);
	shader.SetMatrix4("u_InverseProjection", uniforms.InvProjection);
	shader.SetMatrix4("u_InverseView", uniforms.InvView);
	shader.SetMatrix4("u_PrevProjection", uniforms.PrevProj);
	shader.SetMatrix4("u_PrevView", uniforms.PrevView);
	shader.SetMatrix4("u_PrevInverseProjection", uniforms.InvPrevProj);
	shader.SetMatrix4("u_PrevInverseView", uniforms.InvPrevView);
	shader.SetVector3f("u_ViewerPosition", glm::vec3(uniforms.InvView[3]));
	shader.SetVector3f("u_Incident", glm::vec3(uniforms.InvView[3]));
	shader.SetVector3f("u_SunDirection", SunDirection);
	shader.SetFloat("u_zNear", Camera.GetNearPlane());
	shader.SetFloat("u_zFar", Camera.GetFarPlane());
}

void UnbindEverything() {
	glBindFramebuffer(GL_FRAMEBUFFER, 0);
	glUseProgram(0);
}

void RenderEntityList(const std::vector<Lumen::Entity*> EntityList, GLClasses::Shader& shader) {
	for (auto& e : EntityList) {
		Lumen::RenderEntity(*e, shader);
	}
}

void RenderProbe(Lumen::ProbeMap& probe, int face, glm::vec3 center, const std::vector<Lumen::Entity*> EntityList, GLClasses::Shader& shader) {

	if (face >= 6) {
		throw "What.";
	}

	glEnable(GL_CULL_FACE);
	glEnable(GL_DEPTH_TEST);

	probe.CapturePoints[face] = center;

	const glm::mat4 projection_matrix = glm::perspective(glm::radians(90.0f), 1.0f, 0.1f, 800.0f);

	std::array<glm::mat4, 6> view_matrices =
	{
		glm::lookAt(center, center + glm::vec3(1.0f, 0.0f, 0.0f), glm::vec3(0.0f, -1.0f,  0.0f)),
		glm::lookAt(center, center + glm::vec3(-1.0f, 0.0f, 0.0f), glm::vec3(0.0f, -1.0f,  0.0f)),
		glm::lookAt(center, center + glm::vec3(0.0f, 1.0f, 0.0f), glm::vec3(0.0f,  0.0f,  1.0f)),
		glm::lookAt(center, center + glm::vec3(0.0f,-1.0f, 0.0f), glm::vec3(0.0f,  0.0f, -1.0f)),
		glm::lookAt(center, center + glm::vec3(0.0f, 0.0f, 1.0f), glm::vec3(0.0f, -1.0f,  0.0f)),
		glm::lookAt(center, center + glm::vec3(0.0f, 0.0f,-1.0f), glm::vec3(0.0f, -1.0f,  0.0f))
	};

	probe.BindFace(face, false);
	shader.Use();
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	shader.SetVector3f("u_CapturePosition", center);
	shader.SetMatrix4("u_ViewProjection", projection_matrix * view_matrices[face]);
	shader.SetInteger("u_AlbedoMap", 0);
	shader.SetInteger("u_NormalMap", 1);
	shader.SetInteger("u_RoughnessMap", 2);
	shader.SetInteger("u_MetalnessMap", 3);
	shader.SetInteger("u_MetalnessRoughnessMap", 5);
	RenderEntityList(EntityList, shader);

	glDisable(GL_CULL_FACE);
	glDisable(GL_DEPTH_TEST);
}

void RenderProbeAllFaces(Lumen::ProbeMap& probe, const glm::vec3& center, const std::vector<Lumen::Entity*> EntityList, GLClasses::Shader& shader) 
{
	for (int i = 0; i < 6; i++) {
		RenderProbe(probe, i, center, EntityList, shader);
	}

	return;
}



GLClasses::Framebuffer GBuffers[2] = { GLClasses::Framebuffer(16, 16, { {GL_RGB, GL_RGB, GL_UNSIGNED_BYTE, true, true}, {GL_RGB16F, GL_RGB, GL_FLOAT, true, true}, {GL_RGB, GL_RGB, GL_UNSIGNED_BYTE, false, false}, {GL_RGB16F, GL_RGB, GL_FLOAT, true, true} }, false, true), GLClasses::Framebuffer(16, 16, { {GL_RGB, GL_RGB, GL_UNSIGNED_BYTE, true, true}, {GL_RGB16F, GL_RGB, GL_FLOAT, true, true}, {GL_RGB, GL_RGB, GL_UNSIGNED_BYTE, false, false}, {GL_RGB16F, GL_RGB, GL_FLOAT, true, true} }, false, true) };
GLClasses::Framebuffer LightingPass(16, 16, {GL_RGB16F, GL_RGB, GL_FLOAT, true, true}, false, true);
GLClasses::Framebuffer SpecularIndirect(16, 16, {GL_RGB16F, GL_RGB, GL_FLOAT, true, true}, false, true);
GLClasses::Framebuffer SpecularIndirectTemporalBuffers[2] = { GLClasses::Framebuffer(16, 16, {GL_RGB16F, GL_RGB, GL_FLOAT, true, true}, false, true), GLClasses::Framebuffer(16, 16, {GL_RGB16F, GL_RGB, GL_FLOAT, true, true}, false, true) };
GLClasses::Framebuffer TAABuffers[2] = { GLClasses::Framebuffer(16, 16, {GL_RGB16F, GL_RGB, GL_FLOAT, true, true}, false, true), GLClasses::Framebuffer(16, 16, {GL_RGB16F, GL_RGB, GL_FLOAT, true, true}, false, true) };

void Lumen::StartPipeline()
{
	RayTracerApp app;
	app.Initialize();
	app.SetCursorLocked(true);

	// Scene setup 
	Object Sponza;
	FileLoader::LoadModelFile(&Sponza, "Models/sponza-pbr/Sponza.gltf");
	Entity MainModel(&Sponza);
	MainModel.m_Model = glm::scale(glm::mat4(1.0f), glm::vec3(0.2f));
	MainModel.m_Model = glm::translate(MainModel.m_Model, glm::vec3(0.0f));
	std::vector<Entity*> EntityRenderList = { &MainModel };
	auto& EntityList = EntityRenderList;

	// Data object initialization 
	GLClasses::VertexBuffer ScreenQuadVBO;
	GLClasses::VertexArray ScreenQuadVAO;
	GLClasses::DepthBuffer Shadowmap(3584, 3584);
	GLClasses::Texture BlueNoise;
	GLClasses::CubeTextureMap Skymap;

	Skymap.CreateCubeTextureMap(
		{
		"Res/Skymap/right.bmp",
		"Res/Skymap/left.bmp",
		"Res/Skymap/top.bmp",
		"Res/Skymap/bottom.bmp",
		"Res/Skymap/front.bmp",
		"Res/Skymap/back.bmp"
		}, true
	);

	BlueNoise.CreateTexture("Res/blue_noise.png", false, false);

	// Setup screensized quad for rendering
	{
		unsigned long long CurrentFrame = 0;
		float QuadVertices_NDC[] =
		{
			-1.0f,  1.0f,  0.0f, 1.0f, -1.0f, -1.0f,  0.0f, 0.0f,
			 1.0f, -1.0f,  1.0f, 0.0f, -1.0f,  1.0f,  0.0f, 1.0f,
			 1.0f, -1.0f,  1.0f, 0.0f,  1.0f,  1.0f,  1.0f, 1.0f
		};

		ScreenQuadVAO.Bind();
		ScreenQuadVBO.Bind();
		ScreenQuadVBO.BufferData(sizeof(QuadVertices_NDC), QuadVertices_NDC, GL_STATIC_DRAW);
		ScreenQuadVBO.VertexAttribPointer(0, 2, GL_FLOAT, 0, 4 * sizeof(GLfloat), 0);
		ScreenQuadVBO.VertexAttribPointer(1, 2, GL_FLOAT, 0, 4 * sizeof(GLfloat), (void*)(2 * sizeof(GLfloat)));
		ScreenQuadVAO.Unbind();
	}


	// Create Shaders
	ShaderManager::CreateShaders();
	GLClasses::Shader& GBufferShader = ShaderManager::GetShader("GBUFFER");
	GLClasses::Shader& LightingShader = ShaderManager::GetShader("LIGHTING_PASS");
	GLClasses::Shader& FinalShader = ShaderManager::GetShader("FINAL");
	GLClasses::Shader& ProbeForwardShader = ShaderManager::GetShader("PROBE_FORWARD");
	GLClasses::Shader& ProbeSpecularShader = ShaderManager::GetShader("PROBE_SPECULAR");
	GLClasses::Shader& SpecularTemporalShader = ShaderManager::GetShader("SPECULAR_TEMPORAL");
	GLClasses::Shader& TAAShader = ShaderManager::GetShader("TAA");

	// History
	glm::mat4 PreviousView;
	glm::mat4 PreviousProjection;
	glm::mat4 View;
	glm::mat4 Projection;
	glm::mat4 InverseView;
	glm::mat4 InverseProjection;

	// Probe Setup
	ProbeMap PlayerProbe(192);

	// Temporal jitter
	GenerateJitterStuff();

	while (!glfwWindowShouldClose(app.GetWindow()))
	{
		// Matrices 
		PreviousProjection = Camera.GetProjectionMatrix();
		PreviousView = Camera.GetViewMatrix();

		// App update 
		app.OnUpdate();

		// Update current matrices 
		Projection = Camera.GetProjectionMatrix();
		View = Camera.GetViewMatrix();
		InverseProjection = glm::inverse(Camera.GetProjectionMatrix());
		InverseView = glm::inverse(Camera.GetViewMatrix());

		// Common uniform buffer
		CommonUniforms UniformBuffer = { View, Projection, InverseView, InverseProjection, PreviousProjection, PreviousView, glm::inverse(PreviousProjection), glm::inverse(PreviousView), (int)app.GetCurrentFrame() };

		// Resize buffers (to maintain aspect ratio)
		GBuffers[0].SetSize(app.GetWidth(), app.GetHeight());
		GBuffers[1].SetSize(app.GetWidth(), app.GetHeight());
		LightingPass.SetSize(app.GetWidth(), app.GetHeight());
		TAABuffers[0].SetSize(app.GetWidth(), app.GetHeight());
		TAABuffers[1].SetSize(app.GetWidth(), app.GetHeight());
		SpecularIndirect.SetSize(app.GetWidth() * SpecularIndirectRes, app.GetHeight() * SpecularIndirectRes);
		SpecularIndirectTemporalBuffers[0].SetSize(app.GetWidth() * SpecularIndirectUpsampleRes, app.GetHeight() * SpecularIndirectUpsampleRes);
		SpecularIndirectTemporalBuffers[1].SetSize(app.GetWidth() * SpecularIndirectUpsampleRes, app.GetHeight() * SpecularIndirectUpsampleRes);

		// Shadowmap update 
		if (app.GetCurrentFrame() % 8 == 0)
		{
			// Shadow pass 
			RenderShadowMap(Shadowmap, SunDirection, EntityRenderList, Camera.GetViewProjection());
		}

		// Probe update (single slice)
		RenderProbe(PlayerProbe, app.GetCurrentFrame() % 6, Camera.GetPosition(), EntityList, ProbeForwardShader);

		// Ping pong framebuffers
		bool FrameCheckerStep = app.GetCurrentFrame() % 2 == 0;
		
		// Gbuffer 
		auto& GBuffer = FrameCheckerStep ? GBuffers[0] : GBuffers[1];
		auto& PrevGBuffer = FrameCheckerStep ? GBuffers[1] : GBuffers[0];

		// Specular 
		auto& SpecularTemporal = FrameCheckerStep ? SpecularIndirectTemporalBuffers[0] : SpecularIndirectTemporalBuffers[1];
		auto& PrevSpecularTemporal = FrameCheckerStep ? SpecularIndirectTemporalBuffers[1] : SpecularIndirectTemporalBuffers[0];

		// TAA
		auto& TAATemporal = FrameCheckerStep ? TAABuffers[0] : TAABuffers[1];
		auto& PrevTAATemporal = FrameCheckerStep ? TAABuffers[1] : TAABuffers[0];


		// Render GBuffer
		glEnable(GL_CULL_FACE);
		glEnable(GL_DEPTH_TEST);

		GBuffer.Bind();
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		GBufferShader.Use();
		GBufferShader.SetMatrix4("u_ViewProjection", Camera.GetViewProjection());
		GBufferShader.SetInteger("u_AlbedoMap", 0);
		GBufferShader.SetInteger("u_NormalMap", 1);
		GBufferShader.SetInteger("u_RoughnessMap", 2);
		GBufferShader.SetInteger("u_MetalnessMap", 3);
		GBufferShader.SetInteger("u_MetalnessRoughnessMap", 5);
		GBufferShader.SetMatrix4("u_JitterMatrix", GetTAAJitterMatrix(app.GetCurrentFrame(), GBuffer.GetDimensions()));

		RenderEntityList(EntityRenderList, GBufferShader);
		UnbindEverything();

		// Post processing passes here : 
		glDisable(GL_CULL_FACE);
		glDisable(GL_DEPTH_TEST);

		SpecularIndirect.Bind();
		ProbeSpecularShader.Use();
		ProbeSpecularShader.SetVector3f("u_Incident", Camera.GetPosition());
		ProbeSpecularShader.SetInteger("u_Depth", 0);
		ProbeSpecularShader.SetInteger("u_Normals", 1);
		ProbeSpecularShader.SetInteger("u_PBR", 2);
		ProbeSpecularShader.SetInteger("u_ProbeAlbedo", 4);
		ProbeSpecularShader.SetInteger("u_ProbeDepth", 5);
		ProbeSpecularShader.SetInteger("u_ProbeNormals", 6);
		ProbeSpecularShader.SetInteger("u_EnvironmentMap", 7);
		ProbeSpecularShader.SetInteger("u_Shadowmap", 8);
		ProbeSpecularShader.SetInteger("u_Frame", app.GetCurrentFrame());
		ProbeSpecularShader.SetVector3f("u_SunDirection", SunDirection);
		ProbeSpecularShader.SetMatrix4("u_SunShadowMatrix", GetLightViewProjection(SunDirection));
		ProbeSpecularShader.SetVector2f("u_Jitter", GetTAAJitterSecondary(app.GetCurrentFrame()));
		ProbeSpecularShader.SetVector2f("u_Dimensions", SpecularIndirect.GetDimensions());
		SetCommonUniforms(ProbeSpecularShader, UniformBuffer);

		for (int i = 0; i < 6; i++) {
			std::string name = "u_ProbeCapturePoints[" + std::to_string(i) + "]";
			ProbeSpecularShader.SetVector3f(name.c_str(), PlayerProbe.CapturePoints[i]);
		}

		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, GBuffer.GetDepthBuffer());

		glActiveTexture(GL_TEXTURE1);
		glBindTexture(GL_TEXTURE_2D, GBuffer.GetTexture(1));
		//glBindTexture(GL_TEXTURE_2D, GBuffer.GetTexture(3));

		glActiveTexture(GL_TEXTURE2);
		glBindTexture(GL_TEXTURE_2D, GBuffer.GetTexture(2));

		glActiveTexture(GL_TEXTURE4);
		glBindTexture(GL_TEXTURE_CUBE_MAP, PlayerProbe.m_CubemapTexture);

		glActiveTexture(GL_TEXTURE5);
		glBindTexture(GL_TEXTURE_CUBE_MAP, PlayerProbe.m_DepthCubemap);

		glActiveTexture(GL_TEXTURE6);
		glBindTexture(GL_TEXTURE_CUBE_MAP, PlayerProbe.NormalPBRPackedCubemap);

		glActiveTexture(GL_TEXTURE7);
		glBindTexture(GL_TEXTURE_CUBE_MAP, Skymap.GetID());

		glActiveTexture(GL_TEXTURE8);
		glBindTexture(GL_TEXTURE_2D, Shadowmap.GetDepthTexture());


		ScreenQuadVAO.Bind();
		glDrawArrays(GL_TRIANGLES, 0, 6);
		ScreenQuadVAO.Unbind();

		// Temporal resolve :

		SpecularTemporalShader.Use();
		SpecularTemporal.Bind();

		SpecularTemporalShader.SetInteger("u_Specular", 0);
		SpecularTemporalShader.SetInteger("u_HistorySpecular", 1);
		SpecularTemporalShader.SetInteger("u_Depth", 2);
		SpecularTemporalShader.SetInteger("u_PreviousDepth", 3);
		SetCommonUniforms(SpecularTemporalShader, UniformBuffer);

		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, SpecularIndirect.GetTexture());

		glActiveTexture(GL_TEXTURE1);
		glBindTexture(GL_TEXTURE_2D, PrevSpecularTemporal.GetTexture());

		glActiveTexture(GL_TEXTURE2);
		glBindTexture(GL_TEXTURE_2D, GBuffer.GetDepthBuffer());

		glActiveTexture(GL_TEXTURE3);
		glBindTexture(GL_TEXTURE_2D, PrevGBuffer.GetDepthBuffer());

		ScreenQuadVAO.Bind();
		glDrawArrays(GL_TRIANGLES, 0, 6);
		ScreenQuadVAO.Unbind();


		// Lighting pass : 

		LightingShader.Use();
		LightingPass.Bind();

		LightingShader.SetInteger("u_AlbedoTexture", 0);
		LightingShader.SetInteger("u_NormalTexture", 1);
		LightingShader.SetInteger("u_PBRTexture", 2);
		LightingShader.SetInteger("u_DepthTexture", 3);
		LightingShader.SetInteger("u_ShadowTexture", 4);
		LightingShader.SetInteger("u_BlueNoise", 5);
		LightingShader.SetInteger("u_Skymap", 6);
		LightingShader.SetInteger("u_Probe", 7);
		LightingShader.SetInteger("u_ResolvedSpecular", 8);

		LightingShader.SetMatrix4("u_LightVP", GetLightViewProjection(SunDirection));
		LightingShader.SetVector2f("u_Dims", glm::vec2(app.GetWidth(), app.GetHeight()));

		LightingShader.SetVector3f("u_LightDirection", SunDirection);
		LightingShader.SetVector3f("u_ViewerPosition", Camera.GetPosition());

		SetCommonUniforms(LightingShader, UniformBuffer);
		
		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, GBuffer.GetTexture(0));

		glActiveTexture(GL_TEXTURE1);
		glBindTexture(GL_TEXTURE_2D, GBuffer.GetTexture(1));

		glActiveTexture(GL_TEXTURE2);
		glBindTexture(GL_TEXTURE_2D, GBuffer.GetTexture(2));

		glActiveTexture(GL_TEXTURE3);
		glBindTexture(GL_TEXTURE_2D, GBuffer.GetDepthBuffer());

		glActiveTexture(GL_TEXTURE4);
		glBindTexture(GL_TEXTURE_2D, Shadowmap.GetDepthTexture());

		glActiveTexture(GL_TEXTURE5);
		glBindTexture(GL_TEXTURE_2D, BlueNoise.GetTextureID());

		glActiveTexture(GL_TEXTURE6);
		glBindTexture(GL_TEXTURE_CUBE_MAP, Skymap.GetID());

		glActiveTexture(GL_TEXTURE7);
		glBindTexture(GL_TEXTURE_CUBE_MAP, PlayerProbe.m_CubemapTexture);
		
		glActiveTexture(GL_TEXTURE8);
		glBindTexture(GL_TEXTURE_2D, SpecularTemporal.GetTexture());

		ScreenQuadVAO.Bind();
		glDrawArrays(GL_TRIANGLES, 0, 6);
		ScreenQuadVAO.Unbind();


		// TAA 

		TAAShader.Use();
		TAATemporal.Bind();

		TAAShader.SetInteger("u_CurrentColorTexture", 0);
		TAAShader.SetInteger("u_PreviousColorTexture", 1);
		TAAShader.SetInteger("u_DepthTexture", 2);
		TAAShader.SetInteger("u_PreviousDepthTexture", 3);
		TAAShader.SetBool("u_Enabled", TAA);
		SetCommonUniforms(TAAShader, UniformBuffer);

		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, LightingPass.GetTexture());

		glActiveTexture(GL_TEXTURE1);
		glBindTexture(GL_TEXTURE_2D, PrevTAATemporal.GetTexture());

		glActiveTexture(GL_TEXTURE2);
		glBindTexture(GL_TEXTURE_2D, GBuffer.GetDepthBuffer());

		glActiveTexture(GL_TEXTURE3);
		glBindTexture(GL_TEXTURE_2D, PrevGBuffer.GetDepthBuffer());

		ScreenQuadVAO.Bind();
		glDrawArrays(GL_TRIANGLES, 0, 6);
		ScreenQuadVAO.Unbind();



		// Tonemapper 

		glBindFramebuffer(GL_FRAMEBUFFER, 0);
		glViewport(0, 0, app.GetWidth(), app.GetHeight());

		FinalShader.Use();
		FinalShader.SetInteger("u_MainTexture", 0);

		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, TAATemporal.GetTexture(0));

		ScreenQuadVAO.Bind();
		glDrawArrays(GL_TRIANGLES, 0, 6);
		ScreenQuadVAO.Unbind();

		// Finish : 
		glFinish();
		app.FinishFrame();
		GLClasses::DisplayFrameRate(app.GetWindow(), "Lumen ");

	}
}
