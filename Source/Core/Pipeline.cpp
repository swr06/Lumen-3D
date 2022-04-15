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

#include "BVH/BVHConstructor.h"

#include "ProbeMap.h"

#include "TAAJitter.h"

#include <string>

#include "Voxelizer.h"

#include "BloomRenderer.h"
#include "BloomFBO.h"

static float RENDER_SCALE = 1.0f;
static int FXAA_PASSES = 2;
static float CAS_Amount = 0.375f;
static bool TAAU = false;
static float TAAUConfidenceExponent = 0.75f;

// Camera
Lumen::FPSCamera Camera(90.0f, 800.0f / 600.0f);

static uint32_t ProbePolygons;
static uint32_t MainRenderPolygons;

// Vertical Sync (FPS Cap)
static bool VSync = false;

// Sun Settings 
static float SunTick = 50.0f;
static glm::vec3 SunDirection = glm::vec3(0.1f, -1.0f, 0.1f);

// Flags 
static bool TAA = true;

static const float IndirectRes = 0.5f;
static bool SVGF = true;

// Specular settings
static bool SSCT = false;
static bool RoughSpecular = true;
static bool CheckerboardIndirect = true;

// AO 
const float ScreenspaceOcclusionRes = 0.5f;
static bool DoScreenspaceAO = true;
static bool DoScreenspaceShadow = false;
static float ScreenspaceAOStrength = 0.75f;
static bool ScreenspaceOcclusionCheckerboard = true;

// Volumetrics
const float VolumetricsResolution = 0.5f;
static bool VolumetricLighting = false;
static bool VolumetricsFilter = true;
static bool VolumetricsTemporal = true;
static float SunVolumetricsStrength = 1.0f;

// Bloom
static bool DoBloom = true;

// Update rates
static int ProbeUpdateRate = 1;
static bool FullyDynamicVoxelization = false;

// Spatial
static bool SpatialFiltering = true;

double RoundToNearest(double n, double x) {
	return round(n / x) * x;
}

static float Align(float value, float size)
{
	return std::floor(value / size) * size;
}

static glm::vec3 SnapPosition(glm::vec3 p, float amt) {

	p.x = Align(p.x, amt);
	p.y = Align(p.y, amt);
	p.z = Align(p.z, amt);

	return p;
}

static glm::vec3 SnapPosition(glm::vec3 p, float ax, float ay, float az) {

	p.x = Align(p.x, ax);
	p.y = Align(p.y, ay);
	p.z = Align(p.z, az);

	return p;
}

// Application 
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
		glfwSwapInterval((int)VSync);

		GLFWwindow* window = GetWindow();
		float camera_speed = 0.525f * 3.0f;

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
		ImGui::Text("-- Info --");
		ImGui::Text("Position : %f,  %f,  %f", Camera.GetPosition().x, Camera.GetPosition().y, Camera.GetPosition().z);
		ImGui::Text("Front : %f,  %f,  %f", Camera.GetFront().x, Camera.GetFront().y, Camera.GetFront().z);
		ImGui::Text("VSync : %d", VSync);
		ImGui::Text("Main Polygons Rendered : %d", MainRenderPolygons);
		ImGui::Text("Probe Polygons Rendered : %d", ProbePolygons);

		ImGui::NewLine();
		ImGui::NewLine();

		ImGui::SliderFloat3("Sun Direction (X, Y, Z) : ", &SunDirection[0], -1.0f, 1.0f);

		ImGui::NewLine();
		
		ImGui::SliderInt("Probe Update Rate", &ProbeUpdateRate, 1, 6);
		ImGui::Checkbox("Fully Dynamic Voxelization?", &FullyDynamicVoxelization);

		ImGui::NewLine();

		ImGui::Checkbox("Spatial Filtering?", &SpatialFiltering);

		ImGui::NewLine();

		ImGui::NewLine();
		
		ImGui::Text("Indirect Resolution : %f on each axis", IndirectRes);

		ImGui::NewLine();

		ImGui::Text("Specular Upsample Resolution : %f on each axis", IndirectRes);
		ImGui::Checkbox("Rough Specular?", &RoughSpecular);

		ImGui::NewLine();

		ImGui::Checkbox("Indirect Checkerboarding? (Resolution effectively halved)", &CheckerboardIndirect);
		ImGui::Checkbox("SVGF?", &SVGF);
		ImGui::NewLine();

		ImGui::Text("RTAO/Screenspace shadows resolve resolution : %f on each axis", ScreenspaceOcclusionRes);
		ImGui::Checkbox("Do Screenspace RTAO?", &DoScreenspaceAO);
		ImGui::Checkbox("Do Screenspace Direct Shadows?", &DoScreenspaceShadow);
		ImGui::Checkbox("Checkerboard render?", &ScreenspaceOcclusionCheckerboard);

		if (DoScreenspaceAO)
			ImGui::SliderFloat("Screenspace RTAO Strength", &ScreenspaceAOStrength, 0.1f, 2.0f);

		ImGui::NewLine();
		ImGui::NewLine();
		ImGui::Checkbox("Volumetrics?", &VolumetricLighting);
		ImGui::Text("Volumetrics Raymarch Resolution : %f", VolumetricsResolution);
		ImGui::Checkbox("Temporal Filter Volumetrics?", &VolumetricsTemporal);
		ImGui::Checkbox("Filter Volumetrics?", &VolumetricsFilter);
		ImGui::SliderFloat("Volumetrics Strength", &SunVolumetricsStrength, 0.1f, 4.0f);

		ImGui::NewLine();
		ImGui::NewLine();

		ImGui::Checkbox("Bloom?", &DoBloom);

		ImGui::NewLine();
		ImGui::NewLine();

		ImGui::Text("Antialiasing : ");
		ImGui::NewLine();
		ImGui::SliderFloat("Global Render Resolution", &RENDER_SCALE, 0.25f, 1.0f);
		ImGui::SliderFloat("CAS Amount", &CAS_Amount, 0.0f, 1.0f);
		ImGui::SliderInt("FXAA Passes", &FXAA_PASSES, 0, 4);
		ImGui::Checkbox("TAA", &TAA);

		ImGui::Checkbox("TAA-UPSCALING?", &TAAU);

		if (TAAU) {
			ImGui::SliderFloat("TAA-U Confidence Exponent", &TAAUConfidenceExponent, 0.2f, 8.0f);
		}
	}

	void OnEvent(Lumen::Event e) override
	{
		if (e.type == Lumen::EventTypes::MouseScroll)
		{
			float Sign = e.msy < 0.0f ? 1.0f : -1.0f;
			Camera.SetFov(Camera.GetFov() + 2.0f * Sign);
			Camera.SetFov(glm::clamp(Camera.GetFov(), 1.0f, 89.0f));
		}

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
			Lumen::Voxelizer::RecompileShaders();
		}

		if (e.type == Lumen::EventTypes::KeyPress && e.key == GLFW_KEY_V && this->GetCurrentFrame() > 5)
		{
			VSync = !VSync;
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
	shader.SetMatrix4("u_InversePrevProjection", uniforms.InvPrevProj);
	shader.SetMatrix4("u_InversePrevView", uniforms.InvPrevView);
	shader.SetVector3f("u_ViewerPosition", glm::vec3(uniforms.InvView[3]));
	shader.SetVector3f("u_Incident", glm::vec3(uniforms.InvView[3]));
	shader.SetVector3f("u_SunDirection", SunDirection);
	shader.SetFloat("u_zNear", Camera.GetNearPlane());
	shader.SetFloat("u_zFar", Camera.GetFarPlane());
	shader.SetVector2f("u_HaltonJitter", Lumen::GetTAAJitterSecondary(uniforms.Frame));
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

void RenderProbe(Lumen::ProbeMap& probe, int face, glm::vec3 center, const std::vector<Lumen::Entity*> EntityList, GLClasses::Shader& shader, int env, GLClasses::VertexArray& vao) {

	if (face >= 6) {
		throw "What.";
	}

	//center = SnapPosition(center, 3.0f, 4.0f, 4.0f);
	center = SnapPosition(center, 0.5f);

	glEnable(GL_CULL_FACE);
	glEnable(GL_DEPTH_TEST);

	probe.CapturePoints[face] = center;

	const glm::mat4 projection_matrix = glm::perspective(glm::radians(90.0f), 1.0f, 0.1f, 800.0f);
	const glm::mat4 inverse_projection = glm::inverse(projection_matrix);

	auto& ProbeSkyShader = Lumen::ShaderManager::GetShader("PROBE_SKY");

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
	glClearColor(0., 0., 0., 0.);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	shader.SetVector3f("u_CapturePosition", center);
	shader.SetMatrix4("u_ViewProjection", projection_matrix * view_matrices[face]);
	shader.SetInteger("u_AlbedoMap", 0);
	shader.SetInteger("u_NormalMap", 1);
	shader.SetInteger("u_RoughnessMap", 2);
	shader.SetInteger("u_MetalnessMap", 3);
	shader.SetInteger("u_MetalnessRoughnessMap", 5);
	RenderEntityList(EntityList, shader);

	//glDisable(GL_CULL_FACE);
	//glDisable(GL_DEPTH_TEST);
	//
	//ProbeSkyShader.Use();
	//ProbeSkyShader.SetInteger("u_EnvironmentMap", 0);
	//ProbeSkyShader.SetInteger("u_Mask", 1);
	//ProbeSkyShader.SetMatrix4("u_InverseProjection", inverse_projection);
	//ProbeSkyShader.SetMatrix4("u_InverseView", glm::inverse(view_matrices[face]));
	//
	//glActiveTexture(GL_TEXTURE0);
	//glBindTexture(GL_TEXTURE_CUBE_MAP, env);
	//
	//glActiveTexture(GL_TEXTURE1);
	//glBindTexture(GL_TEXTURE_CUBE_MAP, probe.m_DepthCubemap);
	//
	//vao.Bind();
	//glDrawArrays(GL_TRIANGLES, 0, 6);
	//vao.Unbind();
}

void RenderProbeAllFaces(Lumen::ProbeMap& probe, const glm::vec3& center, const std::vector<Lumen::Entity*> EntityList, GLClasses::Shader& shader, int env, GLClasses::VertexArray& vao)
{
	for (int i = 0; i < 6; i++) {
		RenderProbe(probe, i, center, EntityList, shader, env, vao);
	}

	return;
}


// Geometry buffer (For deferred shading)
GLClasses::Framebuffer GBuffers[2] = { GLClasses::Framebuffer(16, 16, { {GL_RGB16F, GL_RGB, GL_FLOAT, true, true}, {GL_RGB16F, GL_RGB, GL_FLOAT, true, true}, {GL_RGBA, GL_RGBA, GL_UNSIGNED_BYTE, false, false}, {GL_RGB16F, GL_RGB, GL_FLOAT, true, true} }, false, true), GLClasses::Framebuffer(16, 16, { {GL_RGB16F, GL_RGB, GL_FLOAT, true, true}, {GL_RGB16F, GL_RGB, GL_FLOAT, true, true}, {GL_RGBA, GL_RGBA, GL_UNSIGNED_BYTE, false, false}, {GL_RGB16F, GL_RGB, GL_FLOAT, true, true} }, false, true) };

// Lighting 
GLClasses::Framebuffer LightingPass(16, 16, {GL_RGB16F, GL_RGB, GL_FLOAT, true, true}, false, false);


// Downsampled buffers
GLClasses::Framebuffer DownsampledGBuffer(16, 16, { {GL_R32F, GL_RED, GL_FLOAT, true, true}, {GL_RGBA16F, GL_RGBA, GL_FLOAT, true, true}, {GL_RGBA16F, GL_RGBA, GL_FLOAT, true, true} }, false, false);

// AO, Screenspace shadows 
GLClasses::Framebuffer ScreenspaceOcclusion(16, 16, { { GL_RED, GL_RED, GL_UNSIGNED_BYTE, true, true } , { GL_RED, GL_RED, GL_UNSIGNED_BYTE, true, true }}, false, false);
GLClasses::Framebuffer ScreenspaceOcclusionCheckerboardConstruct(16, 16, { { GL_RED, GL_RED, GL_UNSIGNED_BYTE, true, true } , { GL_RED, GL_RED, GL_UNSIGNED_BYTE, true, true }}, false, false);
GLClasses::Framebuffer ScreenspaceOcclusionTemporalBuffers[2] = { GLClasses::Framebuffer(16, 16, { { GL_RED, GL_RED, GL_UNSIGNED_BYTE, true, true } , { GL_RED, GL_RED, GL_UNSIGNED_BYTE, true, true }}, false, false), GLClasses::Framebuffer(16, 16, { { GL_RED, GL_RED, GL_UNSIGNED_BYTE, true, true } , { GL_RED, GL_RED, GL_UNSIGNED_BYTE, true, true } }, false, false) };

// Specular Indirect 
GLClasses::Framebuffer SpecularIndirectBuffers[2]{ GLClasses::Framebuffer(16, 16, { {GL_RGBA16F, GL_RGBA, GL_FLOAT, true, true}, { GL_RED, GL_RED, GL_UNSIGNED_BYTE, true, true }  }, false, false), GLClasses::Framebuffer(16, 16, { {GL_RGBA16F, GL_RGBA, GL_FLOAT, true, true}, { GL_RED, GL_RED, GL_UNSIGNED_BYTE, true, true } }, false, false) };
GLClasses::Framebuffer SpecularIndirectTemporalBuffers[2] = { GLClasses::Framebuffer(16, 16, {{GL_RGBA16F, GL_RGBA, GL_FLOAT, true, true}, {GL_R16F, GL_RED, GL_FLOAT, false, false} }, false, false), GLClasses::Framebuffer(16, 16, {{GL_RGBA16F, GL_RGBA, GL_FLOAT, true, true},{GL_R16F, GL_RED, GL_FLOAT, false, false}}, false, false) };
GLClasses::Framebuffer SpecularIndirectConeTraceInput(16, 16, { GL_RGB16F, GL_RGB, GL_FLOAT, true, true }, false, false);
GLClasses::Framebuffer SpecularIndirectConeTraceInputAlternate(16, 16, { GL_RGB16F, GL_RGB, GL_FLOAT, true, true }, false, false);
GLClasses::Framebuffer SpecularIndirectConeTraced(16, 16, { GL_RGB16F, GL_RGB, GL_FLOAT, true, true }, false, false);

// Indirect diffuse 
GLClasses::Framebuffer DiffuseIndirectTrace(16, 16, { GL_RGBA16F, GL_RGBA, GL_FLOAT, true, true }, false, false);
GLClasses::Framebuffer SVGFVarianceResolve(16, 16, { { GL_RGBA16F, GL_RGBA, GL_FLOAT, true, true }, {GL_R16F, GL_RED, GL_FLOAT, true, true} }, false, false);

// 2 attachments : raw temporal filtered, moments/utility
GLClasses::Framebuffer DiffuseIndirectTemporalBuffers[2] = { GLClasses::Framebuffer(16, 16, {{GL_RGBA16F, GL_RGBA, GL_FLOAT, true, true},{GL_RGBA16F, GL_RGBA, GL_FLOAT, false, false} }, false, false), GLClasses::Framebuffer(16, 16, { {GL_RGBA16F, GL_RGBA, GL_FLOAT, true, true},{GL_RGBA16F, GL_RGBA, GL_FLOAT, false, false}}, false, false) };

// Common indirect buffers 
GLClasses::Framebuffer IndirectCheckerUpscaled(16, 16, { {GL_RGBA16F, GL_RGBA, GL_FLOAT, true, true}, {GL_RGBA16F, GL_RGBA, GL_FLOAT, true, true} }, false, false);

// Volumetrics
GLClasses::Framebuffer VolumetricsBuffers[2] = { GLClasses::Framebuffer(16, 16, { {GL_RGBA16F, GL_RGBA, GL_FLOAT, true, true} }, false, false), GLClasses::Framebuffer(16, 16, { {GL_RGBA16F, GL_RGBA, GL_FLOAT, true, true} }, false, false) };
GLClasses::Framebuffer VolumetricsFiltered(16, 16, { {GL_RGBA16F, GL_RGBA, GL_FLOAT, true, true} }, false, false);

// TAA Buffers
GLClasses::Framebuffer TAABuffers[2] = { GLClasses::Framebuffer(16, 16, {GL_RGBA16F, GL_RGBA, GL_FLOAT, true, true}, false, false), GLClasses::Framebuffer(16, 16, {GL_RGBA16F, GL_RGBA, GL_FLOAT, true, true}, false, false) };
GLClasses::Framebuffer FXAABuffers[2] = { GLClasses::Framebuffer(16, 16, {GL_RGB16F, GL_RGB, GL_FLOAT, true, true}, false, false), GLClasses::Framebuffer(16, 16, {GL_RGB16F, GL_RGB, GL_FLOAT, true, true}, false, false) };
GLClasses::Framebuffer TonemappedFBO(16, 16, { {GL_RGBA16F, GL_RGBA, GL_FLOAT, true, true} }, false, false);

GLClasses::Framebuffer PostProcessCombined(16, 16, { GL_RGB16F, GL_RGB, GL_FLOAT, true, true }, false, false);

// Spatial Buffers 
GLClasses::Framebuffer SpatialFilterBuffers[2] = { GLClasses::Framebuffer(16, 16, {{ GL_RGBA16F, GL_RGBA, GL_FLOAT, true, true },{ GL_RGBA16F, GL_RGBA, GL_FLOAT, true, true }, {GL_R16F, GL_RED, GL_FLOAT, true, true}}, false, false), GLClasses::Framebuffer(16, 16, {{ GL_RGBA16F, GL_RGBA, GL_FLOAT, true, true },{ GL_RGBA16F, GL_RGBA, GL_FLOAT, true, true }, {GL_R16F, GL_RED, GL_FLOAT, true, true}}, false, false) };


int GetUpdateCascade(int Frame) {

	int Rand = rand() % 6000;
	int Rand2 = rand() % 1000;
	int Rand3 = rand() % 1000;

	if (Rand < 3550) {
		if (Rand2 < 500) {
			return 0;
		}

		if (Rand2 < 920) {

			if (Rand3 < 950) {
				return 1;
			}

			return 2;
		}

		return 2;
	}

	if (Rand < 4300) { return 2; }
	if (Rand < 5200) { return 3; }
	if (Rand < 5600) { return 4; }

	return 5;
}




// Main pipeline
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

	const glm::mat4 ZOrientMatrix = glm::mat4(glm::vec4(1.0f, 0.0f, 0.0f, 0.0f), glm::vec4(0.0f, 0.0f, 1.0f, 0.0f), glm::vec4(0.0f, 1.0f, 0.0f, 0.0f), glm::vec4(1.0f));

	Object Cube;
	//FileLoader::LoadModelFile(&SecondaryLargeModel, "Models/Lucy/LucyModel.obj");
	FileLoader::LoadModelFile(&Cube, "Models/cube/Cube.gltf");

	Object RedCube;
	FileLoader::LoadModelFile(&RedCube, "Models/redcube/Cube.gltf");

	Object BlueCube;
	FileLoader::LoadModelFile(&BlueCube, "Models/bluecube/Cube.gltf");

	Object SecondaryModel;
	FileLoader::LoadModelFile(&SecondaryModel, "Models/dragon/dragon.obj");

	Object Suzanne;
	//FileLoader::LoadModelFile(&Suzanne, "Models/suzanne/Suzanne.gltf");
	FileLoader::LoadModelFile(&Suzanne, "Models/suzanneglass/Suzanne.gltf");
	
	Entity SecondaryEntity(&Cube);
	SecondaryEntity.m_Model = glm::scale(glm::mat4(1.0f), glm::vec3(16.0f));
	SecondaryEntity.m_Model = glm::translate(SecondaryEntity.m_Model, glm::vec3(0.0f, 2.0f, 0.0f));

	Entity RedCubeEntity(&RedCube);
	RedCubeEntity.m_Model = glm::scale(glm::mat4(1.0f), glm::vec3(12.0f));
	RedCubeEntity.m_Model = glm::translate(RedCubeEntity.m_Model, glm::vec3(230.0f, 97.0f, 10.0f) / 12.0f);

	Entity BlueCubeEntity(&BlueCube);
	BlueCubeEntity.m_Model = glm::scale(glm::mat4(1.0f), glm::vec3(12.0f));
	BlueCubeEntity.m_Model = glm::translate(BlueCubeEntity.m_Model, glm::vec3(230.0f, 97.0f, -20.0f) / 12.0f);

	Entity RedCubeEntity2(&RedCube);
	RedCubeEntity2.m_Model = glm::scale(glm::mat4(1.0f), glm::vec3(12.0f));
	RedCubeEntity2.m_Model = glm::translate(RedCubeEntity2.m_Model, glm::vec3(-230.0f, 97.0f, 10.0f) / 12.0f);

	Entity BlueCubeEntity2(&BlueCube);
	BlueCubeEntity2.m_Model = glm::scale(glm::mat4(1.0f), glm::vec3(12.0f));
	BlueCubeEntity2.m_Model = glm::translate(BlueCubeEntity2.m_Model, glm::vec3(-230.0f, 97.0f, -20.0f) / 12.0f);



	float EntityScale = 3.5f;
	
	Entity SecondaryEntity0(&SecondaryModel);
	SecondaryEntity0.m_Model = glm::scale(glm::mat4(1.0f), glm::vec3(EntityScale));
	SecondaryEntity0.m_Model = glm::translate(SecondaryEntity0.m_Model, glm::vec3(0.0f, 2.0f, -90.0f) * (1.0f / EntityScale));
	SecondaryEntity0.m_EmissiveAmount = 10.0f;
	
	Entity SecondaryEntity1(&SecondaryModel);
	SecondaryEntity1.m_Model = glm::scale(glm::mat4(1.0f), glm::vec3(EntityScale));
	SecondaryEntity1.m_Model = glm::translate(SecondaryEntity1.m_Model, glm::vec3(0.0f, 2.0f, 90.0f) * (1.0f / EntityScale));
	SecondaryEntity1.m_EmissiveAmount = 10.0f;
	
	Entity SecondaryEntity2(&SecondaryModel);
	SecondaryEntity2.m_Model = glm::scale(glm::mat4(1.0f), glm::vec3(EntityScale));
	SecondaryEntity2.m_Model = glm::translate(SecondaryEntity2.m_Model, glm::vec3(220.0f, 2.0f, 0.0f) * (1.0f / EntityScale));
	SecondaryEntity2.m_EmissiveAmount = 10.0f;
	
	Entity SecondaryEntity3(&SecondaryModel);
	SecondaryEntity3.m_Model = glm::scale(glm::mat4(1.0f), glm::vec3(EntityScale));
	SecondaryEntity3.m_Model = glm::translate(SecondaryEntity3.m_Model, glm::vec3(-220.0f, 2.0f, 0.0f) * (1.0f / EntityScale));
	SecondaryEntity3.m_EmissiveAmount = 10.0f;

	float EntityScale2 = 10.0f;

	Entity SecondaryEntity4(&Suzanne);
	SecondaryEntity4.m_Model = glm::scale(glm::mat4(1.0f), glm::vec3(EntityScale2));
	SecondaryEntity4.m_Model = glm::translate(SecondaryEntity4.m_Model, glm::vec3(50.0f, 24.0f, 0.0f) * (1.0f / (EntityScale2)));
	
	Entity SecondaryEntity5(&Suzanne);
	SecondaryEntity5.m_Model = glm::scale(glm::mat4(1.0f), glm::vec3(EntityScale2));
	SecondaryEntity5.m_Model = glm::translate(SecondaryEntity5.m_Model, glm::vec3(-50, 24.0f, 0.0f) * (1.0f / (EntityScale2)));

	// Construct BVH

	//BVH::Node* SponzaRootBVHNode = BVH::BuildBVH(Sponza);

	// Entity list
	std::vector<Entity*> EntityRenderList = { &MainModel, &SecondaryEntity, &SecondaryEntity0, &SecondaryEntity1, &SecondaryEntity2, &SecondaryEntity3, &RedCubeEntity, &BlueCubeEntity, &RedCubeEntity2, &BlueCubeEntity2, &SecondaryEntity4, &SecondaryEntity5 };
	auto& EntityList = EntityRenderList;

	// Clear CPU side vertex/index data (After bvh construction ofc.) 
	std::vector<Object*> ObjectList = { &Sponza, &Cube, &SecondaryModel };
	
	for (auto& e : ObjectList) {
		e->ClearCPUSideData();
	}
	

	// Data object initialization 
	GLClasses::VertexBuffer ScreenQuadVBO;
	GLClasses::VertexArray ScreenQuadVAO;
	GLClasses::DepthBuffer Shadowmap(4096, 4096);
	GLClasses::Texture BlueNoise;
	GLClasses::CubeTextureMap Skymap;


	glm::vec3 PreviousSunDirection = SunDirection;


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
	GLClasses::Shader& SpecularTraceShader = ShaderManager::GetShader("SPECULAR_TRACE");
	GLClasses::Shader& ScreenspaceConetraceShader = ShaderManager::GetShader("SPECULAR_CONE_TRACE");
	GLClasses::Shader& SpecularTemporalShader = ShaderManager::GetShader("SPECULAR_TEMPORAL");
	GLClasses::Shader& IndirectCheckerboarder = ShaderManager::GetShader("INDIRECT_CHECKER");

	GLClasses::Shader& DiffuseVXTrace = ShaderManager::GetShader("VOXEL_DIFFUSE_GI");

	GLClasses::Shader& TAAShader = ShaderManager::GetShader("TAA");
	GLClasses::Shader& SSOcclusionTraceShader = ShaderManager::GetShader("SCREENSPACE_OCCLUSION_RT");
	GLClasses::Shader& SSOcclusionTemporalShader = ShaderManager::GetShader("SCREENSPACE_OCCLUSION_TEMPORAL");
	GLClasses::Shader& SSOcclusionCheckerboardShader = ShaderManager::GetShader("SCREENSPACE_OCCLUSION_CHECKERBOARD");

	GLClasses::Shader& BasicBlitShader = ShaderManager::GetShader("BLIT");
	GLClasses::Shader& RedOutputShader = ShaderManager::GetShader("RED");
	GLClasses::Shader& ConeTraceConvolutionShader = ShaderManager::GetShader("CONE_TRACE_CONVOLVE");

	GLClasses::Shader& SpatialFilter = ShaderManager::GetShader("SPATIAL");

	GLClasses::Shader& SVGFTemporalShader = ShaderManager::GetShader("SVGF_TEMPORAL");
	GLClasses::Shader& SVGFVarianceShader = ShaderManager::GetShader("SVGF_VARIANCE");

	GLClasses::Shader& GBufferDownsampleShader = ShaderManager::GetShader("GBUFFER_DOWNSAMPLER");
	GLClasses::Shader& VolumetricsShader = ShaderManager::GetShader("VOLUMETRICS");
	GLClasses::Shader& BilateralFilter4x4Shader = ShaderManager::GetShader("BILATERAL_FILTER");

	GLClasses::Shader& PostProcessCombineShader = ShaderManager::GetShader("POST_COMBINE");

	GLClasses::Shader& FXAAShader = ShaderManager::GetShader("FXAA");
	GLClasses::Shader& CASShader = ShaderManager::GetShader("CAS");


	// History
	glm::mat4 PreviousView;
	glm::mat4 PreviousProjection;
	glm::mat4 View;
	glm::mat4 Projection;
	glm::mat4 InverseView;
	glm::mat4 InverseProjection;

	// Probe Setup
	ProbeMap PlayerProbe(192);
	glm::vec3 PlayerProbeCapturePoint;

	// Temporal jitter
	GenerateJitterStuff();

	// SS Conetracing framebuffer
	GLuint GaussianMipFBO = 0;
	glGenFramebuffers(1, &GaussianMipFBO);
	glBindFramebuffer(GL_FRAMEBUFFER, GaussianMipFBO);
	GLenum Buffers[1] = { GL_COLOR_ATTACHMENT0 };
	glDrawBuffers(1, Buffers);
	glBindFramebuffer(GL_FRAMEBUFFER, 0);

	// Bloom Framebuffers 
	const float BloomResolution = 0.2f;
	BloomFBO BloomBuffer(16, 16);
	BloomFBO BloomBufferB(16, 16);
	BloomRenderer::Initialize();

	// Create voxel cascades 

	Voxelizer::CreateVolumes();

	GLClasses::Framebuffer* FinalDenoiseBufferPtr = &SpatialFilterBuffers[0];

	while (!glfwWindowShouldClose(app.GetWindow()))
	{
		RENDER_SCALE = RoundToNearest(RENDER_SCALE, 0.25f);
		CAS_Amount = RoundToNearest(CAS_Amount, 0.125f);

		// Normalize 
		//SunDirection = glm::normalize(SunDirection);

		// Matrices 
		PreviousProjection = Camera.GetProjectionMatrix();
		PreviousView = Camera.GetViewMatrix();

		// App update 
		PreviousSunDirection = SunDirection;
		app.OnUpdate();

		// Update current matrices 
		Projection = Camera.GetProjectionMatrix();
		View = Camera.GetViewMatrix();
		InverseProjection = glm::inverse(Camera.GetProjectionMatrix());
		InverseView = glm::inverse(Camera.GetViewMatrix());

		// Common uniform buffer
		CommonUniforms UniformBuffer = { View, Projection, InverseView, InverseProjection, PreviousProjection, PreviousView, glm::inverse(PreviousProjection), glm::inverse(PreviousView), (int)app.GetCurrentFrame() };

		float TrueResolutionW = app.GetWidth();
		float TrueResolutionH = app.GetHeight();

		float ScaledResolutionW = app.GetWidth() * RENDER_SCALE;
		float ScaledResolutionH = app.GetHeight() * RENDER_SCALE;

		// Resize buffers (to maintain aspect ratio)
		GBuffers[0].SetSize(ScaledResolutionW, ScaledResolutionH);
		GBuffers[1].SetSize(ScaledResolutionW, ScaledResolutionH);

		DownsampledGBuffer.SetSize(ScaledResolutionW / 2, ScaledResolutionH / 2);

		// Light combine (Direct + Indirect)
		LightingPass.SetSize(ScaledResolutionW, ScaledResolutionH);

		// Temporal AA (Full res)
		TAABuffers[0].SetSize(TrueResolutionW, TrueResolutionH);
		TAABuffers[1].SetSize(TrueResolutionW, TrueResolutionH);
		
		FXAABuffers[0].SetSize(TrueResolutionW, TrueResolutionH);
		FXAABuffers[1].SetSize(TrueResolutionW, TrueResolutionH);

		// Post-process combine  (Full res) 
		PostProcessCombined.SetSize(TrueResolutionW, TrueResolutionH);
		TonemappedFBO.SetSize(TrueResolutionW, TrueResolutionH);

		// Specular 
		SpecularIndirectBuffers[0].SetSize(ScaledResolutionW * IndirectRes * (CheckerboardIndirect ? 0.5f : 1.0f), ScaledResolutionH * IndirectRes);
		SpecularIndirectBuffers[1].SetSize(ScaledResolutionW * IndirectRes * (CheckerboardIndirect ? 0.5f : 1.0f), ScaledResolutionH * IndirectRes);
		SpecularIndirectTemporalBuffers[0].SetSize(ScaledResolutionW * IndirectRes, ScaledResolutionH * IndirectRes);
		SpecularIndirectTemporalBuffers[1].SetSize(ScaledResolutionW * IndirectRes, ScaledResolutionH * IndirectRes);

		// Indirect diffuse
		DiffuseIndirectTrace.SetSize(ScaledResolutionW * IndirectRes * (CheckerboardIndirect ? 0.5f : 1.0f), ScaledResolutionH * IndirectRes);
		DiffuseIndirectTemporalBuffers[0].SetSize(ScaledResolutionW * IndirectRes, ScaledResolutionH * IndirectRes);
		DiffuseIndirectTemporalBuffers[1].SetSize(ScaledResolutionW * IndirectRes, ScaledResolutionH * IndirectRes);
		SVGFVarianceResolve.SetSize(ScaledResolutionW * IndirectRes, ScaledResolutionH * IndirectRes);

		// Checkerboard
		IndirectCheckerUpscaled.SetSize(ScaledResolutionW * IndirectRes, ScaledResolutionH * IndirectRes);

		// Denoiser inputs/output buffers 
		SpatialFilterBuffers[0].SetSize(ScaledResolutionW / 2, ScaledResolutionH / 2);
		SpatialFilterBuffers[1].SetSize(ScaledResolutionW / 2, ScaledResolutionH / 2);

		// Buffers for SSCT 
		SpecularIndirectConeTraced.SetSize(ScaledResolutionW * IndirectRes, ScaledResolutionH * IndirectRes);
		SpecularIndirectConeTraceInput.SetSize(ScaledResolutionW * IndirectRes, ScaledResolutionH * IndirectRes);
		SpecularIndirectConeTraceInputAlternate.SetSize(ScaledResolutionW * IndirectRes, ScaledResolutionH * IndirectRes);

		// Volumetrics
		VolumetricsBuffers[0].SetSize(ScaledResolutionW * VolumetricsResolution, ScaledResolutionH * VolumetricsResolution);
		VolumetricsBuffers[1].SetSize(ScaledResolutionW * VolumetricsResolution, ScaledResolutionH * VolumetricsResolution);
		VolumetricsFiltered.SetSize(ScaledResolutionW * VolumetricsResolution, ScaledResolutionH * VolumetricsResolution);

		// Bloom 
		BloomBuffer.SetSize(ScaledResolutionW * BloomResolution, ScaledResolutionH * BloomResolution);
		BloomBufferB.SetSize(ScaledResolutionW * BloomResolution, ScaledResolutionH * BloomResolution);

		// RTAO
		ScreenspaceOcclusion.SetSize(ScaledResolutionW * ScreenspaceOcclusionRes * (ScreenspaceOcclusionCheckerboard ? 0.5f : 1.0f), ScaledResolutionH * ScreenspaceOcclusionRes);
		ScreenspaceOcclusionTemporalBuffers[0].SetSize(ScaledResolutionW * ScreenspaceOcclusionRes, ScaledResolutionH * ScreenspaceOcclusionRes);
		ScreenspaceOcclusionTemporalBuffers[1].SetSize(ScaledResolutionW * ScreenspaceOcclusionRes, ScaledResolutionH * ScreenspaceOcclusionRes);
		ScreenspaceOcclusionCheckerboardConstruct.SetSize(ScaledResolutionW * ScreenspaceOcclusionRes, ScaledResolutionH * ScreenspaceOcclusionRes);

		// Generate mipmaps 
		if (app.GetCurrentFrame() % 8 == 0) {

			// Mip passes

			glBindTexture(GL_TEXTURE_2D, SpecularIndirectConeTraceInput.GetTexture(0));
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
			glGenerateMipmap(GL_TEXTURE_2D);
			glBindTexture(GL_TEXTURE_2D, 0);

			glBindTexture(GL_TEXTURE_2D, SpecularIndirectConeTraceInputAlternate.GetTexture(0));
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
			glGenerateMipmap(GL_TEXTURE_2D);
			glBindTexture(GL_TEXTURE_2D, 0);
		}


		// Shadowmap update 
		if (app.GetCurrentFrame() % 12 == 0 || PreviousSunDirection != SunDirection)
		{
			// Shadow pass 
			RenderShadowMap(Shadowmap, SunDirection, EntityRenderList, Camera.GetViewProjection());
		}

		// Probe update 
		PlayerProbeCapturePoint = Camera.GetPosition();

		Lumen::ResetPolygonCount();

		for (int i = 0; i < ProbeUpdateRate; i++)
		{
			int RenderFace = (((app.GetCurrentFrame() * ProbeUpdateRate) + i) % 6);
			RenderProbe(PlayerProbe, RenderFace, PlayerProbeCapturePoint, EntityList, ProbeForwardShader, Skymap.GetID(), ScreenQuadVAO);
		}

		ProbePolygons = static_cast<uint32_t>(Lumen::QueryPolygonCount());

		// Ping pong framebuffers
		bool FrameCheckerStep = app.GetCurrentFrame() % 2 == 0;
		
		// Gbuffer 
		auto& GBuffer = FrameCheckerStep ? GBuffers[0] : GBuffers[1];
		auto& PrevGBuffer = FrameCheckerStep ? GBuffers[1] : GBuffers[0];

		// Specular 
		auto& SpecularIndirect = FrameCheckerStep ? SpecularIndirectBuffers[0] : SpecularIndirectBuffers[1];
		auto& PrevSpecularIndirect = FrameCheckerStep ? SpecularIndirectBuffers[1] : SpecularIndirectBuffers[0];
		auto& SpecularTemporal = FrameCheckerStep ? SpecularIndirectTemporalBuffers[0] : SpecularIndirectTemporalBuffers[1];
		auto& PrevSpecularTemporal = FrameCheckerStep ? SpecularIndirectTemporalBuffers[1] : SpecularIndirectTemporalBuffers[0];

		// Diffuse
		auto& DiffuseTemporal = FrameCheckerStep ? DiffuseIndirectTemporalBuffers[0] : DiffuseIndirectTemporalBuffers[1];
		auto& PrevDiffuseTemporal = FrameCheckerStep ? DiffuseIndirectTemporalBuffers[1] : DiffuseIndirectTemporalBuffers[0];

		// TAA
		auto& TAATemporal = FrameCheckerStep ? TAABuffers[0] : TAABuffers[1];
		auto& PrevTAATemporal = FrameCheckerStep ? TAABuffers[1] : TAABuffers[0];

		// Screenspace AO + Direct shadows
		auto& SSRTTemporal = FrameCheckerStep ? ScreenspaceOcclusionTemporalBuffers[0] : ScreenspaceOcclusionTemporalBuffers[1];
		auto& PrevSSRTTemporal = FrameCheckerStep ? ScreenspaceOcclusionTemporalBuffers[1] : ScreenspaceOcclusionTemporalBuffers[0];

		// Volumetrics
		auto& Volumetrics = FrameCheckerStep ? VolumetricsBuffers[0] : VolumetricsBuffers[1];
		auto& VolumetricsHistory = FrameCheckerStep ? VolumetricsBuffers[1] : VolumetricsBuffers[0];


		// Update voxel cascades->

		// Update voxel cascades for the first frames when everything is being updated
		if (app.GetCurrentFrame() < 6 || FullyDynamicVoxelization) {
			for (int i = 0; i < 6; i++) {
				Voxelizer::VoxelizeCascade(GetUpdateCascade(i), Camera.GetPosition(), Camera.GetProjectionMatrix(), Camera.GetViewMatrix(), Shadowmap.GetDepthTexture(), GetLightViewProjection(SunDirection), SunDirection, EntityList);
			}
		}

		if (!FullyDynamicVoxelization) {
			Voxelizer::VoxelizeCascade(GetUpdateCascade(app.GetCurrentFrame()), Camera.GetPosition(), Camera.GetProjectionMatrix(), Camera.GetViewMatrix(), Shadowmap.GetDepthTexture(), GetLightViewProjection(SunDirection), SunDirection, EntityList);
		}

		glm::vec3 PrevPosition = glm::vec3(UniformBuffer.InvPrevView[3]);
		bool CameraMoved = Camera.GetPosition() != PrevPosition;

		// Render GBuffer
		glEnable(GL_CULL_FACE);
		glEnable(GL_DEPTH_TEST);

		Lumen::ResetPolygonCount();

		GBuffer.Bind();
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		GBufferShader.Use();
		GBufferShader.SetMatrix4("u_ViewProjection", Camera.GetViewProjection());
		GBufferShader.SetInteger("u_AlbedoMap", 0);
		GBufferShader.SetInteger("u_NormalMap", 1);
		GBufferShader.SetInteger("u_RoughnessMap", 2);
		GBufferShader.SetInteger("u_MetalnessMap", 3);
		GBufferShader.SetInteger("u_MetalnessRoughnessMap", 5);
		GBufferShader.SetMatrix4("u_JitterMatrix", TAA ? GetTAAJitterMatrix(app.GetCurrentFrame(), GBuffer.GetDimensions()) : glm::mat4(1.0f));

		RenderEntityList(EntityRenderList, GBufferShader);
		UnbindEverything();

		MainRenderPolygons = static_cast<uint32_t>(Lumen::QueryPolygonCount());

		// Post processing passes here : 

		// Flags 
		glDisable(GL_CULL_FACE);
		glDisable(GL_DEPTH_TEST);

		// Downsample Depth (4x) to optimize screen space trace 

		GBufferDownsampleShader.Use();
		DownsampledGBuffer.Bind();

		GBufferDownsampleShader.SetInteger("u_Depth", 0);
		GBufferDownsampleShader.SetInteger("u_Normals", 1);
		GBufferDownsampleShader.SetInteger("u_PBR", 2);

		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, GBuffer.GetDepthBuffer());

		glActiveTexture(GL_TEXTURE1);
		glBindTexture(GL_TEXTURE_2D, GBuffer.GetTexture(1));

		glActiveTexture(GL_TEXTURE2);
		glBindTexture(GL_TEXTURE_2D, GBuffer.GetTexture(2));

		ScreenQuadVAO.Bind();
		glDrawArrays(GL_TRIANGLES, 0, 6);
		ScreenQuadVAO.Unbind();

		DownsampledGBuffer.Unbind();


		// Trace specular GI 

		SpecularIndirect.Bind();
		SpecularTraceShader.Use();
		SpecularTraceShader.SetVector3f("u_Incident", Camera.GetPosition());
		SpecularTraceShader.SetInteger("u_Depth", 0);
		SpecularTraceShader.SetInteger("u_Normals", 1);
		SpecularTraceShader.SetInteger("u_PBR", 2);
		SpecularTraceShader.SetInteger("u_ProbeAlbedo", 4);
		SpecularTraceShader.SetInteger("u_ProbeDepth", 5);
		SpecularTraceShader.SetInteger("u_ProbeNormals", 6);
		SpecularTraceShader.SetInteger("u_EnvironmentMap", 7);
		SpecularTraceShader.SetInteger("u_Shadowmap", 8);
		SpecularTraceShader.SetInteger("u_LFNormals", 9);
		SpecularTraceShader.SetInteger("u_Albedos", 11);
		SpecularTraceShader.SetInteger("u_LowResDepth", 12);
		SpecularTraceShader.SetInteger("u_PreviousFrameDiffuse", 15);
		SpecularTraceShader.SetInteger("u_Frame", app.GetCurrentFrame());
		SpecularTraceShader.SetBool("u_RoughSpecular", RoughSpecular);
		SpecularTraceShader.SetBool("u_Checker", CheckerboardIndirect);
		SpecularTraceShader.SetVector3f("u_SunDirection", SunDirection);
		SpecularTraceShader.SetMatrix4("u_SunShadowMatrix", GetLightViewProjection(SunDirection));
		SpecularTraceShader.SetVector2f("u_Jitter", GetTAAJitterSecondary(app.GetCurrentFrame()));
		SpecularTraceShader.SetVector2f("u_Dimensions", SpecularIndirect.GetDimensions());

		for (int i = 0; i < 6; i++) {
			std::string name = "u_ProbeCapturePoints[" + std::to_string(i) + "]";
			SpecularTraceShader.SetVector3f(name.c_str(), PlayerProbe.CapturePoints[i]);
		}

		// Bind voxel volumes and upload voxel volume info
		{
			int t = 0;

			for (int x = 16; x <= 16 + 5; x++) {

				std::string s = "u_VoxelVolumes[" + std::to_string(t) + "]";
				std::string s2 = "u_VoxelRanges[" + std::to_string(t) + "]";
				std::string s3 = "u_VoxelCenters[" + std::to_string(t) + "]";
				SpecularTraceShader.SetInteger(s.c_str(), x);
				SpecularTraceShader.SetFloat(s2.c_str(), Voxelizer::GetVolumeRanges()[t]);
				SpecularTraceShader.SetVector3f(s3.c_str(), Voxelizer::GetVolumeCenters()[t]);

				glActiveTexture(GL_TEXTURE0 + x);
				glBindTexture(GL_TEXTURE_3D, Voxelizer::GetVolumes()[t]);
				
				t++; 
			}
		}

		SetCommonUniforms(SpecularTraceShader, UniformBuffer);

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

		glActiveTexture(GL_TEXTURE9);
		glBindTexture(GL_TEXTURE_2D, GBuffer.GetTexture(3));

		glActiveTexture(GL_TEXTURE11);
		glBindTexture(GL_TEXTURE_2D, GBuffer.GetTexture(0));
		
		glActiveTexture(GL_TEXTURE12);
		glBindTexture(GL_TEXTURE_2D, DownsampledGBuffer.GetTexture());

		glActiveTexture(GL_TEXTURE15);
		glBindTexture(GL_TEXTURE_2D, SpatialFiltering ? FinalDenoiseBufferPtr->GetTexture(1) : DiffuseTemporal.GetTexture(0));

		// 16 -> 21 are used for the voxel volumes!

		ScreenQuadVAO.Bind();
		glDrawArrays(GL_TRIANGLES, 0, 6);
		ScreenQuadVAO.Unbind();




		// Trace diffuse GI

		DiffuseVXTrace.Use();
		DiffuseIndirectTrace.Bind();

		DiffuseVXTrace.SetInteger("u_Depth", 0);
		DiffuseVXTrace.SetInteger("u_LFNormals", 1);
		DiffuseVXTrace.SetInteger("u_BlueNoise", 2);
		DiffuseVXTrace.SetInteger("u_Skymap", 3);
		DiffuseVXTrace.SetInteger("u_ProbeAlbedo", 4);
		DiffuseVXTrace.SetInteger("u_ProbeDepth", 5);
		DiffuseVXTrace.SetInteger("u_ProbeNormals", 6);
		DiffuseVXTrace.SetInteger("u_Shadowmap", 7);

		DiffuseVXTrace.SetBool("u_Checker", CheckerboardIndirect);
		DiffuseVXTrace.SetVector3f("u_SunDirection", SunDirection);
		DiffuseVXTrace.SetMatrix4("u_SunShadowMatrix", GetLightViewProjection(SunDirection));

		for (int i = 0; i < 6; i++) {
			std::string name = "u_ProbeCapturePoints[" + std::to_string(i) + "]";
			DiffuseVXTrace.SetVector3f(name.c_str(), PlayerProbe.CapturePoints[i]);
		}

		SetCommonUniforms(DiffuseVXTrace, UniformBuffer);

		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, GBuffer.GetDepthBuffer());

		glActiveTexture(GL_TEXTURE1);
		glBindTexture(GL_TEXTURE_2D, GBuffer.GetTexture(3));

		glActiveTexture(GL_TEXTURE2);
		glBindTexture(GL_TEXTURE_2D, BlueNoise.GetTextureID());

		glActiveTexture(GL_TEXTURE3);
		glBindTexture(GL_TEXTURE_CUBE_MAP, Skymap.GetID());

		glActiveTexture(GL_TEXTURE4);
		glBindTexture(GL_TEXTURE_CUBE_MAP, PlayerProbe.m_CubemapTexture);

		glActiveTexture(GL_TEXTURE5);
		glBindTexture(GL_TEXTURE_CUBE_MAP, PlayerProbe.m_DepthCubemap);

		glActiveTexture(GL_TEXTURE6);
		glBindTexture(GL_TEXTURE_CUBE_MAP, PlayerProbe.NormalPBRPackedCubemap);

		glActiveTexture(GL_TEXTURE7);
		glBindTexture(GL_TEXTURE_2D, Shadowmap.GetDepthTexture());

		// Bind voxel volumes and upload volume data 
		{
			int temp = 0;
			int temp2 = 14;

			for (int x = 8; x <= 8 + 5; x++) {
				int t = temp;
				int t2 = temp2;
				std::string s = "u_VoxelVolumes[" + std::to_string(t) + "]";
				std::string s2 = "u_VoxelRanges[" + std::to_string(t) + "]";
				std::string s3 = "u_VoxelCenters[" + std::to_string(t) + "]";
				std::string s4 = "u_VoxelVolumesNormals[" + std::to_string(t) + "]";
				DiffuseVXTrace.SetInteger(s.c_str(), x);
				DiffuseVXTrace.SetFloat(s2.c_str(), Voxelizer::GetVolumeRanges()[t]);
				DiffuseVXTrace.SetVector3f(s3.c_str(), Voxelizer::GetVolumeCenters()[t]);
				DiffuseVXTrace.SetInteger(s4.c_str(), temp2);

				glActiveTexture(GL_TEXTURE0 + x);
				glBindTexture(GL_TEXTURE_3D, Voxelizer::GetVolumes()[t]);

				glActiveTexture(GL_TEXTURE0 + temp2);
				glBindTexture(GL_TEXTURE_3D, Voxelizer::GetVolumeNormals()[t]);

				temp++; temp2++;
			}
		}

		ScreenQuadVAO.Bind();
		glDrawArrays(GL_TRIANGLES, 0, 6);
		ScreenQuadVAO.Unbind();



		// Checkerboard reconstruction 

		if (CheckerboardIndirect) {
			IndirectCheckerboarder.Use();
			IndirectCheckerUpscaled.Bind();
			IndirectCheckerboarder.SetInteger("u_InputSpecular", 0);
			IndirectCheckerboarder.SetInteger("u_InputDiffuse", 1);
			IndirectCheckerboarder.SetInteger("u_Depth", 2);
			IndirectCheckerboarder.SetInteger("u_Normals", 3);
			SetCommonUniforms(IndirectCheckerboarder, UniformBuffer);

			glActiveTexture(GL_TEXTURE0);
			glBindTexture(GL_TEXTURE_2D, SpecularIndirect.GetTexture());

			glActiveTexture(GL_TEXTURE1);
			glBindTexture(GL_TEXTURE_2D, DiffuseIndirectTrace.GetTexture());

			glActiveTexture(GL_TEXTURE2);
			glBindTexture(GL_TEXTURE_2D, GBuffer.GetDepthBuffer());

			glActiveTexture(GL_TEXTURE3);
			glBindTexture(GL_TEXTURE_2D, GBuffer.GetTexture(3));

			ScreenQuadVAO.Bind();
			glDrawArrays(GL_TRIANGLES, 0, 6);
			ScreenQuadVAO.Unbind();

			IndirectCheckerUpscaled.Unbind();
		}


		// --Spatio-temporal filtering passes--

		// Temporally resolve lighting ->

		SpecularTemporalShader.Use();
		SpecularTemporal.Bind();
		SpecularTemporalShader.SetInteger("u_Specular", 0);
		SpecularTemporalShader.SetInteger("u_HistorySpecular", 1);
		SpecularTemporalShader.SetInteger("u_Depth", 2);
		SpecularTemporalShader.SetInteger("u_PreviousDepth", 3);
		SpecularTemporalShader.SetInteger("u_Normals", 4);
		SpecularTemporalShader.SetInteger("u_PrevTransversals", 6);
		SpecularTemporalShader.SetInteger("u_PBR", 7);
		SpecularTemporalShader.SetInteger("u_Frames", 8);
		SpecularTemporalShader.SetInteger("u_HitMask", 9);
		SpecularTemporalShader.SetBool("u_RoughSpecular", RoughSpecular);
		SpecularTemporalShader.SetBool("u_SpecularChecker", CheckerboardIndirect);
		SetCommonUniforms(SpecularTemporalShader, UniformBuffer);

		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, CheckerboardIndirect ? IndirectCheckerUpscaled.GetTexture() : SpecularIndirect.GetTexture());

		glActiveTexture(GL_TEXTURE1);
		glBindTexture(GL_TEXTURE_2D, PrevSpecularTemporal.GetTexture());

		glActiveTexture(GL_TEXTURE2);
		glBindTexture(GL_TEXTURE_2D, GBuffer.GetDepthBuffer());

		glActiveTexture(GL_TEXTURE3);
		glBindTexture(GL_TEXTURE_2D, PrevGBuffer.GetDepthBuffer());

		glActiveTexture(GL_TEXTURE4);
		glBindTexture(GL_TEXTURE_2D, GBuffer.GetTexture(3));

		glActiveTexture(GL_TEXTURE6);
		glBindTexture(GL_TEXTURE_2D, PrevSpecularIndirect.GetTexture());

		glActiveTexture(GL_TEXTURE7);
		glBindTexture(GL_TEXTURE_2D, GBuffer.GetTexture(2));

		glActiveTexture(GL_TEXTURE8);
		glBindTexture(GL_TEXTURE_2D, PrevSpecularTemporal.GetTexture(1));
		
		glActiveTexture(GL_TEXTURE9);
		glBindTexture(GL_TEXTURE_2D, SpecularIndirect.GetTexture(1));

		ScreenQuadVAO.Bind();
		glDrawArrays(GL_TRIANGLES, 0, 6);
		ScreenQuadVAO.Unbind();


		// Diffuse temporal

		SVGFTemporalShader.Use();
		DiffuseTemporal.Bind();

		SVGFTemporalShader.SetInteger("u_Current", 0);
		SVGFTemporalShader.SetInteger("u_History", 1);
		SVGFTemporalShader.SetInteger("u_Depth", 2);
		SVGFTemporalShader.SetInteger("u_PreviousDepth", 3);
		SVGFTemporalShader.SetInteger("u_Normals", 4);
		SVGFTemporalShader.SetInteger("u_Utility", 5);
		SVGFTemporalShader.SetInteger("u_PreviousNormals", 6);
		SVGFTemporalShader.SetInteger("u_LowResDepth", 7);
		SVGFTemporalShader.SetInteger("u_LowResNormals", 8);

		SetCommonUniforms(SVGFTemporalShader, UniformBuffer);

		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, CheckerboardIndirect ? IndirectCheckerUpscaled.GetTexture(1) : DiffuseIndirectTrace.GetTexture());

		glActiveTexture(GL_TEXTURE1);
		glBindTexture(GL_TEXTURE_2D, PrevDiffuseTemporal.GetTexture());

		glActiveTexture(GL_TEXTURE2);
		glBindTexture(GL_TEXTURE_2D, GBuffer.GetDepthBuffer());

		glActiveTexture(GL_TEXTURE3);
		glBindTexture(GL_TEXTURE_2D, PrevGBuffer.GetDepthBuffer());

		glActiveTexture(GL_TEXTURE4);
		glBindTexture(GL_TEXTURE_2D, GBuffer.GetTexture(3));

		glActiveTexture(GL_TEXTURE5);
		glBindTexture(GL_TEXTURE_2D, PrevDiffuseTemporal.GetTexture(1));

		glActiveTexture(GL_TEXTURE6);
		glBindTexture(GL_TEXTURE_2D, PrevGBuffer.GetTexture(3));

		glActiveTexture(GL_TEXTURE7);
		glBindTexture(GL_TEXTURE_2D, DownsampledGBuffer.GetTexture());

		glActiveTexture(GL_TEXTURE8);
		glBindTexture(GL_TEXTURE_2D, DownsampledGBuffer.GetTexture(1));

		
		ScreenQuadVAO.Bind();
		glDrawArrays(GL_TRIANGLES, 0, 6);
		ScreenQuadVAO.Unbind();
		
		DiffuseTemporal.Unbind();



		// -- Spatial Filtering Passes -- 
		
		if (SVGF) {

			// SVGF Variance estimation ->

			SVGFVarianceResolve.Bind();
			SVGFVarianceShader.Use();

			SVGFVarianceShader.SetInteger("u_LowResDepth", 0);
			SVGFVarianceShader.SetInteger("u_LowResNormals", 1);
			SVGFVarianceShader.SetInteger("u_Diffuse", 2);
			SVGFVarianceShader.SetInteger("u_Utility", 3);

			SetCommonUniforms(SVGFVarianceShader, UniformBuffer);

			glActiveTexture(GL_TEXTURE0);
			glBindTexture(GL_TEXTURE_2D, DownsampledGBuffer.GetTexture());

			glActiveTexture(GL_TEXTURE1);
			glBindTexture(GL_TEXTURE_2D, DownsampledGBuffer.GetTexture(1));

			glActiveTexture(GL_TEXTURE2);
			glBindTexture(GL_TEXTURE_2D, DiffuseTemporal.GetTexture());

			glActiveTexture(GL_TEXTURE3);
			glBindTexture(GL_TEXTURE_2D, DiffuseTemporal.GetTexture(1));

			ScreenQuadVAO.Bind();
			glDrawArrays(GL_TRIANGLES, 0, 6);
			ScreenQuadVAO.Unbind();

			SVGFVarianceResolve.Unbind();
		}


		// Wavelet filtering 

		if (SpatialFiltering) {

			FinalDenoiseBufferPtr = nullptr;

			const int Passes = 4;
			//const int StepSizes[5] = { 32, 16, 8, 4, 2 };
			const int StepSizes[4] = { 32, 12, 4, 2 };

			for (int Pass = 0; Pass < Passes; Pass++) {

				auto& CurrentBuffer = SpatialFilterBuffers[(Pass % 2 == 0) ? 0 : 1];
				auto& SpatialHistory = SpatialFilterBuffers[(Pass % 2 == 0) ? 1 : 0];

				FinalDenoiseBufferPtr = &CurrentBuffer;

				bool InitialPass = Pass == 0;

				CurrentBuffer.Bind();
				SpatialFilter.Use();

				SpatialFilter.SetInteger("u_Depth", 0);
				SpatialFilter.SetInteger("u_Normals", 1);
				SpatialFilter.SetInteger("u_PBR", 2);
				SpatialFilter.SetInteger("u_Specular", 3);
				SpatialFilter.SetInteger("u_BlueNoise", 4);
				SpatialFilter.SetInteger("u_SpecularFrames", 5);
				SpatialFilter.SetInteger("u_Diffuse", 6);
				SpatialFilter.SetInteger("u_Variance", 7);
				SpatialFilter.SetInteger("u_TemporalUtility", 8);
				SpatialFilter.SetInteger("u_StepSize", StepSizes[Pass]);
				SpatialFilter.SetInteger("u_Pass", Pass);
				SpatialFilter.SetInteger("u_TotalPasses", Passes);
				SpatialFilter.SetBool("u_SVGF", SVGF);

				SetCommonUniforms(SpatialFilter, UniformBuffer);

				glActiveTexture(GL_TEXTURE0);
				glBindTexture(GL_TEXTURE_2D, DownsampledGBuffer.GetTexture());

				glActiveTexture(GL_TEXTURE1);
				glBindTexture(GL_TEXTURE_2D, DownsampledGBuffer.GetTexture(1));

				glActiveTexture(GL_TEXTURE2);
				glBindTexture(GL_TEXTURE_2D, DownsampledGBuffer.GetTexture(2));

				glActiveTexture(GL_TEXTURE3);
				glBindTexture(GL_TEXTURE_2D, InitialPass ? SpecularTemporal.GetTexture() : SpatialHistory.GetTexture(0));

				glActiveTexture(GL_TEXTURE4);
				glBindTexture(GL_TEXTURE_2D, BlueNoise.GetTextureID());

				glActiveTexture(GL_TEXTURE5);
				glBindTexture(GL_TEXTURE_2D, SpecularTemporal.GetTexture(1));

				glActiveTexture(GL_TEXTURE6);
				glBindTexture(GL_TEXTURE_2D, InitialPass ? (SVGF ? SVGFVarianceResolve.GetTexture() : DiffuseTemporal.GetTexture()) : SpatialHistory.GetTexture(1));

				glActiveTexture(GL_TEXTURE7);
				glBindTexture(GL_TEXTURE_2D, InitialPass ? SVGFVarianceResolve.GetTexture(1) : SpatialHistory.GetTexture(2));

				glActiveTexture(GL_TEXTURE8);
				glBindTexture(GL_TEXTURE_2D, DiffuseTemporal.GetTexture(1));

				ScreenQuadVAO.Bind();
				glDrawArrays(GL_TRIANGLES, 0, 6);
				ScreenQuadVAO.Unbind();
			}
		}

		auto& FinalDenoisedBuffer = FinalDenoiseBufferPtr ? *FinalDenoiseBufferPtr : SpatialFilterBuffers[0];


		// Screenspace raytracing -> (SSRTAO & SS Contact Shadows)

		// Screenspace RTAO / Screenspace direct shadows 

		ScreenspaceOcclusion.Bind();
		SSOcclusionTraceShader.Use();

		SSOcclusionTraceShader.SetInteger("u_Depth", 0);
		SSOcclusionTraceShader.SetInteger("u_Normals", 1);
		SSOcclusionTraceShader.SetInteger("u_PBR", 2);
		SSOcclusionTraceShader.SetBool("u_Shadow", DoScreenspaceShadow);
		SSOcclusionTraceShader.SetBool("u_AO", DoScreenspaceAO);
		SSOcclusionTraceShader.SetBool("u_Checkerboard", ScreenspaceOcclusionCheckerboard);
		SSOcclusionTraceShader.SetVector2f("u_Jitter", GetTAAJitterSecondary(app.GetCurrentFrame()));
		SSOcclusionTraceShader.SetVector2f("u_Dimensions", ScreenspaceOcclusion.GetDimensions());
		SetCommonUniforms(SSOcclusionTraceShader, UniformBuffer);

		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, GBuffer.GetDepthBuffer());

		glActiveTexture(GL_TEXTURE1);
		glBindTexture(GL_TEXTURE_2D, GBuffer.GetTexture(3));

		glActiveTexture(GL_TEXTURE2);
		glBindTexture(GL_TEXTURE_2D, GBuffer.GetTexture(2));


		ScreenQuadVAO.Bind();
		glDrawArrays(GL_TRIANGLES, 0, 6);
		ScreenQuadVAO.Unbind();

		ScreenspaceOcclusion.Unbind();

		// Checkerboard reconstruct 
		if (ScreenspaceOcclusionCheckerboard) {

			ScreenspaceOcclusionCheckerboardConstruct.Bind();
			SSOcclusionCheckerboardShader.Use();

			SSOcclusionCheckerboardShader.SetInteger("u_Input", 0);
			SSOcclusionCheckerboardShader.SetInteger("u_InputDirect", 1);
			SSOcclusionCheckerboardShader.SetInteger("u_Depth", 2);
			SSOcclusionCheckerboardShader.SetInteger("u_Normals", 3);
			SetCommonUniforms(SSOcclusionCheckerboardShader, UniformBuffer);

			glActiveTexture(GL_TEXTURE0);
			glBindTexture(GL_TEXTURE_2D, ScreenspaceOcclusion.GetTexture());

			glActiveTexture(GL_TEXTURE1);
			glBindTexture(GL_TEXTURE_2D, ScreenspaceOcclusion.GetTexture(1));

			glActiveTexture(GL_TEXTURE2);
			glBindTexture(GL_TEXTURE_2D, GBuffer.GetDepthBuffer());

			glActiveTexture(GL_TEXTURE3);
			glBindTexture(GL_TEXTURE_2D, GBuffer.GetTexture(3));

			ScreenQuadVAO.Bind();
			glDrawArrays(GL_TRIANGLES, 0, 6);
			ScreenQuadVAO.Unbind();
		}



		// Temporal resolve : 

		SSRTTemporal.Bind();
		SSOcclusionTemporalShader.Use();
		SSOcclusionTemporalShader.SetInteger("u_Current", 0);
		SSOcclusionTemporalShader.SetInteger("u_History", 1);
		SSOcclusionTemporalShader.SetInteger("u_Depth", 2);
		SSOcclusionTemporalShader.SetInteger("u_PreviousDepth", 3);
		SSOcclusionTemporalShader.SetInteger("u_Normals", 4);
		SSOcclusionTemporalShader.SetInteger("u_CurrentDirect", 5);
		SSOcclusionTemporalShader.SetInteger("u_HistoryDirect", 6);
		SetCommonUniforms(SSOcclusionTemporalShader, UniformBuffer);

		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, ScreenspaceOcclusionCheckerboard ? ScreenspaceOcclusionCheckerboardConstruct.GetTexture() : ScreenspaceOcclusion.GetTexture(0));

		glActiveTexture(GL_TEXTURE1);
		glBindTexture(GL_TEXTURE_2D, PrevSSRTTemporal.GetTexture(0));

		glActiveTexture(GL_TEXTURE2);
		glBindTexture(GL_TEXTURE_2D, GBuffer.GetDepthBuffer());

		glActiveTexture(GL_TEXTURE3);
		glBindTexture(GL_TEXTURE_2D, PrevGBuffer.GetDepthBuffer());

		glActiveTexture(GL_TEXTURE4);
		glBindTexture(GL_TEXTURE_2D, GBuffer.GetTexture(3));

		glActiveTexture(GL_TEXTURE5);
		glBindTexture(GL_TEXTURE_2D, ScreenspaceOcclusionCheckerboard ? ScreenspaceOcclusionCheckerboardConstruct.GetTexture(1) : ScreenspaceOcclusion.GetTexture(1));

		glActiveTexture(GL_TEXTURE6);
		glBindTexture(GL_TEXTURE_2D, PrevSSRTTemporal.GetTexture(1));

		ScreenQuadVAO.Bind();
		glDrawArrays(GL_TRIANGLES, 0, 6);
		ScreenQuadVAO.Unbind();

		SSRTTemporal.Unbind();

		// Volumetrics 

		if (VolumetricLighting) {
			Volumetrics.Bind();

			VolumetricsShader.Use();
			VolumetricsShader.SetInteger("u_LowResDepth", 0);
			VolumetricsShader.SetInteger("u_LFNormals", 1);
			VolumetricsShader.SetInteger("u_Shadowmap", 2);
			VolumetricsShader.SetInteger("u_HistoryVolumetrics", 3);
			VolumetricsShader.SetInteger("u_HistoryDepth", 4);
			VolumetricsShader.SetVector3f("u_SunDirection", SunDirection);
			VolumetricsShader.SetMatrix4("u_SunShadowMatrix", GetLightViewProjection(SunDirection));
			VolumetricsShader.SetFloat("u_SunVLStrength", SunVolumetricsStrength);
			VolumetricsShader.SetBool("u_VolumetricsTemporal", VolumetricsTemporal);
			SetCommonUniforms(VolumetricsShader, UniformBuffer);

			glActiveTexture(GL_TEXTURE0);
			glBindTexture(GL_TEXTURE_2D, DownsampledGBuffer.GetTexture());

			glActiveTexture(GL_TEXTURE1);
			glBindTexture(GL_TEXTURE_2D, GBuffer.GetTexture(3));

			glActiveTexture(GL_TEXTURE2);
			glBindTexture(GL_TEXTURE_2D, Shadowmap.GetDepthTexture());

			glActiveTexture(GL_TEXTURE3);
			glBindTexture(GL_TEXTURE_2D, VolumetricsHistory.GetTexture());

			glActiveTexture(GL_TEXTURE4);
			glBindTexture(GL_TEXTURE_2D, PrevGBuffer.GetDepthBuffer());

			ScreenQuadVAO.Bind();
			glDrawArrays(GL_TRIANGLES, 0, 6);
			ScreenQuadVAO.Unbind();

			Volumetrics.Unbind();

			// Filter volumetrics (4x4 bilateral atrous filter)

			BilateralFilter4x4Shader.Use();
			VolumetricsFiltered.Bind();

			BilateralFilter4x4Shader.SetInteger("u_Input", 0);
			BilateralFilter4x4Shader.SetInteger("u_LowResDepth", 1);
			BilateralFilter4x4Shader.SetFloat("u_SigmaB", CameraMoved || !VolumetricsTemporal ? 1.414f : 0.8f);
			SetCommonUniforms(BilateralFilter4x4Shader, UniformBuffer);

			glActiveTexture(GL_TEXTURE0);
			glBindTexture(GL_TEXTURE_2D, Volumetrics.GetTexture());

			glActiveTexture(GL_TEXTURE1);
			glBindTexture(GL_TEXTURE_2D, DownsampledGBuffer.GetTexture());

			ScreenQuadVAO.Bind();
			glDrawArrays(GL_TRIANGLES, 0, 6);
			ScreenQuadVAO.Unbind();

			VolumetricsFiltered.Unbind();
		}



		// Lighting combine pass : 

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
		LightingShader.SetInteger("u_RTAO", 9);
		LightingShader.SetInteger("u_ScreenspaceShadows", 10);
		LightingShader.SetInteger("u_IndirectDiffuse", 17);
		LightingShader.SetInteger("u_Volumetrics", 18);

		LightingShader.SetInteger("u_LFNormals", 20);
		
		LightingShader.SetBool("u_DirectSSShadows", DoScreenspaceShadow);
		LightingShader.SetBool("u_VolumetricsEnabled", VolumetricLighting);
		LightingShader.SetFloat("u_RTAOStrength", ScreenspaceAOStrength);

		LightingShader.SetMatrix4("u_LightVP", GetLightViewProjection(SunDirection));
		LightingShader.SetVector2f("u_Dims", glm::vec2(app.GetWidth(), app.GetHeight()));

		LightingShader.SetVector3f("u_LightDirection", SunDirection);
		LightingShader.SetVector3f("u_ViewerPosition", Camera.GetPosition());

		
		for (int i = 0; i < 6; i++) {
			std::string name = "u_ProbeCapturePoints[" + std::to_string(i) + "]";
			LightingShader.SetVector3f(name.c_str(), PlayerProbe.CapturePoints[i]);
		}

		
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
		glBindTexture(GL_TEXTURE_2D, SSCT ? SpecularIndirectConeTraced.GetTexture() : (SpatialFiltering ? FinalDenoisedBuffer.GetTexture() : SpecularTemporal.GetTexture()));
		
		glActiveTexture(GL_TEXTURE9);
		glBindTexture(GL_TEXTURE_2D, SSRTTemporal.GetTexture());
		
		glActiveTexture(GL_TEXTURE10);
		glBindTexture(GL_TEXTURE_2D, SSRTTemporal.GetTexture(1));
		
		glActiveTexture(GL_TEXTURE17);
		glBindTexture(GL_TEXTURE_2D, SpatialFiltering ? FinalDenoisedBuffer.GetTexture(1) : DiffuseTemporal.GetTexture());

		glActiveTexture(GL_TEXTURE18);
		glBindTexture(GL_TEXTURE_2D, VolumetricsFilter ? VolumetricsFiltered.GetTexture() : Volumetrics.GetTexture());

		glActiveTexture(GL_TEXTURE20);
		glBindTexture(GL_TEXTURE_2D, GBuffer.GetTexture(3));

		// SLOTS 11 - 15 ARE USED FOR VOXEL VOLUME TEXTURES

		// Bind voxel volumes 
		
		int temp__ = 0;

		for (int x = 11; x <= 11 + 5; x++) {
			int t = temp__;
			std::string s = "u_VoxelVolumes[" + std::to_string(t) + "]";
			std::string s2 = "u_VoxelRanges[" + std::to_string(t) + "]";
			std::string s3 = "u_VoxelCenters[" + std::to_string(t) + "]";
			LightingShader.SetInteger(s.c_str(), x);
			LightingShader.SetFloat(s2.c_str(), Voxelizer::GetVolumeRanges()[t]);
			LightingShader.SetVector3f(s3.c_str(), Voxelizer::GetVolumeCenters()[t]);

			glActiveTexture(GL_TEXTURE0 + x);
			glBindTexture(GL_TEXTURE_3D, Voxelizer::GetVolumes()[t]);

			temp__++;
		}

		LightingShader.SetInteger("u_VoxelVolume0", 12);


		ScreenQuadVAO.Bind();
		glDrawArrays(GL_TRIANGLES, 0, 6);
		ScreenQuadVAO.Unbind();

		// Temporal AA Resolve :

		TAAShader.Use();
		TAATemporal.Bind();

		TAAShader.SetInteger("u_CurrentColorTexture", 0);
		TAAShader.SetInteger("u_PreviousColorTexture", 1);
		TAAShader.SetInteger("u_DepthTexture", 2);
		TAAShader.SetInteger("u_PreviousDepthTexture", 3);
		TAAShader.SetVector2f("u_CurrentJitter", GetTAAJitter(app.GetCurrentFrame()));
		TAAShader.SetFloat("u_TAAUConfidenceExponent", TAAUConfidenceExponent);
		TAAShader.SetBool("u_Enabled", TAA);
		TAAShader.SetBool("u_TAAU", TAAU);
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

		// Bloom 
		GLuint BrightTex = 0;

		if (DoBloom) {
			BloomRenderer::RenderBloom(TAATemporal.GetTexture(), GBuffer.GetTexture(2), BloomBuffer, BloomBufferB, BrightTex, true);
		}


		// Combine post process

		PostProcessCombineShader.Use();
		PostProcessCombined.Bind();

		PostProcessCombineShader.SetInteger("u_Texture", 0);
		PostProcessCombineShader.SetInteger("u_BloomMips[0]", 1);
		PostProcessCombineShader.SetInteger("u_BloomMips[1]", 2);
		PostProcessCombineShader.SetInteger("u_BloomMips[2]", 3);
		PostProcessCombineShader.SetInteger("u_BloomMips[3]", 4);
		PostProcessCombineShader.SetInteger("u_BloomMips[4]", 5);
		PostProcessCombineShader.SetInteger("u_BloomBrightTexture", 6);
		PostProcessCombineShader.SetBool("u_BloomEnabled", DoBloom);

		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, TAATemporal.GetTexture(0));

		glActiveTexture(GL_TEXTURE1);
		glBindTexture(GL_TEXTURE_2D, BloomBuffer.m_Mips[0]);

		glActiveTexture(GL_TEXTURE2);
		glBindTexture(GL_TEXTURE_2D, BloomBuffer.m_Mips[1]);

		glActiveTexture(GL_TEXTURE3);
		glBindTexture(GL_TEXTURE_2D, BloomBuffer.m_Mips[2]);

		glActiveTexture(GL_TEXTURE4);
		glBindTexture(GL_TEXTURE_2D, BloomBuffer.m_Mips[3]);

		glActiveTexture(GL_TEXTURE5);
		glBindTexture(GL_TEXTURE_2D, BloomBuffer.m_Mips[4]);

		glActiveTexture(GL_TEXTURE6);
		glBindTexture(GL_TEXTURE_2D, BrightTex);

		ScreenQuadVAO.Bind();
		glDrawArrays(GL_TRIANGLES, 0, 6);
		ScreenQuadVAO.Unbind();

		PostProcessCombined.Unbind();

		// FXAA

		GLClasses::Framebuffer* PrevFXAAFBO = &PostProcessCombined;
		GLClasses::Framebuffer* CurrentFXAAFBO = &FXAABuffers[0];

		for (int x = 0; x < FXAA_PASSES; x++) {

			CurrentFXAAFBO->Bind();

			FXAAShader.Use();
			FXAAShader.SetInteger("u_MainTexture", 0);
			FXAAShader.SetFloat("u_RenderScale", RENDER_SCALE);

			glActiveTexture(GL_TEXTURE0);
			glBindTexture(GL_TEXTURE_2D, PrevFXAAFBO->GetTexture());

			ScreenQuadVAO.Bind();
			glDrawArrays(GL_TRIANGLES, 0, 6);
			ScreenQuadVAO.Unbind();

			CurrentFXAAFBO->Unbind();

			PrevFXAAFBO = CurrentFXAAFBO;

			if (CurrentFXAAFBO == &FXAABuffers[0]) { CurrentFXAAFBO = &FXAABuffers[1]; }
			else if (CurrentFXAAFBO == &FXAABuffers[1]) { CurrentFXAAFBO = &FXAABuffers[0]; } 
		}

		CurrentFXAAFBO = PrevFXAAFBO;//FXAA_PASSES > 0 ? CurrentFXAAFBO : &PostProcessCombined;


		// Tonemap + Gamma correct + Output

		TonemappedFBO.Bind();

		FinalShader.Use();
		FinalShader.SetInteger("u_MainTexture", 0);

		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, CurrentFXAAFBO->GetTexture(0));

		ScreenQuadVAO.Bind();
		glDrawArrays(GL_TRIANGLES, 0, 6);
		ScreenQuadVAO.Unbind();

		TonemappedFBO.Unbind();

		// CAS + convert to srgb + blit to screen

		glBindFramebuffer(GL_FRAMEBUFFER, 0);
		glViewport(0, 0, app.GetWidth(), app.GetHeight());

		CASShader.Use();

		CASShader.SetInteger("u_Texture", 0);
		CASShader.SetFloat("u_SharpenAmount", CAS_Amount);
		CASShader.SetBool("u_Upscaled", RENDER_SCALE < 0.95f);

		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, TonemappedFBO.GetTexture(0));

		ScreenQuadVAO.Bind();
		glDrawArrays(GL_TRIANGLES, 0, 6);
		ScreenQuadVAO.Unbind();



		// Finish :
		app.FinishFrame();
		GLClasses::DisplayFrameRate(app.GetWindow(), "Lumen ");
	}
}
