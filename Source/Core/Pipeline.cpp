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
static bool FXAA = true;

// Specular settings
static const float SpecularIndirectRes = 0.5f;
static bool SSCT = false;
static bool RoughSpecular = true;
static bool CheckerboardIndirect = true;

// Diffuse settings 
static const float DiffuseIndirectRes = 0.5f;

// AO 
const float ScreenspaceOcclusionRes = 0.5f;
static bool DoScreenspaceAO = true;
static bool DoScreenspaceShadow = false;
static float ScreenspaceAOStrength = 0.75f;
static bool ScreenspaceOcclusionCheckerboard = true;

// Probe update rate 
static int ProbeUpdateRate = 1;

// Spatial
static bool SpatialFiltering = true;

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
		float camera_speed = 0.525f * 2.0f;

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

		ImGui::NewLine();

		ImGui::Checkbox("Spatial Filtering?", &SpatialFiltering);

		ImGui::NewLine();

		ImGui::NewLine();
		
		ImGui::Text("Specular Resolution : %f on each axis", SpecularIndirectRes);
		ImGui::Text("Specular Upsample Resolution : %f on each axis", SpecularIndirectRes);
		ImGui::Checkbox("Rough Specular?", &RoughSpecular);
		ImGui::Checkbox("Screenspace cone tracing?", &SSCT);

		ImGui::NewLine();

		ImGui::Text("Indirect Diffuse Resolution : %f on each axis", DiffuseIndirectRes);
		ImGui::Checkbox("Indirect Checkerboarding? (Resolution effectively halved)", &CheckerboardIndirect);
		ImGui::NewLine();

		ImGui::Text("RTAO/Screenspace shadows resolve resolution : %f on each axis", ScreenspaceOcclusionRes);
		ImGui::Checkbox("Do Screenspace RTAO?", &DoScreenspaceAO);
		ImGui::Checkbox("Do Screenspace Direct Shadows?", &DoScreenspaceShadow);
		ImGui::Checkbox("Checkerboard render?", &ScreenspaceOcclusionCheckerboard);

		if (DoScreenspaceAO)
			ImGui::SliderFloat("Screenspace RTAO Strength", &ScreenspaceAOStrength, 0.1f, 2.0f);

		ImGui::NewLine();
		ImGui::NewLine();
		ImGui::Text("Antialiasing : ");
		ImGui::NewLine();
		ImGui::Checkbox("TAA", &TAA);
		ImGui::Checkbox("FXAA 3.11", &FXAA);
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
GLClasses::Framebuffer SpecularIndirectBuffers[2]{ GLClasses::Framebuffer(16, 16, { {GL_RGBA16F, GL_RGBA, GL_FLOAT, true, true} }, false, false), GLClasses::Framebuffer(16, 16, { {GL_RGBA16F, GL_RGBA, GL_FLOAT, true, true} }, false, false) };
GLClasses::Framebuffer SpecularIndirectTemporalBuffers[2] = { GLClasses::Framebuffer(16, 16, {{GL_RGBA16F, GL_RGBA, GL_FLOAT, true, true}, {GL_R16F, GL_RED, GL_FLOAT, false, false} }, false, false), GLClasses::Framebuffer(16, 16, {{GL_RGBA16F, GL_RGBA, GL_FLOAT, true, true},{GL_R16F, GL_RED, GL_FLOAT, false, false}}, false, false) };
GLClasses::Framebuffer SpecularIndirectConeTraceInput(16, 16, { GL_RGB16F, GL_RGB, GL_FLOAT, true, true }, false, false);
GLClasses::Framebuffer SpecularIndirectConeTraceInputAlternate(16, 16, { GL_RGB16F, GL_RGB, GL_FLOAT, true, true }, false, false);
GLClasses::Framebuffer SpecularIndirectConeTraced(16, 16, { GL_RGB16F, GL_RGB, GL_FLOAT, true, true }, false, false);

// Indirect diffuse 
GLClasses::Framebuffer DiffuseIndirectTrace(16, 16, { GL_RGBA16F, GL_RGBA, GL_FLOAT, true, true }, false, false);
GLClasses::Framebuffer DiffuseIndirectTemporalBuffers[2] = { GLClasses::Framebuffer(16, 16, {{GL_RGBA16F, GL_RGBA, GL_FLOAT, true, true},{GL_R16F, GL_RED, GL_FLOAT, false, false} }, false, false), GLClasses::Framebuffer(16, 16, { {GL_RGBA16F, GL_RGBA, GL_FLOAT, true, true},{GL_R16F, GL_RED, GL_FLOAT, false, false} }, false, false) };

// Common indirect buffers 
GLClasses::Framebuffer IndirectCheckerUpscaled(16, 16, { {GL_RGBA16F, GL_RGBA, GL_FLOAT, true, true}, {GL_RGBA16F, GL_RGBA, GL_FLOAT, true, true} }, false, false);

// TAA Buffers
GLClasses::Framebuffer TAABuffers[2] = { GLClasses::Framebuffer(16, 16, {GL_RGB16F, GL_RGB, GL_FLOAT, true, true}, false, false), GLClasses::Framebuffer(16, 16, {GL_RGB16F, GL_RGB, GL_FLOAT, true, true}, false, false) };

// Spatial Buffers 
GLClasses::Framebuffer SpatialFilterBuffers[2] = { GLClasses::Framebuffer(16, 16, {{ GL_RGBA16F, GL_RGBA, GL_FLOAT, true, true },{ GL_RGBA16F, GL_RGBA, GL_FLOAT, true, true }}, false, false), GLClasses::Framebuffer(16, 16, {{ GL_RGBA16F, GL_RGBA, GL_FLOAT, true, true },{ GL_RGBA16F, GL_RGBA, GL_FLOAT, true, true }}, false, false) };


int GetUpdateCascade(int Frame) {

	int Rand = rand() % 6000;
	int Rand2 = rand() % 1000;
	int Rand3 = rand() % 1000;

	if (Rand < 3250) {
		if (Rand2 < 600) {
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

	Object SecondaryLargeModel;
	FileLoader::LoadModelFile(&SecondaryLargeModel, "Models/Lucy/LucyModel.obj");

	Object SecondaryModel;
	FileLoader::LoadModelFile(&SecondaryModel, "Models/dragon/dragon.obj");
	
	Entity SecondaryEntity(&SecondaryLargeModel);
	SecondaryEntity.m_Model = glm::scale(glm::mat4(1.0f), glm::vec3(0.6f));
	SecondaryEntity.m_Model = glm::translate(SecondaryEntity.m_Model, glm::vec3(0.0f, 0.8f, 0.0f));
	SecondaryEntity.m_EntityRoughness = 0.4f;
	SecondaryEntity.m_EntityMetalness = 1.0f;

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

	Entity SecondaryEntity4(&SecondaryModel);
	SecondaryEntity4.m_Model = glm::scale(glm::mat4(1.0f), glm::vec3(EntityScale * 0.5));
	SecondaryEntity4.m_Model = glm::translate(SecondaryEntity4.m_Model, glm::vec3(50.0f, 3.0f, 0.0f) * (1.0f / (EntityScale * 0.5f)));
	SecondaryEntity4.m_EmissiveAmount = 10.0f;

	Entity SecondaryEntity5(&SecondaryModel);
	SecondaryEntity5.m_Model = glm::scale(glm::mat4(1.0f), glm::vec3(EntityScale * 0.5f));
	SecondaryEntity5.m_Model = glm::translate(SecondaryEntity5.m_Model, glm::vec3(-50, 3.0f, 0.0f) * (1.0f / (EntityScale * 0.5f)));
	SecondaryEntity5.m_EmissiveAmount = 10.0f;

	// Construct BVH

	//BVH::Node* SponzaRootBVHNode = BVH::BuildBVH(Sponza);

	// Entity list
	std::vector<Entity*> EntityRenderList = { &MainModel, &SecondaryEntity, &SecondaryEntity0, &SecondaryEntity1, &SecondaryEntity2, &SecondaryEntity3, &SecondaryEntity4, &SecondaryEntity5 };
	auto& EntityList = EntityRenderList;

	// Clear CPU side vertex/index data (After bvh construction ofc.) 
	std::vector<Object*> ObjectList = { &Sponza, &SecondaryLargeModel, &SecondaryModel };
	
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
	GLClasses::Shader& ProbeSpecularShader = ShaderManager::GetShader("PROBE_SPECULAR");
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

	GLClasses::Shader& GBufferDownsampleShader = ShaderManager::GetShader("GBUFFER_DOWNSAMPLER");


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

	// Conetracing framebuffer
	GLuint GaussianMipFBO = 0;
	glGenFramebuffers(1, &GaussianMipFBO);
	glBindFramebuffer(GL_FRAMEBUFFER, GaussianMipFBO);
	GLenum Buffers[1] = { GL_COLOR_ATTACHMENT0 };
	glDrawBuffers(1, Buffers);
	glBindFramebuffer(GL_FRAMEBUFFER, 0);

	// Create voxel cascades 

	Voxelizer::CreateVolumes();

	while (!glfwWindowShouldClose(app.GetWindow()))
	{
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

		// Resize buffers (to maintain aspect ratio)
		GBuffers[0].SetSize(app.GetWidth(), app.GetHeight());
		GBuffers[1].SetSize(app.GetWidth(), app.GetHeight());

		DownsampledGBuffer.SetSize(app.GetWidth() / 2, app.GetHeight() / 2);

		// Light combine (Direct + Indirect)
		LightingPass.SetSize(app.GetWidth(), app.GetHeight());

		// Temporal AA
		TAABuffers[0].SetSize(app.GetWidth(), app.GetHeight());
		TAABuffers[1].SetSize(app.GetWidth(), app.GetHeight());

		// Specular 
		SpecularIndirectBuffers[0].SetSize(app.GetWidth() * SpecularIndirectRes * (CheckerboardIndirect ? 0.5f : 1.0f), app.GetHeight() * SpecularIndirectRes);
		SpecularIndirectBuffers[1].SetSize(app.GetWidth() * SpecularIndirectRes * (CheckerboardIndirect ? 0.5f : 1.0f), app.GetHeight() * SpecularIndirectRes);
		SpecularIndirectTemporalBuffers[0].SetSize(app.GetWidth() * SpecularIndirectRes, app.GetHeight() * SpecularIndirectRes);
		SpecularIndirectTemporalBuffers[1].SetSize(app.GetWidth() * SpecularIndirectRes, app.GetHeight() * SpecularIndirectRes);

		// Indirect diffuse
		DiffuseIndirectTrace.SetSize(app.GetWidth() * DiffuseIndirectRes * (CheckerboardIndirect ? 0.5f : 1.0f), app.GetHeight() * DiffuseIndirectRes);
		DiffuseIndirectTemporalBuffers[0].SetSize(app.GetWidth() * DiffuseIndirectRes, app.GetHeight() * DiffuseIndirectRes);
		DiffuseIndirectTemporalBuffers[1].SetSize(app.GetWidth() * DiffuseIndirectRes, app.GetHeight() * DiffuseIndirectRes);

		IndirectCheckerUpscaled.SetSize(app.GetWidth() * DiffuseIndirectRes, app.GetHeight() * DiffuseIndirectRes);

		// Buffers for SSCT 
		SpecularIndirectConeTraced.SetSize(app.GetWidth() * SpecularIndirectRes, app.GetHeight() * SpecularIndirectRes);
		SpecularIndirectConeTraceInput.SetSize(app.GetWidth() * SpecularIndirectRes, app.GetHeight() * SpecularIndirectRes);
		SpecularIndirectConeTraceInputAlternate.SetSize(app.GetWidth() * SpecularIndirectRes, app.GetHeight() * SpecularIndirectRes);

		// RTAO
		ScreenspaceOcclusion.SetSize(app.GetWidth() * ScreenspaceOcclusionRes * (ScreenspaceOcclusionCheckerboard ? 0.5f : 1.0f), app.GetHeight() * ScreenspaceOcclusionRes);
		ScreenspaceOcclusionTemporalBuffers[0].SetSize(app.GetWidth() * ScreenspaceOcclusionRes, app.GetHeight() * ScreenspaceOcclusionRes);
		ScreenspaceOcclusionTemporalBuffers[1].SetSize(app.GetWidth() * ScreenspaceOcclusionRes, app.GetHeight() * ScreenspaceOcclusionRes);
		ScreenspaceOcclusionCheckerboardConstruct.SetSize(app.GetWidth() * ScreenspaceOcclusionRes, app.GetHeight() * ScreenspaceOcclusionRes);

		// Spatial 

		SpatialFilterBuffers[0].SetSize(app.GetWidth() / 2, app.GetHeight() / 2);
		SpatialFilterBuffers[1].SetSize(app.GetWidth() / 2, app.GetHeight() / 2);

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

		// Update voxel cascades->

		// Update voxel cascades for the first frames when everything is being updated
		if (app.GetCurrentFrame() < 6) {
			std::cout << "\nUpdating all voxel cascades!";
			for (int i = 0; i < 6; i++) {
				Voxelizer::VoxelizeCascade(GetUpdateCascade(i), Camera.GetPosition(), Camera.GetProjectionMatrix(), Camera.GetViewMatrix(), Shadowmap.GetDepthTexture(), GetLightViewProjection(SunDirection), SunDirection, EntityList);
			}
		}

		if (app.GetCurrentFrame() % 1 == 0) {
			Voxelizer::VoxelizeCascade(GetUpdateCascade(app.GetCurrentFrame()), Camera.GetPosition(), Camera.GetProjectionMatrix(), Camera.GetViewMatrix(), Shadowmap.GetDepthTexture(), GetLightViewProjection(SunDirection), SunDirection, EntityList);
		}

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

		// Specular Indirect lighting trace :

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
		ProbeSpecularShader.SetInteger("u_LFNormals", 9);
		ProbeSpecularShader.SetInteger("u_Albedos", 11);
		ProbeSpecularShader.SetInteger("u_LowResDepth", 12);
		ProbeSpecularShader.SetInteger("u_Frame", app.GetCurrentFrame());
		ProbeSpecularShader.SetBool("u_RoughSpecular", RoughSpecular);
		ProbeSpecularShader.SetBool("u_Checker", CheckerboardIndirect);
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

		glActiveTexture(GL_TEXTURE9);
		glBindTexture(GL_TEXTURE_2D, GBuffer.GetTexture(3));

		glActiveTexture(GL_TEXTURE11);
		glBindTexture(GL_TEXTURE_2D, GBuffer.GetTexture(0));
		
		glActiveTexture(GL_TEXTURE12);
		glBindTexture(GL_TEXTURE_2D, DownsampledGBuffer.GetTexture());

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
		DiffuseVXTrace.SetBool("u_Checker", CheckerboardIndirect);

		SetCommonUniforms(DiffuseVXTrace, UniformBuffer);

		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, GBuffer.GetDepthBuffer());

		glActiveTexture(GL_TEXTURE1);
		glBindTexture(GL_TEXTURE_2D, GBuffer.GetTexture(3));

		glActiveTexture(GL_TEXTURE2);
		glBindTexture(GL_TEXTURE_2D, BlueNoise.GetTextureID());

		glActiveTexture(GL_TEXTURE3);
		glBindTexture(GL_TEXTURE_CUBE_MAP, Skymap.GetID());

		{
			int temp = 0;

			for (int x = 8; x <= 8 + 5; x++) {
				int t = temp;
				std::string s = "u_VoxelVolumes[" + std::to_string(t) + "]";
				std::string s2 = "u_VoxelRanges[" + std::to_string(t) + "]";
				std::string s3 = "u_VoxelCenters[" + std::to_string(t) + "]";
				DiffuseVXTrace.SetInteger(s.c_str(), x);
				DiffuseVXTrace.SetFloat(s2.c_str(), Voxelizer::GetVolumeRanges()[t]);
				DiffuseVXTrace.SetVector3f(s3.c_str(), Voxelizer::GetVolumeCenters()[t]);

				glActiveTexture(GL_TEXTURE0 + x);
				glBindTexture(GL_TEXTURE_3D, Voxelizer::GetVolumes()[t]);

				temp++;
			}
		}

		ScreenQuadVAO.Bind();
		glDrawArrays(GL_TRIANGLES, 0, 6);
		ScreenQuadVAO.Unbind();

		// Checkerboard reconstruction ->

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




		// Specular temporal resolve :

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
		SpecularTemporalShader.SetBool("u_RoughSpecular", RoughSpecular);
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
		SVGFTemporalShader.SetInteger("u_Frames", 5);
		SVGFTemporalShader.SetInteger("u_PreviousNormals", 6);

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
		
		ScreenQuadVAO.Bind();
		glDrawArrays(GL_TRIANGLES, 0, 6);
		ScreenQuadVAO.Unbind();
		
		DiffuseTemporal.Unbind();


		// Wavelet filtering 

		GLClasses::Framebuffer* FinalDenoiseBufferPtr = nullptr;

		if (SpatialFiltering) {

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
				SpatialFilter.SetInteger("u_StepSize", StepSizes[Pass]);
				SpatialFilter.SetInteger("u_Pass", Pass);
				SpatialFilter.SetInteger("u_TotalPasses", Passes);

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
				glBindTexture(GL_TEXTURE_2D, InitialPass ? DiffuseTemporal.GetTexture() : SpatialHistory.GetTexture(1));


				ScreenQuadVAO.Bind();
				glDrawArrays(GL_TRIANGLES, 0, 6);
				ScreenQuadVAO.Unbind();
			}
		}

		auto& FinalDenoisedBuffer = FinalDenoiseBufferPtr ? *FinalDenoiseBufferPtr : SpatialFilterBuffers[0];

		// Screenspace raytracing ->

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
		LightingShader.SetBool("u_DirectSSShadows", DoScreenspaceShadow);
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

		// Tonemap + Gamma correct + Output :

		glBindFramebuffer(GL_FRAMEBUFFER, 0);
		glViewport(0, 0, app.GetWidth(), app.GetHeight());

		FinalShader.Use();
		FinalShader.SetInteger("u_MainTexture", 0);
		FinalShader.SetBool("u_FXAA", FXAA);

		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, TAATemporal.GetTexture(0));

		ScreenQuadVAO.Bind();
		glDrawArrays(GL_TRIANGLES, 0, 6);
		ScreenQuadVAO.Unbind();

		// Finish :
		app.FinishFrame();
		GLClasses::DisplayFrameRate(app.GetWindow(), "Lumen ");
	}
}
