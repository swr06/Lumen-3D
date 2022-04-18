#pragma once

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <glad/glad.h>
#include <iostream>
#include "GLClasses/DepthBuffer.h"
#include "GLClasses/Shader.h"
#include "Mesh.h"
#include "Object.h"
#include "Entity.h"

namespace Lumen {
	void InitShadowMapRenderer();
	void RenderShadowMap(GLClasses::DepthBuffer& depthbuffer, const glm::vec3& player_pos, const glm::vec3& sun_dir, const glm::vec2& sunrot, std::vector<Entity*> entites, glm::mat4 m);
	glm::mat4 GetLightViewProjection(const glm::vec3& sun_dir);
}