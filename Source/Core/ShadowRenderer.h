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
#include "Shadowmap.h"

namespace Lumen {
	void RenderShadowMap(Shadowmap& Shadowmap, const glm::vec3& Origin, glm::vec3 SunDirection, const std::vector<Entity*>& Entities, float Distance, glm::mat4& Projection, glm::mat4& View);
}