#pragma once

#include <iostream>

#include <glm/glm.hpp>

#include <glad/glad.h>

#include "GLClasses/ComputeShader.h"

#include "Entity.h"

namespace Lumen {
	namespace Voxelizer {

		void CreateVolumes();

		void VoxelizeCascade(int Cascade, glm::vec3 Position, const glm::mat4& Projection, const glm::mat4& View, GLuint Shadowmap, const glm::mat4& LVP, const glm::vec3&, const std::vector<Entity*>& EntityList);

		GLuint* GetVolumes();
		float* GetVolumeRanges();
		glm::vec3* GetVolumeCenters();

		void RecompileShaders();
	}
}