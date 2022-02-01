#pragma once

#include <iostream>
#include <glm/glm.hpp>
#include <glad/glad.h>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <cmath>

namespace Lumen {

	void GenerateJitterStuff();
	glm::vec2 GetTAAJitter(int CurrentFrame);
	glm::vec2 GetTAAJitterSecondary(int CurrentFrame);
}