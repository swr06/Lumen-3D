#pragma once

#include <iostream>
#include <vector>

#include <glm/glm.hpp>

namespace Lumen {
	namespace BVH {
		typedef uint32_t uint;

		struct Bounds {
			glm::vec3 Min;
			glm::vec3 Max;
		};

		struct Node {
			Bounds bounds;
			bool is_child;

		};

		struct FlattenedNode {

		};
	}
};