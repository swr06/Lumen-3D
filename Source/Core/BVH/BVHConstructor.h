#pragma once

/*

VERY WIP !

*/

#include <iostream>
#include <vector>

#include <cmath>

#include <glm/glm.hpp>

#include "../Utils/Vertex.h"

#include "../Object.h"

namespace Lumen {
	namespace BVH {
		typedef uint32_t uint;

		const float ARBITRARY_MAX = 1000000.0f;
		const float ARBITRARY_MIN = -1000000.0f;

		class Bounds {

		public :

			Bounds() : Min(glm::vec3(ARBITRARY_MAX)), Max(glm::vec3(ARBITRARY_MIN)) {}
			Bounds(const glm::vec3& min, const glm::vec3& max) : Min(glm::vec3(min)), Max(glm::vec3(max)) {}

			glm::vec3 Min;
			glm::vec3 Max;

			glm::vec3 GetCenter() {
				return (Min + Max) / 2.0f;
			}

			glm::vec3 GetExtent() {
				return Max - Min;
			}
		};

		struct Node {
			Bounds NodeBounds;
			uint NodeIndex = 0;
			uint StartIndex = 0;
			uint Length = 0;
			uint LeftChild = 0;
			uint RightChild = 0;
			uint ParentNode = 0;
			bool IsLeftNode = false;
		};

		struct FlattenedNode {

		};

		struct Triangle 
		{
			Vertex v0, v1, v2;

			inline glm::vec3 GetCentroid() const noexcept {
				return (v0.position + v1.position + v2.position) / 3.0f;
			}

			inline Bounds GetBounds() const noexcept {
				glm::vec3 min, max;
				min.x = glm::min(glm::min(v0.position.x, v1.position.x), v2.position.x);
				min.y = glm::min(glm::min(v0.position.y, v1.position.y), v2.position.y);
				min.z = glm::min(glm::min(v0.position.z, v1.position.z), v2.position.z);

				max.x = glm::max(glm::max(v0.position.x, v1.position.x), v2.position.x);
				max.y = glm::max(glm::max(v0.position.y, v1.position.y), v2.position.y);
				max.z = glm::max(glm::max(v0.position.z, v1.position.z), v2.position.z);

				return Bounds(min, max);
			}
		};

		Node* BuildBVH(std::vector<Vertex>& Vertices, std::vector<GLuint>& Indices);
		Node* BuildBVH(Object& object);
	}
};