#include "BVHConstructor.h"

namespace Lumen {
	namespace BVH {

		static std::vector<Triangle> BVHTriangles;
		static std::vector<Triangle> Triangles;
		static std::vector<Node> BVHNodes;

		inline glm::vec3 Vec3Min(const glm::vec3& a, const glm::vec3& b) {
			return glm::vec3(glm::min(a.x, b.x), glm::min(a.y, b.y), glm::min(a.z, b.z));
		}

		inline glm::vec3 Vec3Max(const glm::vec3& a, const glm::vec3& b) {
			return glm::vec3(glm::max(a.x, b.x), glm::max(a.y, b.y), glm::max(a.z, b.z));
		}

		inline bool ShouldBeLeaf(const uint& Length) {
			return Length <= 4;
		}

		int FindLongestAxis(const Bounds& bounds) {
			glm::vec3 Diff = bounds.Max - bounds.Min;

			float Max = glm::max(glm::max(Diff.x, Diff.y), Diff.z);

			if (Max == Diff.x) { return 0; }
			if (Max == Diff.y) { return 1; }
			if (Max == Diff.z) { return 2; }
			
			throw "What";
			
			return 0;
		}

		void BuildNode(Node* node) {

			uint StartIndex = node->StartIndex;
			uint EndIndex = node->Length + StartIndex;
			uint Length = EndIndex - StartIndex;

			bool IsLeaf = ShouldBeLeaf(Length);

			if (IsLeaf) {
				return;
			}

			else {

				// Create 2 new nodes 

				Node& LeftNode = BVHNodes.emplace_back();
				node->LeftChild = BVHNodes.size() - 1;
				LeftNode.NodeIndex = node->LeftChild;

				Node& RightNode = BVHNodes.emplace_back();
				node->RightChild = BVHNodes.size() - 1;
				RightNode.NodeIndex = node->RightChild;

				int LongestAxis = FindLongestAxis(node->NodeBounds);
				float Border = node->NodeBounds.GetCenter()[LongestAxis];

				Bounds LeftNodeBounds;
				Bounds RightNodeBounds;

				// Split based on spatial median 

				uint LeftIndexItr = StartIndex;
				uint RightIndexItr = EndIndex;
				uint FirstTriangleIndex = 0;

				// Split/Sort triangles into 2 sets 
				// One for the left node and one for the right node
				while (LeftIndexItr < RightIndexItr) {

					// Forward iteration
					while (LeftIndexItr < RightIndexItr) {

						Triangle& CurrentTriangle = Triangles[LeftIndexItr];

						if (CurrentTriangle.v0.position[LongestAxis] > Border &&
							CurrentTriangle.v1.position[LongestAxis] > Border &&
							CurrentTriangle.v2.position[LongestAxis] > Border) {

							break;
						}

						FirstTriangleIndex++;
						LeftIndexItr++;
					}

					// Backward iteration
					while (LeftIndexItr < RightIndexItr) {

						Triangle& CurrentTriangle = Triangles[RightIndexItr];

						if (CurrentTriangle.v0.position[LongestAxis] <= Border &&
							CurrentTriangle.v1.position[LongestAxis] <= Border &&
							CurrentTriangle.v2.position[LongestAxis] <= Border) {

							break;
						}

						RightIndexItr--;
					}

					// Swap triangles
					if (LeftIndexItr < RightIndexItr) {

						Triangle Temp = Triangles[LeftIndexItr];
						Triangles[LeftIndexItr] = Triangles[RightIndexItr];
						Triangles[RightIndexItr] = Temp;

					}


				}


				// Set node properties 

				LeftNode.LeftChild = 0;
				LeftNode.RightChild = 0;
				LeftNode.StartIndex = node->StartIndex;
				LeftNode.Length = FirstTriangleIndex;
				LeftNode.ParentNode = node->NodeIndex;
				LeftNode.IsLeftNode = true;

				// Find bounding boxes 

				// Left node bounding box 
				glm::vec3 LeftMin = glm::vec3(10000.0f);
				glm::vec3 LeftMax = glm::vec3(-10000.0f);

				for (int i = LeftNode.StartIndex; i < LeftNode.StartIndex + LeftNode.Length; i++)
				{
					Bounds CurrentBounds = Triangles.at(i).GetBounds();
					LeftMin = Vec3Min(CurrentBounds.Min, LeftMin);
					LeftMax = Vec3Min(CurrentBounds.Max, LeftMax);
				}

				LeftNode.NodeBounds = Bounds(LeftMin, LeftMax);

				// Initialize right node
				RightNode.LeftChild = 0;
				RightNode.RightChild = 0;
				RightNode.StartIndex = node->StartIndex + FirstTriangleIndex;
				RightNode.Length = node->Length - FirstTriangleIndex;
				RightNode.ParentNode = node->NodeIndex;
				RightNode.IsLeftNode = true;


				// Find right bounding box
				glm::vec3 RightMin = glm::vec3(10000.0f);
				glm::vec3 RightMax = glm::vec3(-10000.0f);

				for (int i = RightNode.StartIndex; i < RightNode.StartIndex + RightNode.Length; i++)
				{
					Bounds CurrentBounds = Triangles.at(i).GetBounds();
					RightMin = Vec3Min(CurrentBounds.Min, RightMin);
					RightMax = Vec3Min(CurrentBounds.Max, RightMax);
				}

				RightNode.NodeBounds = Bounds(RightMin, RightMax);


				if (!ShouldBeLeaf(RightNode.Length))
				{
					BuildNode(&RightNode);
				}

				if (!ShouldBeLeaf(LeftNode.Length)) {
					BuildNode(&LeftNode);
				}
			}



		}

		void BuildBVH(std::vector<Vertex>& Vertices, std::vector<GLuint>& Indices)
		{
			// First, generate triangles from vertices to make everything easier to work with 

			for (int x = 0; x < Indices.size(); x += 3)
			{
				Vertex v0 = Vertices.at(Indices.at(x + 0));
				Vertex v1 = Vertices.at(Indices.at(x + 1));
				Vertex v2 = Vertices.at(Indices.at(x + 2));

				Triangle tri = { (v0), (v1), (v2) };
				Triangles.push_back(tri);
			}


			// Create bounding box 

			glm::vec3 InitialMin = glm::vec3(10000.0f);
			glm::vec3 InitialMax = glm::vec3(-10000.0f);

			for (int i = 0; i < Triangles.size(); i++)
			{
				Bounds CurrentBounds = Triangles.at(i).GetBounds();
				InitialMin = Vec3Min(CurrentBounds.Min, InitialMin);
				InitialMax = Vec3Min(CurrentBounds.Max, InitialMax);
			}

		}
	}
}