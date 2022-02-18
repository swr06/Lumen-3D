#include "BVHConstructor.h"

#include <stack>

namespace Lumen {
	namespace BVH {

		static std::vector<Triangle> Triangles;
		static uint64_t TotalIterations = 0;
		static uint32_t LastIndex = 0;

		inline glm::vec3 Vec3Min(const glm::vec3& a, const glm::vec3& b) {
			return glm::vec3(glm::min(a.x, b.x), glm::min(a.y, b.y), glm::min(a.z, b.z));
		}

		inline glm::vec3 Vec3Max(const glm::vec3& a, const glm::vec3& b) {
			return glm::vec3(glm::max(a.x, b.x), glm::max(a.y, b.y), glm::max(a.z, b.z));
		}

		inline bool ShouldBeLeaf(uint Length) {
			return Length <= 256;
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

		void BuildNode(Node* RootNode) {

			std::stack<Node*> NodeStack;

			NodeStack.push(RootNode);

			while (!NodeStack.empty()) {

				std::cout << "\nProcess index : " << LastIndex;

				Node* node = NodeStack.top();
				NodeStack.pop();

				TotalIterations++;
				uint StartIndex = node->StartIndex;
				uint EndIndex = node->Length + StartIndex;
				uint Length = EndIndex - StartIndex;

				bool IsLeaf = ShouldBeLeaf(Length);

				if (IsLeaf) {
					continue;
				}

				else {

					// Create 2 new nodes 

					Node* LeftNodePtr = new Node;
					Node& LeftNode = *LeftNodePtr;
					LastIndex++;
					node->LeftChild = LastIndex;
					LeftNode.NodeIndex = node->LeftChild;

					Node* RightNodePtr = new Node;
					Node& RightNode = *RightNodePtr;
					LastIndex++;
					node->RightChild = LastIndex;
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

					LeftNode.StartIndex = node->StartIndex;
					LeftNode.Length = FirstTriangleIndex;
					LeftNode.LeftChild = 0;
					LeftNode.RightChild = 0;
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
						LeftMax = Vec3Max(CurrentBounds.Max, LeftMax);
					}

					LeftNode.NodeBounds = Bounds(LeftMin, LeftMax);

					// Initialize right node
					RightNode.LeftChild = 0;
					RightNode.RightChild = 0;
					RightNode.StartIndex = node->StartIndex + FirstTriangleIndex;
					RightNode.Length = node->Length - FirstTriangleIndex;
					RightNode.ParentNode = node->NodeIndex;
					RightNode.IsLeftNode = false;


					// Find right bounding box
					glm::vec3 RightMin = glm::vec3(10000.0f);
					glm::vec3 RightMax = glm::vec3(-10000.0f);

					for (int i = RightNode.StartIndex; i < RightNode.StartIndex + RightNode.Length; i++)
					{
						Bounds CurrentBounds = Triangles.at(i).GetBounds();
						RightMin = Vec3Min(CurrentBounds.Min, RightMin);
						RightMax = Vec3Max(CurrentBounds.Max, RightMax);
					}

					RightNode.NodeBounds = Bounds(RightMin, RightMax);


					if (!ShouldBeLeaf(RightNode.Length)) {
						//BuildNode(&RightNode);
						NodeStack.push(RightNodePtr);
					}

					if (!ShouldBeLeaf(LeftNode.Length)) {
						//BuildNode(&LeftNode);
						NodeStack.push(LeftNodePtr);
					}

				}
			}


		}

		Node* BuildBVH(std::vector<Vertex>& Vertices, std::vector<GLuint>& Indices)
		{
			// First, generate triangles from vertices to make everything easier to work with 

			std::cout << "\nGenerating Triangles..";

			for (int x = 0; x < Indices.size(); x += 3)
			{
				Vertex v0 = Vertices.at(Indices.at(x + 0));
				Vertex v1 = Vertices.at(Indices.at(x + 1));
				Vertex v2 = Vertices.at(Indices.at(x + 2));

				Triangle tri = { (v0), (v1), (v2) };
				Triangles.push_back(tri);
			}

			std::cout << "\nGenerated Triangles!";

			std::cout << "\Creating initial node bounding box..";


			// Create bounding box 

			glm::vec3 InitialMin = glm::vec3(10000.0f);
			glm::vec3 InitialMax = glm::vec3(-10000.0f);

			for (int i = 0; i < Triangles.size(); i++)
			{
				Bounds CurrentBounds = Triangles.at(i).GetBounds();
				InitialMin = Vec3Min(CurrentBounds.Min, InitialMin);
				InitialMax = Vec3Max(CurrentBounds.Max, InitialMax);
			}

			Node* RootNodePtr = new Node;
			Node& RootNode = *RootNodePtr;

			RootNode.LeftChild = 0;
			RootNode.RightChild = 0;
			RootNode.StartIndex = 0;
			RootNode.Length = Triangles.size();
			RootNode.ParentNode = 0;
			RootNode.IsLeftNode = false;
			RootNode.NodeBounds = Bounds(InitialMin, InitialMax);

			std::cout << "\nGenerated bounding box!";

			// Recursively construct 
			BuildNode(&RootNode);

			// Output debug stats 

			std::cout << "\n\n\n";
			std::cout << "--BVH Construction Info--";
			std::cout << "Triangle Count : " << Triangles.size();
			std::cout << "Node Count : " << LastIndex;
			std::cout << "\n\n\n";

			return &RootNode;
		}


		Node* BuildBVH(Object& object)
		{
			TotalIterations = 0;
				
			// First, generate triangles from vertices to make everything easier to work with 

			std::cout << "\nGenerating Triangles..";

			for (auto& Mesh : object.m_Meshes) {

				auto& Indices = Mesh.m_Indices;
				auto& Vertices = Mesh.m_Vertices;

				for (int x = 0; x < Indices.size(); x += 3)
				{
					Vertex v0 = Vertices.at(Indices.at(x + 0));
					Vertex v1 = Vertices.at(Indices.at(x + 1));
					Vertex v2 = Vertices.at(Indices.at(x + 2));

					Triangle tri = { (v0), (v1), (v2) };
					Triangles.push_back(tri);
				}
			}

			std::cout << "\nGenerated Triangles!";

			// Create bounding box 
			std::cout << "\Creating initial node bounding box..";

			glm::vec3 InitialMin = glm::vec3(10000.0f);
			glm::vec3 InitialMax = glm::vec3(-10000.0f);

			for (int i = 0; i < Triangles.size(); i++)
			{
				Bounds CurrentBounds = Triangles.at(i).GetBounds();
				InitialMin = Vec3Min(CurrentBounds.Min, InitialMin);
				InitialMax = Vec3Max(CurrentBounds.Max, InitialMax);
			}

			Node* RootNodePtr = new Node;
			Node& RootNode = *RootNodePtr;

			RootNode.LeftChild = 0;
			RootNode.RightChild = 0;
			RootNode.StartIndex = 0;
			RootNode.Length = Triangles.size() - 1;
			RootNode.ParentNode = 0;
			RootNode.IsLeftNode = false;
			RootNode.NodeBounds = Bounds(InitialMin, InitialMax);

			std::cout << "\nGenerated bounding box!";

			// Recursively construct 
			BuildNode(&RootNode);

			// Output debug stats 

			std::cout << "\n\n\n";
			std::cout << "--BVH Construction Info--";
			std::cout << "Triangle Count : " << Triangles.size();
			std::cout << "Node Count : " << Triangles.size();
			std::cout << "\n\n\n";

			return &RootNode;


		}

	}




}