#include "ModelRenderer.h"

#include <glm/glm.hpp>

void Lumen::RenderEntity(Entity& entity, GLClasses::Shader& shader)
{
	Object* object = entity.m_Object;

	shader.SetMatrix4("u_ModelMatrix", entity.m_Model);
	shader.SetMatrix3("u_NormalMatrix", glm::mat3(glm::transpose(glm::inverse(entity.m_Model))));

	int DrawCalls = 0;

	for (auto& e : object->m_Meshes)
	{


		const Mesh* mesh = &e;

		if (mesh->m_AlbedoMap.GetID() != 0)
		{
			mesh->m_AlbedoMap.Bind(0);
		}

		if (mesh->m_NormalMap.GetID() != 0)
		{
			mesh->m_NormalMap.Bind(1);
		}

		else {
			std::cout << "knife";
		}

		if (mesh->m_RoughnessMap.GetID() != 0)
		{
			mesh->m_RoughnessMap.Bind(2);
		}

		if (mesh->m_MetalnessMap.GetID() != 0)
		{
			mesh->m_MetalnessMap.Bind(3);
		}

		if (mesh->m_AmbientOcclusionMap.GetID() != 0)
		{
			mesh->m_AmbientOcclusionMap.Bind(4);
		}

		shader.SetBool("u_UsesGLTFPBR", false);

		if (mesh->TexturePaths[5].size() > 0 && mesh->m_MetalnessRoughnessMap.GetID() > 0) {

			shader.SetBool("u_UsesGLTFPBR", true);
			mesh->m_MetalnessRoughnessMap.Bind(5);
		}

		const GLClasses::VertexArray& VAO = mesh->m_VertexArray;
		VAO.Bind();

		if (mesh->m_Indexed)
		{
			DrawCalls++;
			glDrawElements(GL_TRIANGLES, mesh->m_IndicesCount, GL_UNSIGNED_INT, 0);
		}

		else
		{
			DrawCalls++;
			glDrawArrays(GL_TRIANGLES, 0, mesh->m_VertexCount);
		}

		VAO.Unbind();
	}


	if (std::fmod(glfwGetTime(), 0.5f) < 0.001f)
	{
		std::cout << "\nDRAW CALLS : " << DrawCalls;
	}
}
