#include "ModelFileLoader.h"

#include <assimp/pbrmaterial.h>
#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>
#include <chrono>

#include "MeshOptimizer.h"
#include <string>
#include <vector>
#include <array>

#include <fstream>

#define PACK_U16(lsb, msb) ((uint16_t) ( ((uint16_t)(lsb) & 0xFF) | (((uint16_t)(msb) & 0xFF) << 8) ))

/* Model Loader
Uses the assimp model loading library to load the models. It uses a recursive model to process the meshes and materials
*/

namespace Lumen
{
	namespace FileLoader
	{
		bool is_gltf = false;

		static bool FileExists(const std::string& str) {
			std::ifstream file(str);

			if (file.is_open() && file.good())
			{
				file.close();
				return true;
			}

			return false;

		}

		void LoadMaterialTextures(aiMesh* mesh, aiMaterial* mat, Mesh* _mesh, const std::string& path)
		{
			std::filesystem::path pth(path);

			std::string texture_path = pth.parent_path().string().c_str();
			aiString material_name;
			aiString diffuse_texture;
			aiString specular_texture;
			aiString normal_texture;
			aiString roughness_texture;
			aiString metallic_texture;
			aiString ao_texture;
			_mesh->m_IsGLTF = is_gltf;

			if (mat->GetTexture(AI_MATKEY_GLTF_PBRMETALLICROUGHNESS_BASE_COLOR_TEXTURE, &diffuse_texture) == aiReturn_SUCCESS)
			{
				std::string pth = texture_path + "/" + diffuse_texture.C_Str();

				if (FileExists(pth)) {
					_mesh->TexturePaths[0] = pth;
					_mesh->RawTexturePaths[0] = diffuse_texture.C_Str();
				}
			}

			else if (mat->GetTexture(aiTextureType_DIFFUSE, 0, &diffuse_texture) == aiReturn_SUCCESS)
			{
				std::string pth = texture_path + "/" + diffuse_texture.C_Str();

				if (FileExists(pth) && diffuse_texture.length > 0) {

					_mesh->TexturePaths[0] = pth;
					_mesh->RawTexturePaths[0] = diffuse_texture.C_Str();
				}
			}

			if (mat->GetTexture(aiTextureType_NORMALS, 0, &normal_texture) == aiReturn_SUCCESS)
			{
				std::string pth = texture_path + "/" + normal_texture.C_Str();

				if (FileExists(pth)) {

					_mesh->TexturePaths[1] = pth;
					_mesh->RawTexturePaths[1] = normal_texture.C_Str();

				}
			}


			if (mat->GetTexture(AI_MATKEY_GLTF_PBRMETALLICROUGHNESS_METALLICROUGHNESS_TEXTURE, &metallic_texture) == aiReturn_SUCCESS)
			{
				std::string pth = texture_path + "/" + metallic_texture.C_Str();

				if (FileExists(pth)) {

					_mesh->TexturePaths[5] = pth;
					_mesh->RawTexturePaths[5] = metallic_texture.C_Str();
				}
			}

			else {
				if (mat->GetTexture(aiTextureType_METALNESS, 0, &metallic_texture) == aiReturn_SUCCESS)
				{
					std::string pth = texture_path + "/" + metallic_texture.C_Str();

					if (FileExists(pth)) {

						_mesh->TexturePaths[3] = pth;
						_mesh->RawTexturePaths[3] = metallic_texture.C_Str();
					}
				}

				if (mat->GetTexture(aiTextureType_DIFFUSE_ROUGHNESS, 0, &roughness_texture) == aiReturn_SUCCESS)
				{
					std::string pth = texture_path + "/" + roughness_texture.C_Str();

					if (FileExists(pth)) {

						_mesh->TexturePaths[2] = pth;
						_mesh->RawTexturePaths[2] = roughness_texture.C_Str();
					}
				}
			}

			if (mat->GetTexture(aiTextureType_AMBIENT_OCCLUSION, 0, &ao_texture) == aiReturn_SUCCESS)
			{
				std::string pth = texture_path + "/" + ao_texture.C_Str();

				if (FileExists(pth)) {

					_mesh->TexturePaths[4] = pth;
					_mesh->RawTexturePaths[4] = ao_texture.C_Str();

				}
			}
		}

		void ProcessAssimpMesh(aiMesh* mesh, const aiScene* scene, Object* object, const std::string& pth, const glm::vec4& col, const glm::vec3& reflectivity, glm::vec3 emissive_color)
		{
			Mesh& _mesh = object->GenerateMesh();
			std::vector<Vertex>& vertices = _mesh.m_Vertices;
			std::vector<GLuint>& indices = _mesh.m_Indices;

			for (int i = 0; i < mesh->mNumVertices; i++)
			{
				Vertex vt;
				vt.position = glm::vec3(
					mesh->mVertices[i].x,
					mesh->mVertices[i].y,
					mesh->mVertices[i].z
				);

				glm::vec3 vnormal = glm::vec3(0.0f), vtan = glm::vec3(0.0f);

				if (mesh->HasNormals())
				{
					vnormal = glm::vec3(
						mesh->mNormals[i].x,
						mesh->mNormals[i].y,
						mesh->mNormals[i].z
					);
				}

				if (mesh->mTextureCoords[0])
				{
					vt.texcoords = glm::packHalf2x16(glm::vec2(
						mesh->mTextureCoords[0][i].x,
						mesh->mTextureCoords[0][i].y
					));

					if (mesh->mTangents)
					{
						vtan.x = mesh->mTangents[i].x;
						vtan.y = mesh->mTangents[i].y;
						vtan.z = mesh->mTangents[i].z;
					}
				}

				else
				{
					vt.texcoords = glm::packHalf2x16(glm::vec2(0.0f, 0.0f));
				}

				glm::uvec3 data;
				data.x = glm::packHalf2x16(glm::vec2(vnormal.x, vnormal.y));
				data.y = glm::packHalf2x16(glm::vec2(vnormal.z, vtan.x));
				data.z = glm::packHalf2x16(glm::vec2(vtan.y, vtan.z));
				vt.normal_tangent_data = data;
				vertices.push_back(vt);
			}

			/* Push back the indices */
			for (int i = 0; i < mesh->mNumFaces; i++)
			{
				aiFace face = mesh->mFaces[i];

				for (unsigned int j = 0; j < face.mNumIndices; j++)
				{
					indices.push_back(face.mIndices[j]);
				}
			}

			/* Load material maps
			- Albedo map
			- Specular map
			- Normal map
			*/

			_mesh.m_Name = mesh->mName.C_Str();

			// process materials
			aiMaterial* material = scene->mMaterials[mesh->mMaterialIndex];
			_mesh.m_Color = col;

			_mesh.m_EmissivityColor = emissive_color;

			LoadMaterialTextures(mesh, material, &object->m_Meshes.back(), pth);
		}

		uint32_t mesh_count = 0;

		void ProcessAssimpNode(aiNode* Node, const aiScene* Scene, Object* object, const std::string& pth)
		{
			// Process all the meshes in the node
			// Add the transparent meshes to the transparent mesh queue and add all the opaque ones

			for (int i = 0; i < Node->mNumMeshes; i++)
			{
				mesh_count++;
				aiMesh* mesh = Scene->mMeshes[Node->mMeshes[i]];
				aiMaterial* material = Scene->mMaterials[mesh->mMaterialIndex];
				aiColor4D diffuse_color;
				aiGetMaterialColor(material, AI_MATKEY_COLOR_DIFFUSE, &diffuse_color);

				float transparency = 0.0f;

				if (aiGetMaterialFloat(material, AI_MATKEY_OPACITY, &transparency) == AI_FAILURE)
				{
					transparency = 0.0f;
				}

				aiVector3D _reflectivity;

				if (material->Get(AI_MATKEY_COLOR_REFLECTIVE, _reflectivity) == AI_FAILURE)
				{
					_reflectivity = aiVector3D(0.0f);
				}

				aiVector3D emissivity;

				material->Get(AI_MATKEY_COLOR_EMISSIVE, emissivity);

				

				glm::vec3 reflectivity = glm::vec3(_reflectivity.x, _reflectivity.y, _reflectivity.z);
				glm::vec4 final_color;
				final_color = glm::vec4(diffuse_color.r, diffuse_color.g, diffuse_color.b, diffuse_color.a);
				ProcessAssimpMesh(mesh, Scene, object, pth,
					final_color, reflectivity, glm::vec3(emissivity.x, emissivity.y, emissivity.z));
			}

			for (int i = 0; i < Node->mNumChildren; i++)
			{
				ProcessAssimpNode(Node->mChildren[i], Scene, object, pth);
			}
		}

		void LoadModelFile(Object* object, const std::string& filepath)
		{
			if (filepath.find("glb") != std::string::npos || filepath.find("gltf") != std::string::npos)
			{
				is_gltf = true;
			}

			Assimp::Importer importer;

			const aiScene* Scene = importer.ReadFile
			(
				filepath,
				aiProcess_JoinIdenticalVertices |
				aiProcess_Triangulate |
				aiProcess_CalcTangentSpace |
				aiProcess_GenUVCoords |
				aiProcess_FlipUVs
			);

			if (!Scene || Scene->mFlags & AI_SCENE_FLAGS_INCOMPLETE || !Scene->mRootNode)
			{
				std::stringstream str;
				str << "ERROR LOADING ASSIMP MODEL (" << filepath << ") ||  ASSIMP ERROR : " << importer.GetErrorString();
				Logger::Log(str.str());
				return;
			}

			ProcessAssimpNode(Scene->mRootNode, Scene, object, filepath);

			bool optimize = false;

			if (optimize) {
				PartialOptimize(*object);
			}

			else {
				for (auto& e : object->m_Meshes)
				{
					e.m_AlbedoMap.CreateTexture(e.TexturePaths[0], true, true);
					e.m_NormalMap.CreateTexture(e.TexturePaths[1], false, true);
					e.m_RoughnessMap.CreateTexture(e.TexturePaths[2], false, true);
					e.m_MetalnessMap.CreateTexture(e.TexturePaths[3], false, true);
					e.m_AmbientOcclusionMap.CreateTexture(e.TexturePaths[4], false, true);
					e.m_MetalnessRoughnessMap.CreateTexture(e.TexturePaths[5], false, true);
				}
			}

			object->Buffer();

			mesh_count = 0;
			is_gltf = false;

			return;
		}
	}
}
