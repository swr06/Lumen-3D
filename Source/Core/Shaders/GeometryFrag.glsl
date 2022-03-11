#version 440 core

layout (location = 0) out vec3 o_Albedo;
layout (location = 1) out vec3 o_HFNormal;
layout (location = 2) out vec4 o_PBR;
layout (location = 3) out vec3 o_LFNormal;

uniform sampler2D u_AlbedoMap;
uniform sampler2D u_NormalMap;
uniform sampler2D u_MetalnessMap;
uniform sampler2D u_RoughnessMap;
uniform sampler2D u_MetalnessRoughnessMap;

uniform bool u_UsesGLTFPBR;
uniform bool u_UsesNormalMap;
uniform bool u_UsesRoughnessMap;
uniform bool u_UsesMetalnessMap;
uniform float u_ModelEmission;

uniform vec3 u_EmissiveColor;
uniform bool u_UsesAlbedoTexture;
uniform vec3 u_ModelColor;

uniform float u_EntityRoughness;
uniform float u_EntityMetalness;

in vec2 v_TexCoords;
in vec3 v_FragPosition;
in vec3 v_Normal;
in mat3 v_TBNMatrix;


void main()
{
	const bool Whiteworld = false;

	o_Albedo = Whiteworld ? vec3(1.0f) : (u_UsesAlbedoTexture ? texture(u_AlbedoMap, v_TexCoords).xyz : u_ModelColor);

	o_Albedo += o_Albedo * u_EmissiveColor * u_ModelEmission * 8.0f;

	o_HFNormal = u_UsesNormalMap ? normalize(v_TBNMatrix * (texture(u_NormalMap, v_TexCoords).xyz * 2.0f - 1.0f)) : v_Normal;
	o_LFNormal = v_Normal;

	// https://www.khronos.org/blog/art-pipeline-for-gltf

	if (u_UsesGLTFPBR) {
		//o_PBR.xyz = vec3(texture(u_MetalnessRoughnessMap, v_TexCoords).yx, 1.0f);
		vec4 mapfetch = texture(u_MetalnessRoughnessMap, v_TexCoords);
		o_PBR.xyz = vec3(mapfetch.yz, mapfetch.x);
	}

	else {

		o_PBR.xyz = vec3(u_UsesRoughnessMap ? texture(u_RoughnessMap, v_TexCoords).r : u_EntityRoughness, 
						u_UsesMetalnessMap ? texture(u_MetalnessMap, v_TexCoords).r : u_EntityMetalness, 
						0.0f);

	}

	o_PBR.x = clamp(o_PBR.x, 0.02f, 1.0f);

	o_PBR.w = u_ModelEmission;
}
