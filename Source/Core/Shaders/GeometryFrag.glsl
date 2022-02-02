#version 440 core

layout (location = 0) out vec3 o_Albedo;
layout (location = 1) out vec3 o_HFNormal;
layout (location = 2) out vec3 o_PBR;
layout (location = 3) out vec3 o_LFNormal;

uniform sampler2D u_AlbedoMap;
uniform sampler2D u_NormalMap;
uniform sampler2D u_MetalnessMap;
uniform sampler2D u_RoughnessMap;
uniform sampler2D u_MetalnessRoughnessMap;

uniform bool u_UsesGLTFPBR;

uniform vec3 u_EmissiveColor;

in vec2 v_TexCoords;
in vec3 v_FragPosition;
in vec3 v_Normal;
in mat3 v_TBNMatrix;

void main()
{
	o_Albedo = texture(u_AlbedoMap, v_TexCoords).xyz;


	o_HFNormal = normalize(v_TBNMatrix * (texture(u_NormalMap, v_TexCoords).xyz * 2.0f - 1.0f));
	o_LFNormal = v_Normal;

	if (u_UsesGLTFPBR) {
		o_PBR = vec3(texture(u_MetalnessRoughnessMap, v_TexCoords).yx, 1.0f);
	}

	else {

		o_PBR = vec3(texture(u_RoughnessMap, v_TexCoords).r, 
						texture(u_MetalnessMap, v_TexCoords).r, 
						1.0f);
	}
}
