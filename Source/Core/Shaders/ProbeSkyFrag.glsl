#version 440 core

layout (location = 0) out vec4 o_AlbedoRoughness;
layout (location = 1) out float o_Depth;
layout (location = 2) out vec4 o_NormalMetalness;

uniform samplerCube u_EnvironmentMap;
uniform samplerCube u_Mask;

in vec2 v_TexCoords;

uniform mat4 u_InverseProjection;
uniform mat4 u_InverseView;

vec3 GetRayDirectionAt(vec2 txc)
{
	vec4 clip = vec4(txc * 2.0f - 1.0f, -1.0f, 1.0f);
	vec4 eye = vec4(vec2(u_InverseProjection * clip), -1.0f, 0.0f);
	return vec3(u_InverseView * eye);
}

void main()
{
	discard;

	vec3 rd = normalize(GetRayDirectionAt(v_TexCoords));
	
	float Depth = texture(u_Mask, rd).x;

	if (Depth == 0.0f) {
	
		o_AlbedoRoughness.xyz = pow(texture(u_EnvironmentMap, rd).xyz, vec3(2.0f));
		o_AlbedoRoughness.w = 0.0f;
	
		o_Depth = 10000.0f;
	
		o_NormalMetalness = vec4(0.0f);
	}
	
	else {
		discard;
	}

}
