#version 330 core

layout (location = 0) out vec4 o_Albedo;
layout (location = 1) out vec3 o_Normal;

uniform sampler2D u_AlbedoMap;
uniform sampler2D u_NormalMap;

uniform bool u_UsesNormalMap;

uniform bool u_UsesAlbedoTexture;
uniform vec4 u_ModelColor;

in vec2 v_TexCoords;
in vec3 v_FragPosition;
in vec3 v_Normal;
in mat3 v_TBNMatrix;

uniform float u_Time;

float HASH2SEED = 0.0f;
vec2 hash2() 
{
	return fract(sin(vec2(HASH2SEED += 0.1, HASH2SEED += 0.1)) * vec2(43758.5453123, 22578.1459123));
}

void main()
{
	HASH2SEED = (v_TexCoords.x * v_TexCoords.y) * 64.0;
	HASH2SEED += fract(u_Time) * 64.0f;

	vec4 AlbedoFetch = (u_UsesAlbedoTexture ? texture(u_AlbedoMap, v_TexCoords).xyzw : u_ModelColor);

	float Hash = hash2().x;

	float Alpha = AlbedoFetch.w;

	if (Hash > Alpha) {
		discard;
	}

	o_Albedo = vec4(AlbedoFetch.xyz * AlbedoFetch.w, 1.);
}
