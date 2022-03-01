#version 450 core 

layout(RG8, binding = 0) uniform image3D o_VoxelVolume;

uniform sampler2D u_AlbedoMap;

in vec3 g_WorldPosition;
in vec3 g_VolumePosition;
in vec3 g_Normal;

in vec2 g_UV;

uniform int u_VolumeSize;
uniform float u_CoverageSizeF;

bool InsideCube(vec3 p, float e) { return abs(p.x) < e && abs(p.y) < e && abs(p.z) < e ; } 

void main() {

	float Size = u_CoverageSizeF;
	float HalfExtent = Size / 2.0f;

	vec3 Clipspace = g_WorldPosition / HalfExtent;

	// Clip space ->
	vec3 Voxel = Clipspace;

	// Clip -> Screen space
	Voxel = Voxel * 0.5f + 0.5f;

	if (Voxel == clamp(Voxel, 0.0f, 1.0f)) {
		ivec3 VoxelSpaceCoord = ivec3(Voxel * float(u_VolumeSize));
		vec3 Radiance = textureLod(u_AlbedoMap, g_UV, 2.0f).xyz;
		imageStore(o_VoxelVolume, VoxelSpaceCoord, vec4(Radiance, 1.0f));
	}
}