#version 450 core 

layout(RG8, binding = 0) uniform image3D o_VoxelVolume;

uniform sampler2D u_AlbedoMap;

in vec3 g_WorldPosition;
in vec3 g_VolumePosition;
in vec3 g_Normal;

void main() {

	vec3 TestColor = textureLod(u_AlbedoMap, vec2(0.5f), 6.0f).xyz;

	ivec3 OutputCoord = ivec3(g_VolumePosition);

	imageStore(o_VoxelVolume, OutputCoord, vec4(TestColor, 1.0f));


}