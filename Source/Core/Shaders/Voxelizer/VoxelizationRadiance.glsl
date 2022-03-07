#version 450 core 

//layout(RGBA8, binding = 0) uniform image3D o_VoxelVolume;
layout(RGBA16F, binding = 0) uniform image3D o_VoxelVolume;
layout(r16ui, binding = 1) uniform uimage3D o_VoxelNormals;


uniform sampler2D u_AlbedoMap;

in vec3 g_WorldPosition;
in vec3 g_VolumePosition;
in vec3 g_Normal;

in vec2 g_UV;

uniform int u_VolumeSize;
uniform float u_CoverageSizeF;

uniform sampler2D u_Shadowmap;
uniform mat4 u_LightVP;

uniform float u_ModelEmission;
uniform vec3 u_VoxelGridCenterF;

uniform vec3 u_SunDir;

uniform int u_CascadeNumber;

bool InsideCube(vec3 p, float e) { return abs(p.x) < e && abs(p.y) < e && abs(p.z) < e ; } 

float remap(float x, float a, float b, float c, float d)
{
    return (((x - a) / (b - a)) * (d - c)) + c;
}

float GetLuminance(vec3 color) {
	return dot(color, vec3(0.299, 0.587, 0.114));
}

float Luminance(vec3 RGB )
{
    return dot(RGB, vec3(0.2126f, 0.7152f, 0.0722f));
}

const float ReinhardExp = 2.44002939f;
vec3 Reinhard(vec3 RGB )
{
	RGB *= ReinhardExp;
    return vec3(RGB) / (vec3(1.0f) + Luminance(RGB));
}

vec3 InverseReinhard(vec3 RGB)
{
    return (RGB / (vec3(1.0f) - Luminance(RGB))) / ReinhardExp;
}

// By inigo quilez
uint EncodeNormal(in vec3 nor)
{
    vec3 mor; uint  id;
    if( abs(nor.y) > abs(mor.x) ) { mor = nor.yzx; id = 1u; }
    if( abs(nor.z) > abs(mor.x) ) { mor = nor.zxy; id = 2u; }
    uint is = (mor.x<0.0)?1u:0u;
    vec2 uv = 0.5 + 0.5*mor.yz/abs(mor.x);
    uvec2 iuv = uvec2(round(uv*vec2(127.0,63.0)));
    return iuv.x | (iuv.y<<7u) | (id<<13u) | (is<<15u);
}

uint GetEncodedNormal(vec3 Normal) {

	uint Encoded = EncodeNormal(Normal);
	///float FloatBits = uintBitsToFloat(Encoded);
	return Encoded;
}

vec4 EncodeLighting(vec3 Lighting) {
	Lighting /= 5.0f;
	vec3 MappedLighting = Lighting;
	return vec4(MappedLighting, 1.0f);
}

vec3 SunColor = vec3(12.0f);

vec3 SampleLighting(vec3 Albedo, vec3 WorldPosition, vec3 Normal) {

	float Bias = 0.004f;  

	float Shadow = 0.0f; 
	vec4 ProjectionCoordinates = u_LightVP * vec4(WorldPosition, 1.0f);
	ProjectionCoordinates.xyz = ProjectionCoordinates.xyz / ProjectionCoordinates.w;
	ProjectionCoordinates.xyz = ProjectionCoordinates.xyz * 0.5f + 0.5f;

    float ShadowFetch = texture(u_Shadowmap, ProjectionCoordinates.xy).x;
    Shadow = float(ProjectionCoordinates.z - Bias > ShadowFetch);

	float Lambertian = clamp(dot(Normal, -u_SunDir), 0.0f, 1.0f);

    return Albedo * SunColor * Lambertian * (1.0f - Shadow);


}

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

		float Mip = clamp(float(u_CascadeNumber - 1.0f), 2.0f, 6.0f);
		vec3 Albedo = textureLod(u_AlbedoMap, g_UV, Mip).xyz;

		vec3 Direct = SampleLighting(Albedo, g_WorldPosition + u_VoxelGridCenterF, normalize(g_Normal));
		vec3 Emission = Albedo * u_ModelEmission * 2.2f;
		vec3 Combined = Direct + Emission;

		vec3 HDR = Combined;
		vec4 EncodeHDR = EncodeLighting(HDR);

		uint NormalEncoded = GetEncodedNormal(g_Normal);

		imageStore(o_VoxelVolume, VoxelSpaceCoord, EncodeHDR);
		imageStore(o_VoxelNormals, VoxelSpaceCoord, uvec4(NormalEncoded));
	}
}