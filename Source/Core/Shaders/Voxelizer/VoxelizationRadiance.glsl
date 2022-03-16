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
uniform float u_VoxelRanges[6];

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

uint   packSnorm2x12(vec2 v) { uvec2 d = uvec2(round(2047.5 + v*2047.5)); return d.x|(d.y<<12u); }
uint   packSnorm2x8( vec2 v) { uvec2 d = uvec2(round( 127.5 + v* 127.5)); return d.x|(d.y<< 8u); }
vec2 unpackSnorm2x8( uint d) { return vec2(uvec2(d,d>> 8)& 255u)/ 127.5 - 1.0; }
vec2 unpackSnorm2x12(uint d) { return vec2(uvec2(d,d>>12)&4095u)/2047.5 - 1.0; }

// By inigo quilez

uint EncodeNormal(in vec3 nor)
{
    vec2 v = vec2( atan(nor.z,nor.x)/3.141593, -1.0+2.0*acos(nor.y)/3.141593 );
    return packSnorm2x8(v);
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

vec3 SunColor = vec3(13.5f);

vec3 SampleLighting(vec3 Albedo, vec3 WorldPosition, vec3 Normal) {

	const vec2 Poisson[6] = vec2[6](vec2(-0.613392, 0.617481),  vec2(0.751946, 0.453352),
                                vec2(0.170019, -0.040254),  vec2(0.078707, -0.715323),
                                vec2(-0.299417, 0.791925),  vec2(-0.075838, -0.529344));

	//WorldPosition += Normal * (u_VoxelRanges[u_CascadeNumber] * (1.0f/128.0f) * sqrt(2.));

	float Bias = 0.004f;  

	float Shadow = 0.0f; 
	vec4 ProjectionCoordinates = u_LightVP * vec4(WorldPosition, 1.0f);
	ProjectionCoordinates.xyz = ProjectionCoordinates.xyz / ProjectionCoordinates.w;
	ProjectionCoordinates.xyz = ProjectionCoordinates.xyz * 0.5f + 0.5f;
        
    float Depth = ProjectionCoordinates.z;

    if (true) {
           
        vec2 TexelSize = 1.0f / textureSize(u_Shadowmap, 0).xy;

		int Steps = 5;

        // pcf 
        for (int i = 0 ; i < Steps ; i++) {
            float Fetch = texture(u_Shadowmap, ProjectionCoordinates.xy + Poisson[i] * TexelSize * 2.4f).x;
            Shadow += float(ProjectionCoordinates.z - Bias > Fetch);
        }

        Shadow /= float(Steps);

    }

	float Lambertian = clamp(dot(Normal, -u_SunDir), 0.0f, 1.0f);

    return Albedo * SunColor * 1.3f * Lambertian * (1.0f - Shadow);
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

		vec3 Normal = normalize(g_Normal);

		ivec3 VoxelSpaceCoord = ivec3(Voxel * float(u_VolumeSize));

		float Mip = clamp(float(u_CascadeNumber - 1.0f), 2.0f, 5.0f);
		vec3 RawAlbedo = textureLod(u_AlbedoMap, g_UV, Mip).xyz;
		vec3 Albedo = RawAlbedo;
		Albedo = pow(Albedo, vec3(1.7f)) * 1.8f;

		vec3 Direct = SampleLighting(Albedo, g_WorldPosition + u_VoxelGridCenterF, Normal);
		vec3 Emission = RawAlbedo * u_ModelEmission * 2.2f;
		vec3 Combined = Direct + Emission;

		vec3 HDR = Combined + vec3(Albedo * 0.01f);
		vec4 EncodeHDR = EncodeLighting(HDR);

		uint NormalEncoded = GetEncodedNormal(Normal);

		imageStore(o_VoxelVolume, VoxelSpaceCoord, EncodeHDR);
		imageStore(o_VoxelNormals, VoxelSpaceCoord, uvec4(NormalEncoded));
	}
}