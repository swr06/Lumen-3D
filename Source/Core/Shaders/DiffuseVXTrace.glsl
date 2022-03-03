#version 400 core 
#extension ARB_bindless_texture : require

#define PI 3.14159265359
#define PHI 1.6180339
#define SAMPLES 1

const float TAU = radians(360.0f);
const float PHI2 = sqrt(5.0f) * 0.5f + 0.5f;
const float GOLDEN_ANGLE = TAU / PHI2 / PHI2;

#define Bayer4(a)   (Bayer2(  0.5 * (a)) * 0.25 + Bayer2(a))
#define Bayer8(a)   (Bayer4(  0.5 * (a)) * 0.25 + Bayer2(a))
#define Bayer16(a)  (Bayer8(  0.5 * (a)) * 0.25 + Bayer2(a))
#define Bayer32(a)  (Bayer16( 0.5 * (a)) * 0.25 + Bayer2(a))
#define Bayer64(a)  (Bayer32( 0.5 * (a)) * 0.25 + Bayer2(a))
#define Bayer128(a) (Bayer64( 0.5 * (a)) * 0.25 + Bayer2(a))
#define Bayer256(a) (Bayer128(0.5 * (a)) * 0.25 + Bayer2(a))

float Bayer2(vec2 a) 
{
    a = floor(a);
    return fract(dot(a, vec2(0.5, a.y * 0.75)));
}

layout (location = 0) out vec4 o_Color;

in vec2 v_TexCoords;

uniform sampler3D u_VoxelVolumes[6];
uniform float u_VoxelRanges[6];
uniform vec3 u_VoxelCenters[6];

uniform sampler2D u_Depth;
uniform sampler2D u_LFNormals;

uniform float u_zNear;
uniform float u_zFar;

uniform mat4 u_Projection;
uniform mat4 u_View;
uniform mat4 u_InverseProjection;
uniform mat4 u_InverseView;

uniform float u_Time;
uniform int u_Frame;

uniform bool u_Checker;

uniform sampler2D u_BlueNoise;

// Environment map
uniform samplerCube u_Skymap;

uniform vec2 u_Jitter;


vec3 WorldPosFromDepth(float depth, vec2 txc)
{
    float z = depth * 2.0 - 1.0;
    vec4 ClipSpacePosition = vec4(txc * 2.0 - 1.0, z, 1.0);
    vec4 ViewSpacePosition = u_InverseProjection * ClipSpacePosition;
    ViewSpacePosition /= ViewSpacePosition.w;
    vec4 WorldPos = u_InverseView * ViewSpacePosition;
    return WorldPos.xyz;
}


bool InsideVolume(vec3 p) { float e = 128.0f; return abs(p.x) < e && abs(p.y) < e && abs(p.z) < e ; } 

bool InScreenspace(vec3 x) {

	if (x.x > 0.0f && x.x < 1.0f && x.y > 0.0f && x.y < 1.0f && x.z > 0.0f && x.z < 1.0f) { 
		return true;
	}

	return false;
}


sampler3D GetCascadeVolume(int cascade) {

	if (cascade == 0) { return u_VoxelVolumes[0]; }
	if (cascade == 1) { return u_VoxelVolumes[1]; }
	if (cascade == 2) { return u_VoxelVolumes[2]; }
	if (cascade == 3) { return u_VoxelVolumes[3]; }
	if (cascade == 4) { return u_VoxelVolumes[4]; }
	if (cascade == 5) { return u_VoxelVolumes[5]; }
	{ return u_VoxelVolumes[0]; }
}

vec3 TransformToVoxelSpace(int Volume, vec3 WorldPosition)
{
	WorldPosition = WorldPosition - u_VoxelCenters[Volume];
	float Size = u_VoxelRanges[Volume];
	float HalfExtent = Size / 2.0f;
	vec3 ScaledPos = WorldPosition / HalfExtent;
	vec3 Voxel = ScaledPos;
	Voxel = Voxel * 0.5f + 0.5f;
	return Voxel;
}

bool PositionInVolume(int Volume, vec3 WorldPosition) {

	WorldPosition = WorldPosition - u_VoxelCenters[Volume];
	float Size = u_VoxelRanges[Volume];
	float HalfExtent = Size / 2.0f;
	vec3 ScaledPos = WorldPosition / HalfExtent;
	vec3 Voxel = ScaledPos;
	Voxel = Voxel * 0.5f + 0.5f;
	return InScreenspace(Voxel);
}

float remap(float x, float a, float b, float c, float d)
{
    return (((x - a) / (b - a)) * (d - c)) + c;
}

float GetLuminance(vec3 color) {
	return dot(color, vec3(0.299, 0.587, 0.114));
}

float Luminance(vec3 x) { return GetLuminance(x) ; }

const float ReinhardExp = 2.44002939f;

vec3 InverseReinhard(vec3 RGB)
{
    return (RGB / (vec3(1.0f) - Luminance(RGB))) / ReinhardExp;
}

vec4 DecodeVolumeLighting(const vec4 Lighting) {

	vec3 RemappedLighting = InverseReinhard(Lighting.xyz);
	RemappedLighting *= 2.0f;
	float AlphaMask = 1.0f - pow(1.0f - Lighting.w, 5.0f);
	return vec4(RemappedLighting.xyz, AlphaMask);
}

int GetCascadeNumber(vec3 P) {

	for (int Cascade = 0 ; Cascade < 6 ; Cascade++) {
		
		if (PositionInVolume(Cascade, P)) {
			return Cascade;
		}

	}

	return 5;
}

vec4 RaymarchCascades(vec3 WorldPosition, vec3 Normal, vec3 Direction, float ConeWidth, float LowDiscrepHash, const int Steps) 
{
	const float Diagonal = sqrt(2.0f);

	float HashScale = mix(1.0f, 2.0f, LowDiscrepHash / 2.0f);

	vec3 RayOrigin = WorldPosition + Normal * Diagonal * 1.0f;

	// Figure out the best cascade to start from
	int CurrentCascade = GetCascadeNumber(RayOrigin);

	float Distance = (u_VoxelRanges[CurrentCascade] / 128.0f) * Diagonal;
	RayOrigin += Direction * Distance;

	vec3 RayPosition = RayOrigin + Direction * Distance * 2.0f * LowDiscrepHash;
	float StepSize = u_VoxelRanges[CurrentCascade] / 128.0f;

	vec4 TotalGI = vec4(0.0f);
	TotalGI.w = 1.0f;

	float SkyVisibility = 1.0f;

	for (int Step = 1 ; Step < Steps ; Step++) {

		if (TotalGI.w < 0.06f) {
			break;
		}

		vec3 Voxel = TransformToVoxelSpace(CurrentCascade, RayPosition);

		sampler3D Volume = GetCascadeVolume(CurrentCascade);

		if (InScreenspace(Voxel)) {

			//float NormalizedStep = float(Step) / float(Steps);
			//float ConeDistanceEstimate = (exp2(NormalizedStep * 4.0f) - 0.9f) / 8.0f;
			//float ConeRadiusSize = ConeDistanceEstimate * 4.0f; 
			//float NaiveMip = clamp(mix(0.0f, 6.0f, NormalizedStep * NormalizedStep));
			float Mip = 0.0f; //ConeRadiusSize;

			vec4 SampleLighting = texture3DLod(Volume, Voxel, Mip).xyzw;
			vec4 Decoded = DecodeVolumeLighting(SampleLighting);
			TotalGI.xyz += Decoded.xyz * Decoded.w;
			TotalGI.w *= 1.0f - SampleLighting.w;
			SkyVisibility *= 1.0f - SampleLighting.w;
		}

		else {

			if (CurrentCascade == 5 || CurrentCascade > 5) {
				break;
			}

			CurrentCascade++;

			StepSize = u_VoxelRanges[CurrentCascade] / 128.0f;
			StepSize *= 1.1f;
		}
		
		RayPosition += Direction * StepSize * HashScale;

	}

	SkyVisibility = pow(SkyVisibility, 64.0f);
	vec3 SkySample = pow(texture(u_Skymap, Direction).xyz, vec3(1.4f));
	vec3 SkyRadiance = SkySample * SkyVisibility * 3.7f;

	float AO = clamp(distance(WorldPosition, RayPosition) / 48.0f, 0.0f, 1.0f);

	return vec4(TotalGI.xyz + SkyRadiance, AO);
}

float HASH2SEED = 0.0f;
vec2 hash2() 
{
	return fract(sin(vec2(HASH2SEED += 0.1, HASH2SEED += 0.1)) * vec2(43758.5453123, 22578.1459123));
}

vec3 CosWeightedHemisphere(const vec3 n, vec2 r) 
{
	float PI2 = 2.0f * 3.1415926f;
	vec3  uu = normalize(cross(n, vec3(0.0,1.0,1.0)));
	vec3  vv = cross(uu, n);
	float ra = sqrt(r.y);
	float rx = ra * cos(PI2 * r.x); 
	float ry = ra * sin(PI2 * r.x);
	float rz = sqrt(1.0 - r.y);
	vec3  rr = vec3(rx * uu + ry * vv + rz * n );
    return normalize(rr);
}


void main() {

	HASH2SEED = (v_TexCoords.x * v_TexCoords.y) * 64.0;
	HASH2SEED += fract(u_Time) * 64.0f;

    float BayerHash = fract(fract(mod(float(u_Frame) + float(0.) * 2., 384.0f) * (1.0 / PHI)) + Bayer32(gl_FragCoord.xy));

	const bool UseBlueNoise = true;

	vec2 Hash; 

	if (UseBlueNoise) {

		int n = u_Frame % 1024;
		vec2 off = fract(vec2(n * 12664745, n * 9560333) / 16777216.0) * 1024.0;
		ivec2 TextureSize = textureSize(u_BlueNoise, 0);
		ivec2 SampleTexelLoc = ivec2(gl_FragCoord.xy + ivec2(floor(off))) % TextureSize;
		Hash = texelFetch(u_BlueNoise, SampleTexelLoc, 0).xy;
	
	}

	else {
		
		Hash = hash2();
	}

	// Calculate pixel 
	ivec2 Pixel = ivec2(gl_FragCoord.xy);

    if (u_Checker) {
        Pixel.x *= 2;
	    bool IsCheckerStep = Pixel.x % 2 == int(Pixel.y % 2 == (u_Frame % 2));
        Pixel.x += int(IsCheckerStep);
    }

    Pixel += ivec2(u_Jitter * 2.0f);
	
    ivec2 HighResPixel = Pixel * 2;
    vec2 HighResUV = vec2(HighResPixel) / textureSize(u_Depth, 0).xy;

	// Fetch 
    float Depth = texelFetch(u_Depth, HighResPixel, 0).x;
	vec3 WorldPosition = WorldPosFromDepth(Depth, HighResUV);
    vec3 Normal = normalize(texelFetch(u_LFNormals, HighResPixel, 0).xyz); 

	vec3 PlayerPosition = u_InverseView[3].xyz;
	vec3 Lo = normalize(PlayerPosition - WorldPosition);
	vec3 R = normalize(reflect(-Lo, Normal));

	vec3 Direction = CosWeightedHemisphere(Normal, Hash.xy);

	vec4 Diffuse = RaymarchCascades(WorldPosition, Normal, Direction, 1.0f, BayerHash, 200);
	o_Color.xyz = Diffuse.xyz;
	o_Color.w = Diffuse.w;
}