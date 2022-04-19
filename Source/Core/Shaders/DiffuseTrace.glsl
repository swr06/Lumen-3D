#version 430 core 
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

// Voxel cascade info
uniform sampler3D u_VoxelVolumes[6];
uniform usampler3D u_VoxelVolumesNormals[6];
uniform float u_VoxelRanges[6];
uniform vec3 u_VoxelCenters[6];

// GBuffer
uniform sampler2D u_Depth;
uniform sampler2D u_LFNormals;
uniform sampler2D u_LowResDepth;
uniform sampler2D u_Albedo;
uniform sampler2D u_PBR;


uniform float u_zNear;
uniform float u_zFar;
uniform mat4 u_Projection;
uniform mat4 u_View;
uniform mat4 u_InverseProjection;
uniform mat4 u_InverseView;
uniform mat4 u_ViewProjection;

uniform float u_Time;
uniform int u_Frame;
uniform bool u_Checker;
uniform bool u_SSRTGI;
uniform bool u_HQ_SSRTGI;
uniform bool u_AmplifySunGI;

uniform sampler2D u_BlueNoise;

uniform bool u_DiffuseSecondBounce;

// Environment map
uniform samplerCube u_Skymap;
uniform vec2 u_Jitter;

// Camera-centric probe data 
uniform vec3 u_ProbeCapturePoints[6];
uniform samplerCube u_ProbeAlbedo;
uniform samplerCube u_ProbeDepth;
uniform samplerCube u_ProbeNormals;

// Direct 
uniform sampler2D u_Shadowmap;
uniform mat4 u_SunShadowMatrix;
uniform vec3 u_SunDirection;

// Pseudo whitenoise rng 
float HASH2SEED = 0.0f;
vec2 hash2() 
{
	return fract(sin(vec2(HASH2SEED += 0.1, HASH2SEED += 0.1)) * vec2(43758.5453123, 22578.1459123));
}

// Math functions 
float remap(float x, float a, float b, float c, float d)
{
    return (((x - a) / (b - a)) * (d - c)) + c;
}

float GetLuminance(vec3 color) {
	return dot(color, vec3(0.299, 0.587, 0.114));
}

float Luminance(vec3 x) { return GetLuminance(x) ; } // wrapper 

mat4 SaturationMatrix(float saturation)
{
  vec3 luminance = vec3(0.3086, 0.6094, 0.0820);
  float oneMinusSat = 1.0 - saturation;
  vec3 red = vec3(luminance.x * oneMinusSat);
  red += vec3(saturation, 0, 0);
  vec3 green = vec3(luminance.y * oneMinusSat);
  green += vec3(0, saturation, 0);
  vec3 blue = vec3(luminance.z * oneMinusSat);
  blue += vec3(0, 0, saturation);
  return mat4(red, 0,
    green, 0,
    blue, 0,
    0, 0, 0, 1);
}

vec3 ProjectToClipSpace(vec3 WorldPos) 
{
	vec4 ProjectedPosition = u_ViewProjection * vec4(WorldPos, 1.0f);
	ProjectedPosition.xyz /= ProjectedPosition.w;
	return ProjectedPosition.xyz;
}

bool SSRayValid(vec2 x) {
	float bias = 0.0001f;
	if (x.x > bias && x.x < 1.0f - bias && x.y > bias && x.y < 1.0f - bias) {
		return true;
	}

	return false;
}

float LinearizeDepth(float depth)
{
	return (2.0 * u_zNear) / (u_zFar + u_zNear - depth * (u_zFar - u_zNear));
}

vec3 ProjectToScreenSpace(vec3 WorldPos) 
{
	vec4 ProjectedPosition = u_ViewProjection * vec4(WorldPos, 1.0f);
	ProjectedPosition.xyz /= ProjectedPosition.w;
	ProjectedPosition.xyz = ProjectedPosition.xyz * 0.5f + 0.5f;
	return ProjectedPosition.xyz;
}


// Space conversions ->
vec3 WorldPosFromDepth(float depth, vec2 txc)
{
    float z = depth * 2.0 - 1.0;
    vec4 ClipSpacePosition = vec4(txc * 2.0 - 1.0, z, 1.0);
    vec4 ViewSpacePosition = u_InverseProjection * ClipSpacePosition;
    ViewSpacePosition /= ViewSpacePosition.w;
    vec4 WorldPos = u_InverseView * ViewSpacePosition;
    return WorldPos.xyz;
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

bool LiesInsideVolume(vec3 vp) { return abs(vp.x) < 128.0f && abs(vp.y) < 128.0f  && abs(vp.z) < 128.0f ; } 

bool InScreenspace(vec3 x) {

	if (x.x > 0.0f && x.x < 1.0f && x.y > 0.0f && x.y < 1.0f && x.z > 0.0f && x.z < 1.0f) { 
		return true;
	}

	return false;
}

const vec3 FACE_NORMALS[6] = vec3[](vec3(1.0, 0.0, 0.0),vec3(-1.0, 0.0, 0.0),vec3(0.0, 1.0, 0.0),vec3(0.0, -1.0, 0.0),vec3(0.0, 0.0, 1.0),vec3(0.0, 0.0, -1.0));
const vec3 OTHER_NORMALS[6] = vec3[](vec3(1.0, 0.0, 0.0),vec3(0.0, 1.0, 0.0),vec3(0.0, 0.0, 1.0),vec3(0.0, -1.0, 0.0),vec3(0.0, 0.0, 1.0),vec3(0.0, 0.0, -1.0));

// MOTHERFUCKINGGG issue where you can't reference uniform texture arrays with a variable index (so this ugly fucking hack needed to be implemented)
sampler3D GetCascadeVolume(int cascade) {

	if (cascade == 0) { return u_VoxelVolumes[0]; }
	if (cascade == 1) { return u_VoxelVolumes[1]; }
	if (cascade == 2) { return u_VoxelVolumes[2]; }
	if (cascade == 3) { return u_VoxelVolumes[3]; }
	if (cascade == 4) { return u_VoxelVolumes[4]; }
	if (cascade == 5) { return u_VoxelVolumes[5]; }
	{ return u_VoxelVolumes[0]; }
}

usampler3D GetCascadeVolumeN(int cascade) {

	if (cascade == 0) { return u_VoxelVolumesNormals[0]; }
	if (cascade == 1) { return u_VoxelVolumesNormals[1]; }
	if (cascade == 2) { return u_VoxelVolumesNormals[2]; }
	if (cascade == 3) { return u_VoxelVolumesNormals[3]; }
	if (cascade == 4) { return u_VoxelVolumesNormals[4]; }
	if (cascade == 5) { return u_VoxelVolumesNormals[5]; }
	{ return u_VoxelVolumesNormals[0]; }
}

// Returns whether a world position P lies inside a cascade
bool PositionInVolume(int Volume, vec3 WorldPosition) {

	WorldPosition = WorldPosition - u_VoxelCenters[Volume];
	float Size = u_VoxelRanges[Volume];
	float HalfExtent = Size / 2.0f;
	vec3 ScaledPos = WorldPosition / HalfExtent;
	vec3 Voxel = ScaledPos;
	Voxel = Voxel * 0.5f + 0.5f;
	return InScreenspace(Voxel);
}

const float ReinhardExp = 1.0f;

vec3 InverseReinhard(vec3 RGB){
    return (RGB / (vec3(1.0f) - Luminance(RGB))) / ReinhardExp;
}

vec4 DecodeVolumeLighting(const vec4 Lighting) {
	return vec4(Lighting.xyz * 5.0f, Lighting.w);
}

uint   packSnorm2x12(vec2 v) { uvec2 d = uvec2(round(2047.5 + v*2047.5)); return d.x|(d.y<<12u); }
uint   packSnorm2x8( vec2 v) { uvec2 d = uvec2(round( 127.5 + v* 127.5)); return d.x|(d.y<< 8u); }
vec2 unpackSnorm2x8( uint d) { return vec2(uvec2(d,d>> 8)& 255u)/ 127.5 - 1.0; }
vec2 unpackSnorm2x12(uint d) { return vec2(uvec2(d,d>>12)&4095u)/2047.5 - 1.0; }

vec3 DecodeNormal(uint data)
{
    vec2 v = unpackSnorm2x8(data);
    v.y = 0.5+0.5*v.y; v *= 3.141593;
    return normalize( vec3( sin(v.y)*cos(v.x), cos(v.y), sin(v.y)*sin(v.x) ));
}


// Voxel raytracing ->

const bool CONE_TRACE = false;

int GetCascadeNumber(vec3 P, int MinCascade) {

	for (int Cascade = int(MinCascade) ; Cascade < 6 ; Cascade++) {
		
		if (PositionInVolume(Cascade, P)) {
			return Cascade;
		}

	}

	return 5;
}

// Voxel volume raytracer 
vec4 RaymarchCascades(vec3 WorldPosition, vec3 Normal, vec3 Direction, float Aperature, float LowDiscrepHash, const int Steps, out vec4 Intersection, int MinCascade, out vec3 VoxelNormal, bool CalcNormal) 
{
	const float Diagonal = sqrt(2.0f);

	float HashScale = mix(1.0f, 2.0f, LowDiscrepHash / 2.0f);

	vec3 RayOrigin = WorldPosition + Normal * Diagonal * 1.0f;

	// Figure out the best cascade to start from
	int CurrentCascade = GetCascadeNumber(RayOrigin, MinCascade);

	float Distance = (u_VoxelRanges[CurrentCascade] / 128.0f) * Diagonal;
	RayOrigin += Direction * Distance;

	vec3 RayPosition = RayOrigin + Direction * Distance * 2.0f * LowDiscrepHash;
	float StepSize = u_VoxelRanges[CurrentCascade] / 128.0f;

	vec4 TotalGI = vec4(0.0f);
	TotalGI.w = 1.0f;

	float SkyVisibility = 1.0f;

	bool IntersectionFound = true;

	VoxelNormal = vec3(0.0f, 0.0f, 1.0f);
	vec3 VoxelPosition = vec3(0.0f);

	for (int Step = 1 ; Step < Steps ; Step++) {

		if (TotalGI.w < 0.06f) {
			break;
		}

		VoxelPosition = TransformToVoxelSpace(CurrentCascade, RayPosition);

		VoxelPosition = floor(VoxelPosition * 128.0f) / 128.0f;

		sampler3D Volume = GetCascadeVolume(CurrentCascade);

		if (InScreenspace(VoxelPosition)) {

			vec4 SampleLighting = texture3DLod(Volume, VoxelPosition, 0.).xyzw;
			vec4 Decoded = DecodeVolumeLighting(SampleLighting);
			 
			TotalGI.xyz += Decoded.xyz * Decoded.w;
			TotalGI.w *= 1.0f - SampleLighting.w;
			SkyVisibility *= 1.0f - SampleLighting.w;
		}

		else {

			if (CurrentCascade == 5 || CurrentCascade > 5) {
				
				IntersectionFound = false;

				break;
			}

			CurrentCascade++;

			StepSize = u_VoxelRanges[CurrentCascade] / 128.0f;
		}
		
		RayPosition += Direction * StepSize * HashScale;
	}

	// Normal calculation 
	uint EncodedNormal = texelFetch(GetCascadeVolumeN(CurrentCascade), ivec3(floor(VoxelPosition * 128.0f)), 0).x;
	VoxelNormal = DecodeNormal(EncodedNormal);

	SkyVisibility = pow(SkyVisibility, 24.0f);
	vec3 SkySample = pow(texture(u_Skymap, Direction).xyz, vec3(2.0f));
	float LambertSky = pow(clamp(dot(VoxelNormal, vec3(0.0f, 1.0f, 0.0f)), 0.0f, 1.0f), 1.0f);
	vec3 SkyRadiance = SkySample * 1.4f * LambertSky * SkyVisibility;

	float AO = clamp(distance(WorldPosition, RayPosition) / 52.0f, 0.0f, 1.0f);

	AO = pow(AO, 1.5f);

	Intersection = vec4(RayPosition, float(IntersectionFound));

	return vec4(SkyRadiance + TotalGI.xyz, AO);
}



// Samples direct lighting for a point 
vec3 SampleRadiance(vec3 P, vec3 N, vec3 A, vec3 E, bool AmplifySunGI) {
	const mat4 Satmat = SaturationMatrix(2.2f);
	const vec3 SUN_COLOR = u_AmplifySunGI ? vec3(10.0f) : vec3(8.0f);
	vec4 ProjectionCoordinates = u_SunShadowMatrix * vec4(P + N * 0.5f, 1.0f);
	ProjectionCoordinates.xyz = ProjectionCoordinates.xyz / ProjectionCoordinates.w;
    ProjectionCoordinates.xyz = ProjectionCoordinates.xyz * 0.5f + 0.5f;
	float SimpleFetch = texture(u_Shadowmap, ProjectionCoordinates.xy).x;
    float Shadow = float(ProjectionCoordinates.z - (0.0045f) > SimpleFetch);
	float Lambertian = clamp(max(0.0f, dot(N, -u_SunDirection))*1.5f,0.0f,1.0f); 
    vec3 Direct = Lambertian * SUN_COLOR * (1.0f - Shadow) * A;
	Direct = AmplifySunGI ? vec3(Satmat * vec4(Direct, 1.0f)) : Direct;
	return Direct + E;
}




// Camera centric probe raytracing -> -> //

// Gets the ID of the face a direction faces in a virtual unit cube (direction has to be normalized)
int GetFaceID(vec3 Direction)
{
    vec3 AbsoluteDirection = abs(Direction);
    float Index = 0.0f;
	if(AbsoluteDirection.z >= AbsoluteDirection.x && AbsoluteDirection.z >= AbsoluteDirection.y){
		Index = Direction.z < 0.0 ? 5.0 : 4.0;
	}

	else if(AbsoluteDirection.y >= AbsoluteDirection.x){
		Index = Direction.y < 0.0 ? 3.0 : 2.0;
	}

	else{
		Index = Direction.x < 0.0 ? 1.0 : 0.0;
	}

    return int(Index);
}

// Probe (per face) capture position 
vec3 GetCapturePoint(vec3 Direction) {
    return u_ProbeCapturePoints[clamp(GetFaceID(Direction),0,5)];
}

// Raytraces the camera centric probe 
vec4 RaytraceProbe(const vec3 WorldPosition, vec3 Direction, float Hash, int Steps, int BinarySteps, float ErrorTolerance, out vec3 oNormal) 
{
	const float Distance = 384.0f;

    float StepSize = Distance / float(Steps);

    vec3 ReflectionVector = Direction; 
    
    vec3 RayPosition = WorldPosition + ReflectionVector * Hash;
    vec3 RayOrigin = RayPosition;
    vec3 PreviousSampleDirection = ReflectionVector;
    bool FoundHit = false; 
    float ExpStep = 1.03f;

    float PrevRayDepth = 0.0f;
    float PrevStepSampleDepth = 0.0f;

    // Find intersection with geometry 
    // Todo : Account for geometrical thickness?
    for (int CurrentStep = 0; CurrentStep < Steps ; CurrentStep++) 
    {
        if (CurrentStep > Steps / 3) {
            StepSize *= ExpStep;
        }

        vec3 CapturePoint = GetCapturePoint(PreviousSampleDirection);
        vec3 Diff = RayPosition - CapturePoint;
        float L = length(Diff);
        vec3 SampleDirection = Diff / L;
        PreviousSampleDirection = SampleDirection;
        float ProbeDepth = (texture(u_ProbeDepth, SampleDirection).x * 128.0f);
        float DepthError = abs(ProbeDepth - L);
        float ThresholdCurr = abs(L - PrevRayDepth); 
        ThresholdCurr = (clamp(ThresholdCurr * 7.0f, 0.01f, Distance) * ErrorTolerance);

        bool AccurateishHit = DepthError < ThresholdCurr;

		// We found a hit! break and integrate lighting 
        if (L > ProbeDepth && AccurateishHit) {
             FoundHit = true;
             break;
        }

		// Break if ray is behind depth buffer 
		if (L > ProbeDepth) { break; }

        PrevRayDepth = L;
        PrevStepSampleDepth = ProbeDepth;
        RayPosition += ReflectionVector * StepSize;
        
    }

    if (FoundHit) 
    {
        vec3 FinalBinaryRefinePos = RayPosition;

        float BR_StepSize = StepSize / 2.0f;
        FinalBinaryRefinePos = FinalBinaryRefinePos - ReflectionVector * BR_StepSize;

		float LlastDepth = 0.0f;

        for (int BinaryRefine = 0 ; BinaryRefine < BinarySteps; BinaryRefine++) 
        {
            BR_StepSize /= 2.0f;
            vec3 CapturePoint = GetCapturePoint(PreviousSampleDirection);
            vec3 Diff = FinalBinaryRefinePos - CapturePoint;
            float L = length(Diff);
            vec3 BR_SampleDirection = Diff / L;
            PreviousSampleDirection = BR_SampleDirection;
            float Depth = (texture(u_ProbeDepth, BR_SampleDirection).x * 128.0f); LlastDepth = Depth;
            float RaySign = (Depth < L) ? -1.0f : 1.0f;
            FinalBinaryRefinePos += ReflectionVector * BR_StepSize * RaySign;
        }

		// Sample data at intersection point 
        vec3 CapturePoint = GetCapturePoint(PreviousSampleDirection);
        vec3 FinalSampleDirection = PreviousSampleDirection;
        float DepthFetch = LlastDepth;
        vec4 AlbedoFetch = textureLod(u_ProbeAlbedo, FinalSampleDirection, 0.0f);
        vec4 NormalFetch = textureLod(u_ProbeNormals, FinalSampleDirection, 0.0f);

		// Integrate lighting 
        vec3 IntersectionPos = (CapturePoint + DepthFetch * FinalSampleDirection);
        vec3 IntersectNormal = NormalFetch.xyz;
        vec3 IntersectAlbedo = AlbedoFetch.xyz;
		vec3 IntersectEmissive = mix(AlbedoFetch.xyz,vec3(Luminance(AlbedoFetch.xyz)),0.1f) * NormalFetch.w * 16.0f;
		oNormal = IntersectNormal.xyz;
		return vec4(SampleRadiance(IntersectionPos, IntersectNormal, IntersectAlbedo, IntersectEmissive, u_AmplifySunGI),
					clamp(distance(WorldPosition, RayPosition) / 42.0f, 0.0f, 1.0f)) ;
    }

	oNormal = vec3(-1.);
    return vec4(vec3(0.0f), 1.0f);
}


// Screenspace raytracer 
vec4 ScreenspaceRaytrace(const vec3 Origin, vec3 Direction, float ThresholdMultiplier, float Hash, int Steps, int BinarySteps, out vec4 IntersectionPos, out vec3 IntersectionNormal)
{
	IntersectionPos = vec4(-1.0f);

    const float Distance = u_HQ_SSRTGI ? 192.0f : 128.0f;

    float StepSize = float(Distance) / float(Steps);
	vec3 RayPosition = Origin + Direction * Hash * 1.2; 
	vec2 FinalUV = vec2(-1.0f);

	//float ExpStep = clamp((log(500.0f) / log(float(Steps))) * 0.75f, 1.05f, 6.0f);// mix(1.075f, 1.5f, float(Hash));
    float ExpStep = u_HQ_SSRTGI ? mix(1.02f, 1.05f, clamp(hash2().x * 1.33333f, 0.0f, 1.0f)) : mix(1.03f, 1.065f, clamp(hash2().x * 1.33333f, 0.0f, 1.0f)); 

	for(int CurrentStep = 0; CurrentStep < Steps; CurrentStep++) 
	{
		float Threshold = StepSize * ThresholdMultiplier * 2.0f;
		
		vec3 ProjectedRayScreenspace = ProjectToClipSpace(RayPosition); 
		
		if(abs(ProjectedRayScreenspace.x) > 1.0f || abs(ProjectedRayScreenspace.y) > 1.0f || abs(ProjectedRayScreenspace.z) > 1.0f) 
		{
			return vec4(vec3(0.0f), -1.0f);
		}
		
		ProjectedRayScreenspace.xyz = ProjectedRayScreenspace.xyz * 0.5f + 0.5f; 

		if (!SSRayValid(ProjectedRayScreenspace.xy))
		{
			return vec4(vec3(0.0f), -1.0f);
		}
		
		float DepthAt = texture(u_LowResDepth, ProjectedRayScreenspace.xy).x; 
		float CurrentRayDepth = LinearizeDepth(ProjectedRayScreenspace.z); 
		float Error = abs(LinearizeDepth(DepthAt) - CurrentRayDepth);
		
        // Intersected!
		if (Error < Threshold && ProjectedRayScreenspace.z > DepthAt) 
		{
			// Binary search for best intersection point along ray step 
            bool DoBinaryRefinement = true;

            vec3 FinalProjected = vec3(0.0f);
            float FinalDepth = 0.0f;

            if (DoBinaryRefinement) {
			    vec3 BinaryStepVector = (Direction * StepSize) / 2.0f;
                RayPosition -= (Direction * StepSize) / 2.0f;
			    
                for (int BinaryStep = 0 ; BinaryStep < BinarySteps ; BinaryStep++) {
			    		
			    	BinaryStepVector /= 2.0f;
			    	vec3 Projected = ProjectToClipSpace(RayPosition); 
			    	Projected = Projected * 0.5f + 0.5f;
                    FinalProjected = Projected;
                    float Fetch = texture(u_LowResDepth, Projected.xy).x;
                    FinalDepth = Fetch;
			    	float BinaryDepthAt = LinearizeDepth(Fetch); 
			    	float BinaryRayDepth = LinearizeDepth(Projected.z); 

			    	if (BinaryDepthAt < BinaryRayDepth) {
			    		RayPosition -= BinaryStepVector;

			    	}

			    	else {
			    		RayPosition += BinaryStepVector;

			    	}

			    }
            }

            else {
                
                FinalProjected = ProjectToScreenSpace(RayPosition);
                FinalDepth = texture(u_Depth, FinalProjected.xy).x;
            }


            // Generate gbuffer data and return 
            vec3 Position = WorldPosFromDepth(FinalDepth, FinalProjected.xy);
            vec3 Normal = texture(u_LFNormals, FinalProjected.xy).xyz;
            vec3 Albedo = texture(u_Albedo, FinalProjected.xy).xyz;
            vec3 Emissive = mix(Albedo, vec3(Luminance(Albedo)), 0.2f) * texture(u_PBR, FinalProjected.xy).w * 24.0f;

			IntersectionPos = vec4(Position, 1.0f);
			IntersectionNormal = Normal;

			//return vec4(Albedo, 1 + 1.0f);

			float AO = clamp(distance(IntersectionPos.xyz, Origin) / 52.0f, 0.0f, 1.0f);
			AO = pow(AO, 1.5f);
			return vec4(SampleRadiance(Position, Normal, Albedo, Emissive, u_AmplifySunGI) * 2.5f, AO + 1.0f);
		}

		if (ProjectedRayScreenspace.z > DepthAt) {
			return vec4(vec3(0.0f), -1.0f);
		}

        // Step 
		RayPosition += StepSize * Direction; 

        StepSize *= ExpStep;
	}

	return vec4(vec3(0.0f), -1.0f);

}


// Lambertian ray generation
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

// Lambert brdf function 
float LambertianBRDF(float N, float D) {
	return max(dot(N, D), 0.0f) / PI;
}

float InverseSchlick(float f0, float VoH) {
    return 1.0 - clamp(f0 + (1.0f - f0) * pow(1.0f - VoH, 5.0f), 0.0f, 1.0f);
}

// Hammon diffuse BRDF 
float DiffuseHammon(vec3 normal, vec3 viewDir, vec3 lightDir, float roughness)
{
    float nDotL = max(dot(normal, lightDir), 0.0f);

    if (nDotL <= 0.0) 
	{
		return 0.0f; 
	}

	const float InvPI = 1.0f/PI;

    float nDotV = max(dot(normal, viewDir), 0.0f);
    float lDotV = max(dot(lightDir, viewDir), 0.0f);
    vec3 halfWay = normalize(viewDir + lightDir);
    float nDotH = max(dot(normal, halfWay), 0.0f);
    float facing = lDotV * 0.5f + 0.5f;
    float singleRough = facing * (0.9f - 0.4f * facing) * ((0.5f + nDotH) * rcp(max(nDotH, 0.02)));
    float singleSmooth = 1.05f * InverseSchlick(0.0f, nDotL) * InverseSchlick(0.0f, max(nDotV, 0.0f));
    float single = clamp(mix(singleSmooth, singleRough, roughness) * rcp(PI), 0.0f, 1.0f);
    float multi = 0.1159f * roughness;
    return clamp((multi + single) * nDotL, 0.0f, 1.0f); // approximate for multi scattering as well
}

vec3 RayBRDFHammon(vec3 N, vec3 I, vec3 D, float Roughness, vec3 Lighting)
{
	vec3 Albedo = vec3(0.5f); 
    vec3 BRDF = Albedo * DiffuseHammon(N, -I, D, Roughness);
    return BRDF;
}

bool IsSky(float NonLinearDepth) {
    if (NonLinearDepth > 0.9999992f || NonLinearDepth == 1.0f) {
        return true;
	}

    return false;
}

// Trace settings 

// First bounce settings 
const bool VX_FIRST_BOUNCE = true;

// Second bounce settings 
const bool DO_SECOND_BOUNCE = true;
const bool VX_SECOND_BOUNCE = true;



void main() {

	HASH2SEED = (v_TexCoords.x * v_TexCoords.y) * 64.0;
	HASH2SEED += fract(u_Time) * 64.0f;

	// Animate by golden ratio 
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
	
	// 1/2 res on each axis
    ivec2 HighResPixel = Pixel * 2;
    vec2 HighResUV = vec2(HighResPixel) / textureSize(u_Depth, 0).xy;

	// Fetch 
    float Depth = texelFetch(u_Depth, HighResPixel, 0).x;

	 if (IsSky(Depth)) {
        o_Color = vec4(vec3(1.0f), 1.0f);
        return;
    }

	vec3 WorldPosition = WorldPosFromDepth(Depth, HighResUV);
    vec3 Normal = normalize(texelFetch(u_LFNormals, HighResPixel, 0).xyz); 

	vec3 PlayerPosition = u_InverseView[3].xyz;
	vec3 Incident = normalize(PlayerPosition - WorldPosition);
	vec3 R = normalize(reflect(-Incident, Normal));

	// Raytrace and accumulate radians ->

	vec4 TotalRadiance = vec4(0.0f);
	vec3 CosineHemisphereDirection = CosWeightedHemisphere(Normal, Hash.xy);

	vec4 IntersectionPositionb = vec4(0.0f);
	vec3 IntersectionNormalb = vec3(0.0f);
	vec4 DiffuseVX = vec4(vec3(0.0f), 1.0f);
	
	// Voxel trace the first bounce? (enabled by default)
	if (VX_FIRST_BOUNCE) {

		// Raytrace in screen space?
		if (u_SSRTGI) {
			vec4 SSRadiance = ScreenspaceRaytrace(WorldPosition + Normal * 1.5, CosineHemisphereDirection, 0.0015f, BayerHash, u_HQ_SSRTGI ? 64 : 32, 6, IntersectionPositionb, IntersectionNormalb);	
		
			// Did the screen trace find an intersection? If so, add lighting for that point
			if (SSRadiance.w > 0.0f) {
				DiffuseVX = vec4(SSRadiance.xyz, SSRadiance.w - 1.0f);
			}

			// RT in voxel space
			else {
				DiffuseVX = RaymarchCascades(WorldPosition, Normal, CosineHemisphereDirection, 1.0f, BayerHash, 236, IntersectionPositionb, 1, IntersectionNormalb, true);
			}

		}

		// RT in voxel space
		else {
			DiffuseVX = RaymarchCascades(WorldPosition, Normal, CosineHemisphereDirection, 1.0f, BayerHash, 256, IntersectionPositionb, 1, IntersectionNormalb, true);
		}
	}

	else {
		DiffuseVX = RaytraceProbe(WorldPosition.xyz + Normal * 3.f, CosineHemisphereDirection, BayerHash, 64, 32, 0.625f, IntersectionNormalb) * 2.;
	}


	TotalRadiance += DiffuseVX;

	Hash = hash2();

	if (DO_SECOND_BOUNCE && IntersectionPositionb.w > 0.6f && u_DiffuseSecondBounce) {
		
		Hash = hash2();
		vec3 SecondCosineHemisphereDirection = CosWeightedHemisphere(IntersectionNormalb, Hash.xy);

		// Calculate ray throughput 
		float CosinePDF = clamp(dot(IntersectionNormalb.xyz, SecondCosineHemisphereDirection) / PI, 0.000001f, 1.0f);

		// I use Hammon brdf as the ray brdf. We divide the ray brdf by ray sampling brdf (lambert) to get the bounce ray weight 
		// this actually isn't perfectly correct, I assume a constant albedo which is wrong. But it's not feasible to store the albedo since we have limited memory 
		vec3 SecondBounceWeight = clamp(vec3(1.0f) * (RayBRDFHammon(IntersectionNormalb, CosineHemisphereDirection, SecondCosineHemisphereDirection, 0.99f, TotalRadiance.xyz) / CosinePDF), 0.0f, PI); 

		vec3 SecondaryBounceRadiance = vec3(0.0f); 

		if (u_DiffuseSecondBounce) {
			if (!VX_SECOND_BOUNCE) {

				vec3 In;

				// Raytrace probe ->
				vec4 ProbeTrace = RaytraceProbe(IntersectionPositionb.xyz + IntersectionNormalb * 2.f, SecondCosineHemisphereDirection, BayerHash, 32, 12, 0.5f, In);
				
				// Is the probe hit valid?
				if (ProbeTrace.w > 0.5f) {

					SecondaryBounceRadiance.xyz = ProbeTrace.xyz;

				}
			}

			else {
				
				// Raytrace voxel volume ->
				
				vec4 oPos = vec4(0.0f);
				vec3 Nn = vec3(0.0f);

				vec4 SecondBounceVX = RaymarchCascades(IntersectionPositionb.xyz, IntersectionNormalb, SecondCosineHemisphereDirection, 1.0f, BayerHash, 92, oPos, 3, Nn, false);

				SecondaryBounceRadiance.xyz = SecondBounceVX.xyz;

			}
		}

		// Account for secondary bounce 

		TotalRadiance.xyz += SecondaryBounceRadiance * SecondBounceWeight;
	}

	o_Color = TotalRadiance;
    o_Color.xyz = max(o_Color.xyz, 0.0f);

	if (isnan(o_Color.x) || isnan(o_Color.y) || isnan(o_Color.z) || isnan(o_Color.w) || isinf(o_Color.x) || isinf(o_Color.y) || isinf(o_Color.z) || isinf(o_Color.w)) {
        o_Color.xyzw = vec4(0.0f);
    }
} 