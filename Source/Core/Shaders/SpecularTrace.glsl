// Screenspace (clip/view space), probe space tracers implemented 
// This was a fucking nightmare to implement right 
// P.S : This file has a bunch of code that isn't really needed to function solely because I used it as sort of a "playground" for custom screenspace/probespace/voxelspace raytracers 

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

layout (location = 0) out vec4 o_SpecularIndirect;

// Intersection transversal, used as an input to the temporal denoiser 
// Can also be used to aid the spatial denoiser
layout (location = 1) out float o_Transversal;

in vec2 v_TexCoords;

uniform float u_Time;
uniform int u_Frame;
uniform vec2 u_Jitter;
uniform vec2 u_Dimensions;

uniform bool u_RoughSpecular;

uniform mat4 u_Projection;
uniform mat4 u_View;
uniform mat4 u_InverseProjection;
uniform mat4 u_InverseView;

uniform vec3 u_SunDirection;

uniform samplerCube u_ProbeAlbedo;
uniform samplerCube u_ProbeDepth;
uniform samplerCube u_ProbeNormals;

uniform samplerCube u_EnvironmentMap;

uniform vec3 u_ProbeCapturePoints[6];

uniform vec3 u_Incident;

uniform sampler2D u_Depth;
uniform sampler2D u_LowResDepth;
uniform sampler2D u_Normals;
uniform sampler2D u_LFNormals;
uniform sampler2D u_PBR;
uniform sampler2D u_Albedos;

// Voxel volumes
uniform sampler3D u_VoxelVolumes[6];
uniform usampler3D u_VoxelVolumesNormals[6];
uniform float u_VoxelRanges[6];
uniform vec3 u_VoxelCenters[6];

uniform float u_zNear;
uniform float u_zFar;

uniform sampler2D u_Shadowmap;
uniform mat4 u_SunShadowMatrix;
uniform mat4 u_ViewProjection;
uniform bool u_Checker;

uniform mat4 u_PrevView;
uniform mat4 u_PrevProjection;
uniform sampler2D u_PreviousFrameDiffuse;



// Contains all data necessary to integrate lighting for a point 
struct GBufferData {
   vec3 Position;
   float Depth;
   vec3 Normal;
   vec3 Albedo;
   vec3 Emission;
   vec3 Data;
   bool ValidMask;
   bool Approximated;
   bool SSR;
   vec3 Direction;
   float SkyAmount;
};

float HASH2SEED = 0.0f;
vec2 hash2() 
{
	return fract(sin(vec2(HASH2SEED += 0.1, HASH2SEED += 0.1)) * vec2(43758.5453123, 22578.1459123));
}

vec3 Saturation(vec3 Color, float Adjustment)
{
    const vec3 LuminosityCoefficients = vec3(0.2125f, 0.7154f, 0.0721f);
    vec3 Luminosity = vec3(dot(Color, LuminosityCoefficients));
    return mix(Luminosity, Color, Adjustment);
}

float Luminance(vec3 rgb)
{
    const vec3 W = vec3(0.2125, 0.7154, 0.0721);
    return dot(rgb, W);
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

vec3 ProjectToLastFrame(vec3 WorldPos) 
{
	vec4 ProjectedPosition = u_PrevProjection * u_PrevView * vec4(WorldPos, 1.0f);
	ProjectedPosition.xyz /= ProjectedPosition.w;
	ProjectedPosition.xyz = ProjectedPosition.xyz * 0.5f + 0.5f;
	return ProjectedPosition.xyz;
}

vec3 ProjectToClipSpace(vec3 WorldPos) 
{
	vec4 ProjectedPosition = u_ViewProjection * vec4(WorldPos, 1.0f);
	ProjectedPosition.xyz /= ProjectedPosition.w;
	return ProjectedPosition.xyz;
}

vec3 WorldPosFromDepth(float depth, vec2 txc)
{
    float z = depth * 2.0 - 1.0;
    vec4 ClipSpacePosition = vec4(txc * 2.0 - 1.0, z, 1.0);
    vec4 ViewSpacePosition = u_InverseProjection * ClipSpacePosition;
    ViewSpacePosition /= ViewSpacePosition.w;
    vec4 WorldPos = u_InverseView * ViewSpacePosition;
    return WorldPos.xyz;
}

bool SSRayValid(vec2 x) {
	float bias = 0.0001f;
	if (x.x > bias && x.x < 1.0f - bias && x.y > bias && x.y < 1.0f - bias) {
		return true;
	}

	return false;
}

vec3 GGX_VNDF(vec3 N, float roughness, vec2 Xi)
{
    float alpha = roughness * roughness;
    float alpha2 = alpha * alpha;
	
    float phi = 2.0 * PI * Xi.x;
    float cosTheta = sqrt((1.0 - Xi.y) / (1.0 + (alpha2 - 1.0) * Xi.y));
    float sinTheta = sqrt(1.0 - cosTheta * cosTheta);
	
    vec3 H;
    H.x = cos(phi) * sinTheta;
    H.y = sin(phi) * sinTheta;
    H.z = cosTheta;
	
    vec3 up = abs(N.z) < 0.999 ? vec3(0.0, 0.0, 1.0) : vec3(1.0, 0.0, 0.0);
    vec3 tangent = normalize(cross(up, N));
    vec3 bitangent = cross(N, tangent);
	
    vec3 sampleVec = tangent * H.x + bitangent * H.y + N * H.z;
    return normalize(sampleVec);
} 

vec3 SampleMicrofacet(vec3 N, float R) {

    R = max(R, 0.01f);
	float NearestDot = -100.0f;
	vec3 BestDirection = N;

	for (int i = 0 ; i < 4 ; i++) 
    {
		vec2 Xi = hash2() * vec2(0.9f, 0.85f);
        
        vec3 ImportanceSampled = GGX_VNDF(N, R, Xi);
		float d = dot(ImportanceSampled, N);

		if (d > NearestDot) {
			BestDirection = ImportanceSampled;
			NearestDot = d;
		}
	}

    if (dot(N, BestDirection) < 0.0f) {
        return N;
    }

	return BestDirection;
}

vec3 SampleMicrofacetBayer(vec3 N, float R, vec2 Xi) {

    R = max(R, 0.01f);

    Xi *= vec2(0.7f, 0.5f);

    vec3 ImportanceSampled = GGX_VNDF(N, R, Xi);

    if (dot(N, ImportanceSampled) < 0.0f) {
        return N;
    }

	return ImportanceSampled;
}


int GetFaceID(vec3 Direction)
{
    vec3 AbsoluteDirection = abs(Direction);
    float Index = 0.0f;

	if(AbsoluteDirection.z >= AbsoluteDirection.x && AbsoluteDirection.z >= AbsoluteDirection.y)
	{
		Index = Direction.z < 0.0 ? 5.0 : 4.0;
	}

	else if(AbsoluteDirection.y >= AbsoluteDirection.x)
	{
		Index = Direction.y < 0.0 ? 3.0 : 2.0;
	}

	else
	{
		Index = Direction.x < 0.0 ? 1.0 : 0.0;
	}

    return int(Index);
}

vec3 GetCapturePoint(vec3 Direction) {
    return u_ProbeCapturePoints[clamp(GetFaceID(Direction),0,5)];
}

float DistanceSqr(vec3 A, vec3 B)
{
    vec3 C = A - B;
    return dot(C, C);
}

vec3 CosWeightedHemisphere(const vec3 n) 
{
  	vec2 r = vec2(0.0f);
	r = vec2(hash2());
	float PI2 = 2.0f * PI;
	vec3  uu = normalize(cross(n, vec3(0.0,1.0,1.0)));
	vec3  vv = cross(uu, n);
	float ra = sqrt(r.y);
	float rx = ra * cos(PI2 * r.x); 
	float ry = ra * sin(PI2 * r.x);
	float rz = sqrt(1.0 - r.y);
	vec3  rr = vec3(rx * uu + ry * vv + rz * n );
    return normalize(rr);
}



const float EmissiveDesat = 0.925f;
const float EmissionStrength = 3.5f;
const bool Skylighting = true;

GBufferData Raytrace(vec3 WorldPosition, vec3 Direction, float ErrorTolerance, float Hash, int Steps, int BinarySteps) {
   
    // Settings 

    const float Distance = 384.0f;

    float StepSize = Distance / float(Steps);

    vec3 ReflectionVector = Direction; 
    
    vec3 RayPosition = WorldPosition + ReflectionVector * Hash;
    vec3 RayOrigin = RayPosition;

    vec3 PreviousSampleDirection = ReflectionVector;

    bool FoundHit = false;

    // Exponential stepping 
    float ExpStep = mix(1.01f, 1.015f, Hash);

    float SkyAmount = 0.0f;

    float PrevRayDepth = 0.0f;
    float PrevStepSampleDepth = 0.0f;

    // Find intersection with geometry 
    // Todo : Account for geometrical thickness?
    for (int CurrentStep = 0; CurrentStep < Steps ; CurrentStep++) 
    {
        StepSize *= ExpStep;

        vec3 CapturePoint = GetCapturePoint(PreviousSampleDirection);
        vec3 Diff = RayPosition - CapturePoint;
        float L = length(Diff);
        vec3 SampleDirection = Diff / L;
        PreviousSampleDirection = SampleDirection;
        float ProbeDepth = (texture(u_ProbeDepth, SampleDirection).x * 128.0f);

        if (ProbeDepth < 0.00000001f || ProbeDepth == 0.0f) {
            SkyAmount += 1.0f;
            ProbeDepth = 10000.0f;
        }

        float DepthError = abs(ProbeDepth - L);

        float ThresholdCurr = abs(L - PrevRayDepth); 
        ThresholdCurr = (clamp(ThresholdCurr * 8.0f, 0.5f, Distance) * ErrorTolerance);

        bool AccurateishHit = DepthError < ThresholdCurr;

        if (L > ProbeDepth && AccurateishHit) {
             FoundHit = true;
             break;
        }

        if (L > ProbeDepth) {
            break;
        }

        PrevRayDepth = L;
        PrevStepSampleDepth = ProbeDepth;
        RayPosition += ReflectionVector * StepSize;
        
    }

    if (FoundHit) 
    {
        // Do a basic ssr-style binary search along intersection step and find best intersection point 

        const bool DoBinaryRefinement = true;

        vec3 FinalBinaryRefinePos = RayPosition;

        if (DoBinaryRefinement) {

            float BR_StepSize = StepSize / 2.0f;
            FinalBinaryRefinePos = FinalBinaryRefinePos - ReflectionVector * BR_StepSize;

            for (int BinaryRefine = 0 ; BinaryRefine < BinarySteps; BinaryRefine++) 
            {
                BR_StepSize /= 2.0f;

                vec3 CapturePoint = GetCapturePoint(PreviousSampleDirection);
                vec3 Diff = FinalBinaryRefinePos - CapturePoint;
                float L = length(Diff);
                vec3 BR_SampleDirection = Diff / L;
                PreviousSampleDirection = BR_SampleDirection;

                float Depth = (texture(u_ProbeDepth, BR_SampleDirection).x * 128.0f);

                
                if (Depth < 0.00000001f || Depth == 0.0f) {
                    SkyAmount += 1.0f;
                    Depth = 10000.0f;
                }

                float RaySign = (Depth < L) ? -1.0f : 1.0f;
                FinalBinaryRefinePos += ReflectionVector * BR_StepSize * RaySign;
            }
        }

        vec3 CapturePoint = GetCapturePoint(PreviousSampleDirection);
        vec3 FinalVector = FinalBinaryRefinePos - CapturePoint;
        float FinalLength = length(FinalVector);
        vec3 FinalSampleDirection = FinalVector / FinalLength;


        float DepthFetch = textureLod(u_ProbeDepth, FinalSampleDirection, 0.0f).x * 128.0f;

        vec4 AlbedoFetch = textureLod(u_ProbeAlbedo, FinalSampleDirection, 0.0f) * 1.0f;

        vec4 NormalFetch = textureLod(u_ProbeNormals, FinalSampleDirection, 0.0f);

        float Fade = clamp(pow(exp(-(abs(DepthFetch-FinalLength))), 1.5f), 0.0f, 1.0f);
        float FadeStrong = clamp(pow(exp(-(abs(DepthFetch-FinalLength))), 4.0f), 0.0f, 1.0f);

        GBufferData ReturnValue;
        ReturnValue.Position = (CapturePoint + DepthFetch * FinalSampleDirection);
        ReturnValue.Normal = NormalFetch.xyz;
        ReturnValue.Albedo = AlbedoFetch.xyz * mix(0.0f, 0.5f, Fade);
        ReturnValue.Data = vec3(AlbedoFetch.w, 0.0f, 1.0f);

        // Emission

        ReturnValue.Emission = Saturation(AlbedoFetch.xyz, EmissiveDesat) * NormalFetch.w * EmissionStrength * FadeStrong;

        ReturnValue.ValidMask = true;
        ReturnValue.Approximated = false;

        ReturnValue.SSR = false;

        ReturnValue.Depth = DepthFetch;
        ReturnValue.Direction = PreviousSampleDirection;
        ReturnValue.SkyAmount = 0.0f;

        return ReturnValue;

    }


    // No hit found, return black color

    GBufferData ReturnValue;
    ReturnValue.Position = vec3(RayOrigin) + ReflectionVector * 120.0f;
    ReturnValue.Normal = vec3(0.0f);
    ReturnValue.Albedo = vec3(0.0f);
    ReturnValue.Data = vec3(0.0f, 0.0f, 1.0f);
    ReturnValue.ValidMask = false;
    ReturnValue.Approximated = true;
    ReturnValue.SSR = false;
    ReturnValue.Depth = -1.;
    ReturnValue.Emission = vec3(0.0f);
    ReturnValue.Direction = PreviousSampleDirection;
    ReturnValue.SkyAmount = SkyAmount;

    return ReturnValue;
}

float LogXy(float base, float x) {
    return log(x) / log(base);
}

GBufferData ScreenspaceRaytrace(vec3 Origin, vec3 Direction, float ThresholdMultiplier, float Hash, int Steps, int BinarySteps)
{
    const float Distance = 196.0f;

    GBufferData ReturnValue;
    ReturnValue.Position = Origin;
    ReturnValue.Normal = vec3(0.0f);
    ReturnValue.Albedo = vec3(0.0f);
    ReturnValue.Data = vec3(0.0f, 0.0f, 1.0f);
    ReturnValue.ValidMask = false;
    ReturnValue.Approximated = true;
    ReturnValue.SSR = true;
    ReturnValue.Depth = -1.;
    ReturnValue.Emission = vec3(0.0f);
    ReturnValue.SkyAmount = 0.0f;

    float StepSize = float(Distance) / float(Steps);
	vec3 RayPosition = Origin + Direction * Hash; 
	vec2 FinalUV = vec2(-1.0f);

    float ExpStep = 1.035f; //clamp((log(500.0f) / log(float(Steps))) * 0.75f, 1.05f, 6.0f);// mix(1.075f, 1.5f, float(Hash));

	for(int CurrentStep = 0; CurrentStep < Steps; CurrentStep++) 
	{
        float ToleranceStep = mix(2.2f, 6.0f, pow(float(CurrentStep) / float(Steps), 1.5f));

		float Threshold = StepSize * ThresholdMultiplier * ToleranceStep;
		
		vec3 ProjectedRayScreenspace = ProjectToClipSpace(RayPosition); 
		
		if(abs(ProjectedRayScreenspace.x) > 1.0f || abs(ProjectedRayScreenspace.y) > 1.0f || abs(ProjectedRayScreenspace.z) > 1.0f) 
		{
			return ReturnValue;
		}
		
		ProjectedRayScreenspace.xyz = ProjectedRayScreenspace.xyz * 0.5f + 0.5f; 

		if (!SSRayValid(ProjectedRayScreenspace.xy))
		{
			return ReturnValue;
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
            ReturnValue.Position = WorldPosFromDepth(FinalDepth, FinalProjected.xy);
            ReturnValue.Normal = texture(u_Normals, FinalProjected.xy).xyz;
            ReturnValue.Albedo = texture(u_Albedos, FinalProjected.xy).xyz / 2.0f;
            vec4 PBR = texture(u_PBR, FinalProjected.xy).xyzw;
            ReturnValue.Data = vec3(0.0f, 0.0f, 1.0f);
            ReturnValue.Emission = Saturation(ReturnValue.Albedo, EmissiveDesat) * PBR.w * EmissionStrength;
            ReturnValue.ValidMask = true;
            ReturnValue.Approximated = false;
            ReturnValue.SSR = true;
            ReturnValue.Depth = FinalDepth;

            return ReturnValue;
		}

        // Step 
		RayPosition += StepSize * Direction; 

        StepSize *= ExpStep;
	}

	return ReturnValue;

}


GBufferData ScreenspaceTrace_Clip(vec3 Origin, vec3 Direction, float ThresholdMultiplier, float Hash, int Steps)
{
    const int BinarySteps = 16;

    GBufferData ReturnValue;
    ReturnValue.Position = Origin;
    ReturnValue.Normal = vec3(0.0f);
    ReturnValue.Albedo = vec3(0.0f);
    ReturnValue.Data = vec3(0.0f, 0.0f, 1.0f);
    ReturnValue.ValidMask = false;
    ReturnValue.Approximated = true;
    ReturnValue.SSR = true;
    ReturnValue.Depth = -1.;
    ReturnValue.Emission = vec3(0.0f);
    ReturnValue.SkyAmount = 0.0f;


    vec3 ScreenOrigin = ProjectToScreenSpace(Origin);

    float StepSize = 1.0f / float(Steps);

    vec3 ScreenDirection = normalize(ProjectToScreenSpace(Origin + Direction) - ScreenOrigin);

    vec3 RayPosition = ScreenOrigin + ScreenDirection * Hash * StepSize;


    for (int Step = 0 ; Step < Steps ; Step++) {

        if (!SSRayValid(RayPosition.xy)) {

            return ReturnValue;

        }

        float Tolerance = StepSize * mix(pow(float(Step) / float(Steps), 2.0f), 2.0f, 7.0f); // * (ThresholdMultiplier * 1000.0f);


        float Depth = texture(u_Depth, RayPosition.xy).x;
        float LinearDepth = LinearizeDepth(Depth);
        float Error = abs(LinearizeDepth(RayPosition.z) - LinearDepth);

        if (Error < Tolerance && RayPosition.z > Depth)
        {
            // Binary search the best intersection location

            vec3 BinaryStepVector = (ScreenDirection * StepSize) / 2.0f;
            RayPosition -= BinaryStepVector;

            for(int BinaryStep = 0; BinaryStep < BinarySteps; BinaryStep++)
            {
                float DepthError = texture(u_Depth, RayPosition.xy).r - RayPosition.z;
                float Sign = sign(DepthError);
                RayPosition += Sign * BinaryStepVector;
                BinaryStepVector /= 2.0f;
            }

            // Generate gbuffer data and return 

            float FinalDepth = Depth;
            vec2 FinalProjected = RayPosition.xy;

            ReturnValue.Position = WorldPosFromDepth(FinalDepth, FinalProjected.xy);
            ReturnValue.Normal = texture(u_Normals, FinalProjected.xy).xyz;
            ReturnValue.Albedo = texture(u_Albedos, FinalProjected.xy).xyz / 2.0f;
            vec4 PBR = texture(u_PBR, FinalProjected.xy).xyzw;
            ReturnValue.Data = vec3(0.0f, 0.0f, 1.0f);
            ReturnValue.Emission = Saturation(ReturnValue.Albedo, EmissiveDesat) * PBR.w * EmissionStrength;
            ReturnValue.ValidMask = true;
            ReturnValue.Approximated = false;
            ReturnValue.SSR = true;
            ReturnValue.Depth = FinalDepth;

            return ReturnValue;
        }


        RayPosition += ScreenDirection * StepSize;
    }

	return ReturnValue;

}

// Reproject intersection to previous frame and try to estimate indirect lighting 
vec3 ReprojectIndirect(vec3 P, vec3 A, float R, vec3 Bp) 
{

    vec3 Projected = ProjectToLastFrame(P);

    // Todo : Add distance check
    if (SSRayValid(Projected.xy)) {
        vec4 Fetch = texture(u_PreviousFrameDiffuse, Projected.xy).xyzw;
        return Fetch.xyz * 1.0f * A * pow(Fetch.w,2.5f);//mix(2.5f, 3.5f, clamp(R * 2.0f,0.,1.)));
    }

    vec3 ProjectedBase = ProjectToLastFrame(Bp);

    if (SSRayValid(ProjectedBase.xy)) {
        return texture(u_PreviousFrameDiffuse, ProjectedBase.xy).xyz * 0.9f * A;
    }

    return texture(u_PreviousFrameDiffuse, v_TexCoords).xyz * 0.5f * A;

    return texture(u_EnvironmentMap, vec3(0.0f, 1.0f, 0.0f)).xyz * 0.2f * A * mix(1.0f, 2.5f, float(R<0.2f));
}


const vec3 SUN_COLOR = vec3(8.0f);

// Integrates lighting for a point
vec3 IntegrateLighting(GBufferData Hit, vec3 Direction, const bool FilterShadow, float R, float D, vec3 Origin) {
    
    if (!Hit.ValidMask) {
       return pow(texture(u_EnvironmentMap, Hit.Direction).xyz, vec3(1.2f)) * clamp(Hit.SkyAmount, 0.0f, 1.0f) * 2.7f * float(Skylighting);
    }

    // Sky 
    if (Hit.Depth > 1498.0f && !Hit.SSR) {
        
        return pow(Hit.Albedo, vec3(1.0f / 1.3f)) * 1.4f;

    }

    Hit.Normal = normalize(Hit.Normal);

    float Shadow = 0.0f; 

    // Approximated hits mean that the hits have a very high error
    bool DoShadowMap = !Hit.Approximated;

    if (DoShadowMap) {

        const vec2 Poisson[6] = vec2[6](vec2(-0.613392, 0.617481),  vec2(0.751946, 0.453352),
                                        vec2(0.170019, -0.040254),  vec2(0.078707, -0.715323),
                                        vec2(-0.299417, 0.791925),  vec2(-0.075838, -0.529344));
        
        vec4 ProjectionCoordinates = u_SunShadowMatrix * vec4(Hit.Position + Hit.Normal * 0.5f, 1.0f);
	    ProjectionCoordinates.xyz = ProjectionCoordinates.xyz / ProjectionCoordinates.w;
        ProjectionCoordinates.xyz = ProjectionCoordinates.xyz * 0.5f + 0.5f;
        
        float Depth = ProjectionCoordinates.z;
	    float Bias = 0.0045f;  

        if (FilterShadow) {
           
            
            vec2 TexelSize = 1.0f / textureSize(u_Shadowmap, 0).xy;

            // pcf 
            for (int i = 0 ; i < 4 ; i++) {
                float Fetch = texture(u_Shadowmap, ProjectionCoordinates.xy + Poisson[i] * TexelSize * 2.4f).x;
                Shadow += float(ProjectionCoordinates.z - Bias > Fetch);
            }

            Shadow /= 4.0f;

        }

        else {

            float SimpleFetch = texture(u_Shadowmap, ProjectionCoordinates.xy).x;
            Shadow = float(ProjectionCoordinates.z - Bias > SimpleFetch);
        }

        Shadow = 1.0f - Shadow;
    }

    // Lambert BRDF  
    // Todo : Switch to hammon diffuse brdf (ignore specular brdf to reduce variance)
    float Lambertian = max(0.0f, dot(Hit.Normal, -u_SunDirection));
    vec3 Direct = Lambertian * SUN_COLOR * 0.75f * Shadow * Hit.Albedo;
    vec3 Indirect = ReprojectIndirect(Hit.Position, Hit.Albedo, R, Origin);
    return Direct + Indirect + Hit.Emission;
}



vec4 DecodeVolumeLighting(const vec4 Lighting) {

	return vec4(Lighting.xyz * 5.0f, Lighting.w);
}

vec4 ResolveVoxelRadiance(vec4 Lighting, vec3 P, vec3 Bp, float R) {
    
    Lighting = DecodeVolumeLighting(Lighting);

    const float Al = Luminance(vec3(1.0f) * vec3(0.01f));
    float L = Luminance(Lighting.xyz);

    // No lighting here, use approximate indirect 
    if (L > 0.000001f && L < Al) {

        vec3 Albedo = Lighting.xyz / 0.01f;
        return vec4(ReprojectIndirect(P, Albedo, R, Bp), Lighting.w) * 0.6f;
    }

    return Lighting;
}

bool InScreenspace(vec3 x) {

	if (x.x > 0.0f && x.x < 1.0f && x.y > 0.0f && x.y < 1.0f && x.z > 0.0f && x.z < 1.0f) { 
		return true;
	}

	return false;
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

bool LiesInsideVolume(vec3 vp) { return abs(vp.x) < 128.0f && abs(vp.y) < 128.0f  && abs(vp.z) < 128.0f ; } 

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

int GetCascadeNumber(vec3 P, int MinCascade) {

	for (int Cascade = int(MinCascade) ; Cascade < 6 ; Cascade++) {
		
		if (PositionInVolume(Cascade, P)) {
			return Cascade;
		}

	}

	return 5;
}

vec4 RaymarchCascades(vec3 WorldPosition, vec3 Normal, vec3 Direction, float Aperature, float LowDiscrepHash, const int Steps, int MinCascade, out vec3 RayPositionO) 
{
	const float Diagonal = sqrt(2.0f);

	float HashScale = mix(1.0f, 2.0f, LowDiscrepHash / 1.41f);

	vec3 RayOrigin = WorldPosition;

	int CurrentCascade = GetCascadeNumber(RayOrigin, MinCascade);

	float Distance = (u_VoxelRanges[CurrentCascade] / 128.0f) * Diagonal;
	RayOrigin += Direction * Distance;

	vec3 RayPosition = RayOrigin + Direction * Distance * 2.0f * LowDiscrepHash;
	float StepSize = u_VoxelRanges[CurrentCascade] / 128.0f;

	vec4 TotalGI = vec4(0.0f);
	TotalGI.w = 1.0f;

	float SkyVisibility = 1.0f;

	bool IntersectionFound = true;

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
	SkyVisibility = pow(SkyVisibility, 24.0f);
	vec3 SkySample = pow(texture(u_EnvironmentMap, Direction).xyz, vec3(2.0f));
	float LambertSky = 1.; //pow(clamp(dot(VoxelNormal, vec3(0.0f, 1.0f, 0.0f)), 0.0f, 1.0f), 1.0f);
	vec3 SkyRadiance = SkySample * SkyVisibility * 0.8f;
    RayPositionO = RayPosition;
	return vec4(SkyRadiance + TotalGI.xyz, IntersectionFound ? distance(WorldPosition, RayPosition) : 64.0f);
}


// Temporal upscale offsets  
ivec2 UpscaleOffsets2x2[] = ivec2[](
	ivec2(1, 1),
	ivec2(1, 0),
	ivec2(0, 0),
	ivec2(0, 1));

const ivec2[16] UpscaleOffsets4x4 = ivec2[16](
    ivec2(0, 0),
    ivec2(2, 0),
    ivec2(0, 2),
    ivec2(2, 2),
    ivec2(1, 1),
    ivec2(3, 1),
    ivec2(1, 3),
    ivec2(3, 3),
    ivec2(1, 0),
    ivec2(3, 0),
    ivec2(1, 2),
    ivec2(3, 2),
    ivec2(0, 1),
    ivec2(2, 1),
    ivec2(0, 3),
    ivec2(2, 3)
);

// Returns whether the pixel is part of the sky or not 
bool IsSky(float NonLinearDepth) {
    if (NonLinearDepth > 0.9999992f || NonLinearDepth == 1.0f) {
        return true;
	}

    return false;
}


// Trace settings 

// Do voxel raytracing when probe fails?
const bool VX_TRACE = false;



void main() {
    
    vec2 TexCoordJittered = v_TexCoords;

    HASH2SEED = (TexCoordJittered.x * TexCoordJittered.y) * 64.0;

    // Animate noise for temporal integration
	HASH2SEED += fract(u_Time) * 64.0f;

    // Calculate pixel 
    ivec2 Pixel = ivec2(gl_FragCoord.xy);

    if (u_Checker) {
        Pixel.x *= 2;
	    bool IsCheckerStep = Pixel.x % 2 == int(Pixel.y % 2 == (u_Frame % 2));
        Pixel.x += int(IsCheckerStep);
    }

    // Jitter for temporal super sampling 
    //Pixel += ivec2(UpscaleOffsets4x4[u_Frame % 16]);
    Pixel += ivec2(u_Jitter * 3.0f);

    // Constant resolution (0.5x)
    ivec2 HighResPixel = Pixel * 2;

    vec2 HighResUV = vec2(HighResPixel) / textureSize(u_Depth, 0).xy;

    // GBuffer fetches 
    float Depth = texelFetch(u_Depth, HighResPixel, 0).x;

	// Sky check
    if (IsSky(Depth)) {
        o_SpecularIndirect.xyz = vec3(0.0f);
        o_Transversal = 256.0f;
        o_SpecularIndirect.w = o_Transversal;
        return;
    }

    // Sample gbuffers 
	vec3 WorldPosition = WorldPosFromDepth(Depth, HighResUV);
    vec3 Normal = normalize(texelFetch(u_Normals, HighResPixel, 0).xyz); 
    vec3 LFNormal = normalize(texelFetch(u_LFNormals, HighResPixel, 0).xyz); 
    vec3 PBR = texelFetch(u_PBR, HighResPixel, 0).xyz;
    vec3 Incident = normalize(WorldPosition - u_Incident);

    float Roughness = u_RoughSpecular ? PBR.x : 0.1f;

    // Intersection tolerance (Rougher surfaces can have less accurate reflections with a less noticable quality loss)
    float Tolerance = 1.0f;//mix(1.0, 1.8f, pow(Roughness, 1.5f));
    float ToleranceSS = 0.001f;//mix(0.00145f, 0.04, Roughness * Roughness);

    // Sample integration
    vec3 TotalRadiance = vec3(0.0f);
    float AverageTransversal = 0.0f;
    float TotalWeight = 0.0f;

    // Bias roughness (to reduce noise)
    const float RoughnessBias = 0.98f;
    // Epic remapping -> //float BiasedRoughness = clamp(pow((Roughness * RoughnessBias) + 1.0f, 2.0f) / 8.0f, 0.0f, 1.0f); 
    float BiasedRoughness = max(Roughness * 0.925, 0.07f);//pow(Roughness, 1.25f) * RoughnessBias; 
   
    bool FilterShadowMap = Roughness <= 0.5 + 0.01f;

    bool DoScreenspaceTrace = true;

    // Screenspace tracing is only done if roughness < 0.425
    int SSSteps = Roughness < 0.2f ? 64 : (Roughness < 0.3f ? 40 : 24);
    int SSBinarySteps = Roughness < 0.2f ? 16 : 12;

    // Calculate steps 
    int ProbeSteps = int(mix(128.0f, 48.0f, BiasedRoughness < 0.25f ? pow(BiasedRoughness,1.5f) : BiasedRoughness));
    int ProbeBinarySteps = (BiasedRoughness <= 0.51f) ? 16 : 12;


    for (int Sample ; Sample < SAMPLES ; Sample++) {

        
        float BayerHash = fract(fract(mod(float(u_Frame) + float(Sample) * 2., 384.0f) * (1.0 / PHI)) + Bayer32(gl_FragCoord.xy));

        // Sample microfacet normal from VNDF
        vec3 Microfacet;

        if (u_RoughSpecular) {
            Microfacet = SampleMicrofacet(Normal, BiasedRoughness);
        }
            
        else {
            Microfacet = LFNormal;
        }

        const float Bias_n = 2.75f;


        // Raytrace 

        vec3 Direction = normalize(reflect(Incident, Microfacet));


        GBufferData Intersection;

        if (!VX_TRACE) {


            if (DoScreenspaceTrace && BiasedRoughness <= 0.125f)
            {
                // Trace in screen space 
                Intersection = ScreenspaceRaytrace(WorldPosition + LFNormal * Bias_n, Direction, ToleranceSS, BayerHash, SSSteps, SSBinarySteps);
                
                // If that fails, trace in probe space  
                if (!Intersection.ValidMask) {

                    Intersection = Raytrace(WorldPosition + LFNormal * Bias_n, Direction, Tolerance, BayerHash, ProbeSteps, ProbeBinarySteps);
                }
            }
            
            else {
                Intersection = Raytrace(WorldPosition + LFNormal * Bias_n, Direction, Tolerance, BayerHash, ProbeSteps, ProbeBinarySteps);
            }

            // Integrate lighting for hit point
            vec3 CurrentRadiance = IntegrateLighting(Intersection, Direction, FilterShadowMap, BiasedRoughness, LinearizeDepth(Depth), WorldPosition);

            // Sum up radiance 
            TotalRadiance += CurrentRadiance;

            float CurrentTransversal = Intersection.ValidMask ? distance(Intersection.Position, WorldPosition) : 256.0f;
            AverageTransversal += CurrentTransversal;
        }

        else {

            Intersection = Raytrace(WorldPosition + LFNormal * Bias_n, Direction, Tolerance, BayerHash, ProbeSteps, ProbeBinarySteps);
            
            if (!Intersection.ValidMask) {

                vec3 VXIntersection;
                vec4 VXRadiance = RaymarchCascades(WorldPosition + LFNormal * 1.41414f, Normal, Direction, 1.0f, BayerHash, 192, 0, VXIntersection);
                VXRadiance = ResolveVoxelRadiance(VXRadiance,VXIntersection,WorldPosition,Roughness);
                float CurrentTransversal = VXRadiance.w;
                AverageTransversal += CurrentTransversal;
                TotalRadiance += VXRadiance.xyz;
            }

            else {

                vec3 CurrentRadiance = IntegrateLighting(Intersection, Direction, FilterShadowMap, BiasedRoughness, LinearizeDepth(Depth), WorldPosition);
                TotalRadiance += CurrentRadiance;

                float CurrentTransversal = Intersection.ValidMask ? distance(Intersection.Position, WorldPosition) : 256.0f;
                AverageTransversal += CurrentTransversal;
            }


        }

        // Add weight 
        TotalWeight += 1.0f;
    }

    TotalRadiance /= max(TotalWeight, 1.0f);
    AverageTransversal /= max(TotalWeight, 1.0f);

    o_SpecularIndirect.xyz = TotalRadiance;
    o_Transversal = AverageTransversal / 64.0f;

    // Nan/inf check
    if (isnan(o_SpecularIndirect.x) || isnan(o_SpecularIndirect.y) || isnan(o_SpecularIndirect.z) || isinf(o_SpecularIndirect.x) || isinf(o_SpecularIndirect.y) || isinf(o_SpecularIndirect.z)) {
        o_SpecularIndirect.xyz = vec3(0.0f);
    }

    if (isnan(o_Transversal) || isinf(o_Transversal)) {
        o_Transversal = 0.0f;
        o_SpecularIndirect.w = 0.0f;
    }

    o_SpecularIndirect.w = o_Transversal;
}