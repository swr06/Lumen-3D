#version 400 core

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

layout (location = 0) out vec3 o_Color;

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
uniform sampler2D u_Normals;
uniform sampler2D u_LFNormals;
uniform sampler2D u_PBR;
uniform sampler2D u_Albedos;

uniform float u_zNear;
uniform float u_zFar;

uniform sampler2D u_Shadowmap;
uniform mat4 u_SunShadowMatrix;
uniform mat4 u_ViewProjection;
uniform bool u_Checker;

// Contains all data necessary to integrate lighting for a point 
struct GBufferData {
   vec3 Position;
   vec3 Normal;
   vec3 Albedo;
   vec3 Data;
   bool ValidMask;
   bool Approximated;
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


const float EmissiveDesat = 0.15f;
const float EmissionStrength = 96.0f;

GBufferData Raytrace(vec3 Incident, vec3 WorldPosition, vec3 Normal, vec3 LFNormal, float ErrorTolerance, float Hash, out vec3 MicrofacetReflected) {
   
    // If enabled, the raytracer returns the nearest hit which is approximated using a weighting factor 
    const bool FALLBACK_ON_BEST_STEP = true;


    // Bias 
    WorldPosition = (WorldPosition + LFNormal * 2.0f);

    // Settings 
    const float Distance = 512.0f;
    const int Steps = 160;
    const int BinarySteps = 20;


    float StepSize = Distance / float(Steps);
    float UnditheredStepSize = StepSize;
    StepSize *= mix(Hash * 1.0f, 1.0f, 0.75f);

    vec3 ReflectionVector = normalize(reflect(Incident, Normal)); //CosWeightedHemisphere(LFNormal); 
    vec3 ReflectionVectorLF = normalize(reflect(Incident, LFNormal));
    
    MicrofacetReflected = ReflectionVector;

    vec3 RayPosition = WorldPosition + ReflectionVectorLF * StepSize * 0.5f;
    vec3 RayOrigin = RayPosition;

    vec3 TraceColor = vec3(0.0f);

    vec3 PreviousSampleDirection = ReflectionVector;

    bool FoundHit = false;

    // Exponential stepping 
    int ExponentialStepStart = Steps - (Steps / 4);
    float ExpStep = mix(1.0f, 1.7f, mix(Hash, 1.0f, 0.2f));

    // Approximate hits when we can't find an accurate intersection with the geometry 
    vec3 BestSampleDirection = RayPosition;
    float BestRayWeight = -10.0f;
    
    // Find intersection with geometry 
    // Todo : Account for geometrical thickness?
    for (int CurrentStep = 0; CurrentStep < Steps ; CurrentStep++) 
    {
        if (CurrentStep > ExponentialStepStart) {
            StepSize *= ExpStep;
        }

        vec3 CapturePoint = GetCapturePoint(PreviousSampleDirection);
        vec3 Diff = RayPosition - CapturePoint;
        float L = length(Diff);
        vec3 SampleDirection = Diff / L;
        PreviousSampleDirection = SampleDirection;
        vec3 SamplePosition = SampleDirection * (texture(u_ProbeDepth, SampleDirection).x * 128.0f);
        vec3 SampleWorldPosition = SamplePosition + CapturePoint;
        float Error = DistanceSqr(SampleWorldPosition, RayPosition); 

        float IntersectionThresh = StepSize * 1.5f;//mix(StepSize * 0.6f, StepSize * 1.5f, clamp(float(CurrentStep) * (1.0f / 12.0f), 0.0f, 1.0f));

        if (Error < IntersectionThresh * ErrorTolerance) {
             FoundHit = true;
             break;
        }

        // Compute ray weighting factor 
        float RayWeight = pow(1.0f / float(CurrentStep + 1.0f), 6.0f) * pow(1.0f / max(Error, 0.00000001f), 8.725f);

        // Weight rays
        if (RayWeight > BestRayWeight) {
            BestSampleDirection = SampleDirection;
            BestRayWeight = RayWeight;
        }

        RayPosition += ReflectionVector * StepSize;
        
    }


    // Binary refine intersection point 

    if (FoundHit) 
    {
        const bool DoBinaryRefinement = true;

        vec3 FinalBinaryRefinePos = RayPosition;

        if (DoBinaryRefinement) {

            float BR_StepSize = UnditheredStepSize / 2.0f;
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

                if (Depth < L) 
                {
                    FinalBinaryRefinePos -= ReflectionVector * BR_StepSize;
                }

                else 
                {
                    FinalBinaryRefinePos += ReflectionVector * BR_StepSize;
                }
            }
        }

        vec3 CapturePoint = GetCapturePoint(PreviousSampleDirection);
        vec3 FinalSampleDirection = normalize(FinalBinaryRefinePos - CapturePoint);

        // Return 
        vec4 AlbedoFetch = textureLod(u_ProbeAlbedo, FinalSampleDirection, 0.0f);
        vec4 NormalFetch = textureLod(u_ProbeNormals, FinalSampleDirection, 0.0f);
        float DepthFetch = textureLod(u_ProbeDepth, FinalSampleDirection, 0.0f).x * 128.0f;

        GBufferData ReturnValue;
        ReturnValue.Position = (CapturePoint + DepthFetch * FinalSampleDirection);
        ReturnValue.Normal = NormalFetch.xyz;
        ReturnValue.Albedo = AlbedoFetch.xyz;
        ReturnValue.Data = vec3(AlbedoFetch.w, 0.0f, 1.0f);

        // Emission
        ReturnValue.Albedo += Saturation(ReturnValue.Albedo, EmissiveDesat) * NormalFetch.w * EmissionStrength;

        ReturnValue.ValidMask = true;
        ReturnValue.Approximated = false;

        return ReturnValue;

    }

    // Best error 
    if (FALLBACK_ON_BEST_STEP) {

        vec3 CapturePoint = GetCapturePoint(BestSampleDirection);
        vec3 FinalSampleDirection = BestSampleDirection;
        vec4 AlbedoFetch = textureLod(u_ProbeAlbedo, FinalSampleDirection, 0.0f);
        vec4 NormalFetch = textureLod(u_ProbeNormals, FinalSampleDirection, 0.0f);
        float DepthFetch = textureLod(u_ProbeDepth, FinalSampleDirection, 0.0f).x * 128.0f;

        GBufferData ReturnValue;
        ReturnValue.Position = (CapturePoint + DepthFetch * FinalSampleDirection);
        ReturnValue.Normal = NormalFetch.xyz;
        ReturnValue.Albedo = AlbedoFetch.xyz;
        ReturnValue.Data = vec3(AlbedoFetch.w, 0.0f, 1.0f);

        // Emission
        ReturnValue.Albedo += Saturation(ReturnValue.Albedo, EmissiveDesat) * NormalFetch.w * EmissionStrength;
        ReturnValue.ValidMask = true;
        ReturnValue.Approximated = true;

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

    return ReturnValue;
}

GBufferData ScreenspaceRaytrace(vec3 Incident, vec3 Origin, vec3 Normal, vec3 LFNormal, float ThresholdMultiplier, float Hash, out vec3 MicrofacetReflected)
{
    const float RayDistance = 38.0f;
    const int Steps = 16;

    GBufferData ReturnValue;
    ReturnValue.Position = Origin;
    ReturnValue.Normal = vec3(0.0f);
    ReturnValue.Albedo = vec3(0.0f);
    ReturnValue.Data = vec3(0.0f, 0.0f, 1.0f);
    ReturnValue.ValidMask = false;
    ReturnValue.Approximated = true;

    vec3 Direction = normalize(reflect(Incident, Normal));
    MicrofacetReflected = Direction;
	vec3 StepVector = (Direction * RayDistance) / float(Steps); 
	vec3 RayPosition = Origin + StepVector * Hash; 
	vec2 FinalUV = vec2(-1.0f);

	for(int CurrentStep = 0; CurrentStep < Steps; CurrentStep++) 
	{
		float Threshold = length(StepVector) * ThresholdMultiplier;
		
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
		
		float DepthAt = texture(u_Depth, ProjectedRayScreenspace.xy).x; 
		float CurrentRayDepth = LinearizeDepth(ProjectedRayScreenspace.z); 
		float Error = abs(LinearizeDepth(DepthAt) - CurrentRayDepth);
		
        // Intersected!
		if (Error < Threshold && ProjectedRayScreenspace.z > DepthAt) 
		{
			// Basic Binary refinement : 

            bool DoBinaryRefinement = true;

            vec3 FinalProjected = vec3(0.0f);
            float FinalDepth = 0.0f;

            if (DoBinaryRefinement) {
			    vec3 BinaryStepVector = StepVector / 2.0f;

                RayPosition -= BinaryStepVector;

			    for (int BinaryStep = 0 ; BinaryStep < 7 ; BinaryStep++) {
			    		
			    	BinaryStepVector /= 2.0f;

			    	vec3 Projected = ProjectToClipSpace(RayPosition); 
			    	Projected = Projected * 0.5f + 0.5f;
                    FinalProjected = Projected;
                    float Fetch = texture(u_Depth, Projected.xy).x;
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
            ReturnValue.Albedo = texture(u_Albedos, FinalProjected.xy).xyz * 1.0f;
            vec4 PBR = texture(u_PBR, FinalProjected.xy).xyzw;
            ReturnValue.Data = vec3(PBR.x, 0.0f, 1.0f);
            ReturnValue.Albedo += Saturation(ReturnValue.Albedo, EmissiveDesat) * PBR.w * EmissionStrength;
            ReturnValue.ValidMask = true;
            ReturnValue.Approximated = false;

            return ReturnValue;
		}

        // Step 
		RayPosition += StepVector; 

	}

	return ReturnValue;

}

const vec3 SUN_COLOR = vec3(6.9f, 6.9f, 10.0f);

vec3 IntegrateLighting(GBufferData Hit, vec3 Direction) {
    
    if (!Hit.ValidMask) {
        return vec3(0.0f);
       // return texture(u_EnvironmentMap, Direction).xyz * 0.3f;
    }

    float Shadow = 0.0f; 

    // Approximated hits mean that the hits have a very high error
    bool DoShadowMap = !Hit.Approximated;

    if (DoShadowMap) {

        const vec2 Poisson[6] = vec2[6](vec2(-0.613392, 0.617481),  vec2(0.751946, 0.453352),
                                        vec2(0.170019, -0.040254),  vec2(0.078707, -0.715323),
                                        vec2(-0.299417, 0.791925),  vec2(-0.075838, -0.529344));
        
        vec4 ProjectionCoordinates = u_SunShadowMatrix * vec4(Hit.Position + Hit.Normal * 0.0f, 1.0f);
	    ProjectionCoordinates.xyz = ProjectionCoordinates.xyz / ProjectionCoordinates.w;
        ProjectionCoordinates.xyz = ProjectionCoordinates.xyz * 0.5f + 0.5f;
        float Depth = ProjectionCoordinates.z;
	    float Bias = 0.004f;  
        float Fetch = texture(u_Shadowmap, ProjectionCoordinates.xy).x;
        
        vec2 TexelSize = 1.0f / textureSize(u_Shadowmap, 0).xy;

        // pcf 
        for (int i = 0 ; i < 6 ; i++) {
            float Fetch = texture(u_Shadowmap, ProjectionCoordinates.xy + Poisson[i] * TexelSize * 2.0f).x;
            Shadow += float(ProjectionCoordinates.z - Bias > Fetch);
        }

        Shadow /= 6.0f;
        Shadow = 1.0f - Shadow;
    }

    // Lambert BRDF  
    // Todo : Switch to hammon diffuse brdf (ignore specular brdf to reduce variance)
    float Lambertian = max(0.0f, dot(Hit.Normal, -u_SunDirection));
    vec3 Direct = Lambertian * SUN_COLOR * 0.07f * Shadow * Hit.Albedo;
    vec3 Ambient = texture(u_EnvironmentMap, vec3(0.0f, 1.0f, 0.0f)).xyz * 0.2f * Hit.Albedo;
    return Direct + Ambient;
}



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

bool IsSky(float NonLinearDepth) {
    if (NonLinearDepth > 0.99998f) {
        return true;
	}

    return false;
}

void main() {
    
    // For temporal super sampling 
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

    Pixel += UpscaleOffsets4x4[u_Frame % 16];

    // Constant resolution (0.5x)
    ivec2 HighResPixel = Pixel * 2;

    vec2 HighResUV = vec2(HighResPixel) / textureSize(u_Depth, 0).xy;

    // GBuffer fetches 
    float Depth = texelFetch(u_Depth, HighResPixel, 0).x;

	// Sky check
    if (IsSky(Depth)) {
        o_Color = vec3(0.0f);
        o_Transversal = 32.0f;
        return;
    }

	vec3 WorldPosition = WorldPosFromDepth(Depth, HighResUV);
    vec3 Normal = texelFetch(u_Normals, HighResPixel, 0).xyz; 
    vec3 LFNormal = texelFetch(u_LFNormals, HighResPixel, 0).xyz; 
    vec3 PBR = texelFetch(u_PBR, HighResPixel, 0).xyz;

    float Roughness = PBR.x;

    float Tolerance = mix(1.0, 2.5f, !u_RoughSpecular ? 0.0f : pow(PBR.x, 1.5f));
    float ToleranceSS = mix(0.00175, 0.005f, !u_RoughSpecular ? 0.0f : pow(PBR.x, 1.5f));

    // Sample integration
    vec3 TotalRadiance = vec3(0.0f);
    float AverageTransversal = 0.0f;
    float TotalWeight = 0.0f;

    const float RoughnessBias = 0.85f;

    vec3 Incident = normalize(WorldPosition - u_Incident);

    bool DoScreenspaceTrace = true;


    for (int Sample ; Sample < SAMPLES ; Sample++) {

        float BayerHash = fract(fract(mod(float(u_Frame) + float(Sample) * 2., 384.0f) * (1.0 / PHI)) + Bayer32(gl_FragCoord.xy));

        // Sample microfacet normal from VNDF
        vec3 Microfacet;

        if (u_RoughSpecular) {
            Microfacet = SampleMicrofacet(Normal, PBR.x * RoughnessBias);
        }
            
        else {
            Microfacet = Normal;
        }

        vec3 ReflectedDirection = vec3(0.0f);
        
        // Raytrace!

        GBufferData Intersection;

        if (DoScreenspaceTrace) {
            // Trace in screen space 
            Intersection = ScreenspaceRaytrace(Incident, WorldPosition + LFNormal * 0.8f, Microfacet, LFNormal, ToleranceSS, BayerHash, ReflectedDirection);
            
            // If that fails, trace in probe space  
            if (!Intersection.ValidMask) {
                Intersection = Raytrace(Incident, WorldPosition, Microfacet, LFNormal, Tolerance, BayerHash, ReflectedDirection);
            }
        }

        else {
            Intersection = Raytrace(Incident, WorldPosition, Microfacet, LFNormal, Tolerance, BayerHash, ReflectedDirection);
        }

        // Integrate lighting 
        vec3 CurrentRadiance = IntegrateLighting(Intersection, ReflectedDirection);

        // Sum up radiance 
        TotalRadiance += CurrentRadiance;

        // Store transversal for filtering/reprojection
        float CurrentTransversal = Intersection.ValidMask ? distance(Intersection.Position, WorldPosition) : 64.0f;
        CurrentTransversal /= 64.0f;
        AverageTransversal += CurrentTransversal;

        // Add weight 
        TotalWeight += 1.0f;
    }

    TotalRadiance /= max(TotalWeight, 1.0f);
    AverageTransversal /= max(TotalWeight, 1.0f);

    o_Color = TotalRadiance;
    o_Transversal = AverageTransversal;

    // Nan/inf check
    if (isnan(o_Color.x) || isnan(o_Color.y) || isnan(o_Color.z) || isinf(o_Color.x) || isinf(o_Color.y) || isinf(o_Color.z)) {
        o_Color = vec3(0.0f);
    }

    if (isnan(o_Transversal) || isinf(o_Transversal)) {
        o_Transversal = 0.0f;
    }
}