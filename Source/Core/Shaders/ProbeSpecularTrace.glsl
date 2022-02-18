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

uniform float u_zNear;
uniform float u_zFar;

uniform sampler2D u_Shadowmap;
uniform mat4 u_SunShadowMatrix;
uniform mat4 u_ViewProjection;
uniform bool u_Checker;

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

	for (int i = 0 ; i < 3 ; i++) 
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
const float EmissionStrength = 9.5f;

GBufferData Raytrace(vec3 WorldPosition, vec3 Direction, float ErrorTolerance, float Hash) {
   
    // If enabled, the raytracer returns the nearest hit which is approximated using a weighting factor 
    const bool FALLBACK_ON_BEST_STEP = false;

    // Settings 
    const float Distance = 384.0f;
    const int Steps = 180;
    const int BinarySteps = 16;

    float StepSize = Distance / float(Steps);
    float UnditheredStepSize = StepSize;

    vec3 ReflectionVector = Direction; //CosWeightedHemisphere(LFNormal); 
    
    vec3 RayPosition = WorldPosition + ReflectionVector * Hash;
    vec3 RayOrigin = RayPosition;

    vec3 TraceColor = vec3(0.0f);

    vec3 PreviousSampleDirection = ReflectionVector;

    bool FoundHit = false;

    // Exponential stepping 
    int ExponentialStepStart = Steps - (Steps / 4);
    float ExpStep = 1.04f;//mix(1.0f, 1.4f, mix(Hash, 1.0f, 0.2f));

    // Approximate hits when we can't find an accurate intersection with the geometry 
    vec3 BestSampleDirection = RayPosition;
    float BestRayWeight = -10.0f;

    float PrevRayDepth = 0.0f;
    float PrevStepSampleDepth = 0.0f;

    int NumExpSteps = 0;

    // Find intersection with geometry 
    // Todo : Account for geometrical thickness?
    for (int CurrentStep = 0; CurrentStep < Steps ; CurrentStep++) 
    {
        if (CurrentStep > Steps / 4) {
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
        bool DistanceHitWeight = DepthError < 8.0f;
        //bool DistanceBasedHit = DepthError < ErrorTolerance * length(StepSize) * 2.5f;

        if (L > ProbeDepth && AccurateishHit) {
             FoundHit = true;
             break;
        }

        if (FALLBACK_ON_BEST_STEP) {
            // Compute ray weighting factor 
            float RayWeight = pow(1.0f / float(CurrentStep + 1.0f), 6.0f) * pow(1.0f / max(abs(ProbeDepth - L), 0.00000001f), 8.725f);

            // Weight rays
            if (RayWeight > BestRayWeight) {
                BestSampleDirection = SampleDirection;
                BestRayWeight = RayWeight;
            }
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
                float RaySign = (Depth < L) ? -1.0f : 1.0f;
                FinalBinaryRefinePos += ReflectionVector * BR_StepSize * RaySign;
            }
        }

        vec3 CapturePoint = GetCapturePoint(PreviousSampleDirection);
        vec3 FinalVector = FinalBinaryRefinePos - CapturePoint;
        float FinalLength = length(FinalVector);
        vec3 FinalSampleDirection = FinalVector / FinalLength;


        float DepthFetch = textureLod(u_ProbeDepth, FinalSampleDirection, 0.0f).x * 128.0f;
        vec4 AlbedoFetch = textureLod(u_ProbeAlbedo, FinalSampleDirection, 0.0f);
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
        ReturnValue.Emission = Saturation(ReturnValue.Albedo, EmissiveDesat) * NormalFetch.w * 8.0f * EmissionStrength;
        ReturnValue.ValidMask = true;
        ReturnValue.Approximated = true;
        ReturnValue.SSR = false;

        ReturnValue.Depth = -1.;

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

    return ReturnValue;
}

GBufferData ScreenspaceRaytrace(vec3 Origin, vec3 Direction, float ThresholdMultiplier, float Hash, int Steps)
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

    float StepSize = float(Distance) / float(Steps);
	vec3 RayPosition = Origin + Direction * Hash; 
	vec2 FinalUV = vec2(-1.0f);

    float ExpStep = 1.05f;// mix(1.075f, 1.5f, float(Hash));

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
			    
                for (int BinaryStep = 0 ; BinaryStep < 16 ; BinaryStep++) {
			    		
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



const vec3 SUN_COLOR = vec3(6.9f, 6.9f, 10.0f);

vec3 IntegrateLighting(GBufferData Hit, vec3 Direction, const bool FilterShadow) {
    
    if (!Hit.ValidMask) {
        return vec3(0.0f);
       // return texture(u_EnvironmentMap, Direction).xyz * 0.3f;
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
    vec3 Direct = Lambertian * SUN_COLOR * 0.07f * Shadow * Hit.Albedo * 3.0f;
    vec3 FakeIndirect = texture(u_EnvironmentMap, vec3(0.0f, 1.0f, 0.0f)).xyz * 0.2f * Hit.Albedo;
    return Direct + FakeIndirect + Hit.Emission;
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
    Pixel += ivec2(u_Jitter * 2.0f);

    // Constant resolution (0.5x)
    ivec2 HighResPixel = Pixel * 2;

    vec2 HighResUV = vec2(HighResPixel) / textureSize(u_Depth, 0).xy;

    // GBuffer fetches 
    float Depth = texelFetch(u_Depth, HighResPixel, 0).x;

	// Sky check
    if (IsSky(Depth)) {
        o_SpecularIndirect.xyz = vec3(0.0f);
        o_Transversal = 64.0f;
        o_SpecularIndirect.w = o_Transversal;
        return;
    }

    // Sample gbuffers 
	vec3 WorldPosition = WorldPosFromDepth(Depth, HighResUV);
    vec3 Normal = normalize(texelFetch(u_Normals, HighResPixel, 0).xyz); 
    vec3 LFNormal = normalize(texelFetch(u_LFNormals, HighResPixel, 0).xyz); 
    vec3 PBR = texelFetch(u_PBR, HighResPixel, 0).xyz;
    vec3 Incident = normalize(WorldPosition - u_Incident);

    float Roughness = u_RoughSpecular ? PBR.x : 0.0f;

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
    float BiasedRoughness = pow(Roughness, 1.25f) * RoughnessBias; 
   
    bool FilterShadowMap = Roughness <= 0.5 + 0.01f;

    bool DoScreenspaceTrace = true;

    int SSSteps = BiasedRoughness <= 0.69f + 0.01f ? (BiasedRoughness <= 0.1f ? 56 : 42) : (BiasedRoughness > 0.825f ? 24 : 32); // nice


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

        vec3 Direction = normalize(reflect(Incident, Microfacet));

        // Raytrace!

        GBufferData Intersection;

        const float Bias_n = 2.75f;

        if (DoScreenspaceTrace) {
            // Trace in screen space 
            Intersection = ScreenspaceRaytrace(WorldPosition + LFNormal * Bias_n, Direction, ToleranceSS, BayerHash, SSSteps);
            
            // If that fails, trace in probe space  
            if (!Intersection.ValidMask) {

                Intersection = Raytrace(WorldPosition + LFNormal * Bias_n, Direction, Tolerance, BayerHash);
            }
        }

        else {
            Intersection = Raytrace(WorldPosition + LFNormal * Bias_n, Direction, Tolerance, BayerHash);
        }

        // Integrate lighting 
        vec3 CurrentRadiance = IntegrateLighting(Intersection, Direction, FilterShadowMap);

        // Sum up radiance 
        TotalRadiance += CurrentRadiance;

        // Store transversal for filtering/reprojection
        float CurrentTransversal = Intersection.ValidMask ? distance(Intersection.Position, WorldPosition) : 64.0f;
        AverageTransversal += CurrentTransversal;

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