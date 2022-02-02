#version 400 core

#define PI 3.14159265359
#define PHI 1.6180339
#define SAMPLES 1

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
layout (location = 1) out float o_Transversal;

in vec2 v_TexCoords;

uniform float u_Time;
uniform int u_Frame;
uniform vec2 u_Jitter;
uniform vec2 u_Dimensions;

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
uniform sampler2D u_PBR;

uniform sampler2D u_Shadowmap;
uniform mat4 u_SunShadowMatrix;

struct GBufferData {
   vec3 Position;
   vec3 Normal;
   vec3 Albedo;
   vec3 Data;
   bool ValidMask;
};

float HASH2SEED = 0.0f;
vec2 hash2() 
{
	return fract(sin(vec2(HASH2SEED += 0.1, HASH2SEED += 0.1)) * vec2(43758.5453123, 22578.1459123));
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

vec3 ImportanceSampleGGX(vec3 N, float roughness, vec2 Xi)
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

    R *= 0.75f;
    R = max(R, 0.05f);
	float NearestDot = -100.0f;
	vec3 BestDirection;

	for (int i = 0 ; i < 4 ; i++) 
    {
		vec2 Xi = hash2() * vec2(0.9f, 0.8f);
        
        vec3 ImportanceSampled = ImportanceSampleGGX(N, R, Xi);
		float d = dot(ImportanceSampled, N);

		if (d > NearestDot) {
			BestDirection = ImportanceSampled;
			NearestDot = d;
		}
	}

	return BestDirection;
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

GBufferData Raytrace(vec3 WorldPosition, vec3 Normal, float Depth, float Hash, out vec3 MicrofacetReflected) {
    const float Distance = 384.0f;

    const int Steps = 192;
    const int BinarySteps = 8;

    float StepSize = Distance / float(Steps);
    float UnditheredStepSize = StepSize;

    StepSize *= mix(Hash, 1.0f, 0.8f);

    vec3 ViewDirection = normalize(WorldPosition - u_Incident);
    vec3 ReflectionVector = normalize(reflect(ViewDirection, Normal)); 
    MicrofacetReflected = ReflectionVector;
    vec3 RayPosition = (WorldPosition + Normal * 1.0f) + ReflectionVector * StepSize * 2.0f;
    vec3 RayOrigin = RayPosition;

    float BestError = 100000.0f;
    vec3 TraceColor = vec3(0.0f);
    vec3 BestPosition = RayPosition;

    vec3 PreviousSampleDirection = ReflectionVector;

    bool FoundHit = false;
    
    // Find intersection with geometry 
    // You probably want to take geometrical thickness into account here as well 
    for (int CurrentStep = 0; CurrentStep < Steps ; CurrentStep++) 
    {
        vec3 CapturePoint = GetCapturePoint(PreviousSampleDirection);
        vec3 SampleDirection = normalize(RayPosition - CapturePoint);
        PreviousSampleDirection = SampleDirection;
        vec3 SamplePosition = SampleDirection * (texture(u_ProbeDepth, SampleDirection).x * 128.0f);
        vec3 SampleWorldPosition = SamplePosition + CapturePoint;
        float Error = DistanceSqr(SampleWorldPosition, RayPosition); 

        if (Error < 2.82f) {
             BestError = Error;
             BestPosition = RayPosition;
             FoundHit = true;
             break;
        }

        RayPosition += ReflectionVector * StepSize;
        
    }

    vec3 IntersectionPosition = vec3(0.0f);

    // Binary refinement ->

    if (FoundHit) 
    {
        float BR_StepSize = UnditheredStepSize / 2.0f;
        vec3 FinalBinaryRefinePos = BestPosition;

        for (int BinaryRefine = 0 ; BinaryRefine < BinarySteps; BinaryRefine++) 
        {
            vec3 CapturePoint = GetCapturePoint(PreviousSampleDirection);
            vec3 BR_SampleDirection = normalize(FinalBinaryRefinePos - CapturePoint);
            PreviousSampleDirection = BR_SampleDirection;
            vec3 BR_SamplePosition = BR_SampleDirection * (texture(u_ProbeDepth, BR_SampleDirection).x * 128.0f);
            vec3 BR_SampleWorldPosition = BR_SamplePosition + CapturePoint;
            IntersectionPosition = BR_SampleWorldPosition;

            if (DistanceSqr(FinalBinaryRefinePos, BR_SampleWorldPosition) < BestError) 
            {
                FinalBinaryRefinePos -= ReflectionVector * BR_StepSize;
            }

            else 
            {
                FinalBinaryRefinePos += ReflectionVector * BR_StepSize;
            }

            BR_StepSize /= 2.0f;
        }

        vec3 CapturePoint = GetCapturePoint(PreviousSampleDirection);
        vec3 FinalSampleDirection = normalize(FinalBinaryRefinePos - CapturePoint);

        // Return 
        vec4 AlbedoFetch = texture(u_ProbeAlbedo, FinalSampleDirection);
        vec4 NormalFetch = texture(u_ProbeNormals, FinalSampleDirection);

        GBufferData ReturnValue;
        ReturnValue.Position = IntersectionPosition;
        ReturnValue.Normal = NormalFetch.xyz;
        ReturnValue.Albedo = AlbedoFetch.xyz;
        ReturnValue.Data = vec3(AlbedoFetch.w, NormalFetch.w, 1.0f);
        ReturnValue.ValidMask = true;
        return ReturnValue;

    }

     GBufferData ReturnValue;
     ReturnValue.Position = vec3(RayOrigin) + ReflectionVector * 700.0f;
     ReturnValue.Normal = vec3(0.0f);
     ReturnValue.Albedo = vec3(0.0f);
     ReturnValue.Data = vec3(0.0f, 0.0f, 1.0f);
     ReturnValue.ValidMask = false;

     return ReturnValue;
}



const vec3 SUN_COLOR = vec3(6.9f, 6.9f, 10.0f);

vec3 BRDF(GBufferData Hit, vec3 Direction) {
    
    if (!Hit.ValidMask) {
        return texture(u_EnvironmentMap, Direction).xyz * 0.6f;
    }

    const vec2 Poisson[6] = vec2[6](vec2(-0.613392, 0.617481),  vec2(0.751946, 0.453352),
                                    vec2(0.170019, -0.040254),  vec2(0.078707, -0.715323),
                                    vec2(-0.299417, 0.791925),  vec2(-0.075838, -0.529344));
    
    vec4 ProjectionCoordinates = u_SunShadowMatrix * vec4(Hit.Position, 1.0f);
	ProjectionCoordinates.xyz = ProjectionCoordinates.xyz / ProjectionCoordinates.w;
    ProjectionCoordinates.xyz = ProjectionCoordinates.xyz * 0.5f + 0.5f;
    float Depth = ProjectionCoordinates.z;
	float Bias = 0.0008f;  
    float Fetch = texture(u_Shadowmap, ProjectionCoordinates.xy).x;
    float Shadow = 0.0f; 
    
    vec2 TexelSize = 1.0f / textureSize(u_Shadowmap, 0).xy;

    // pcf 
    for (int i = 0 ; i < 6 ; i++) {
        float Fetch = texture(u_Shadowmap, ProjectionCoordinates.xy + Poisson[i] * TexelSize * 2.0f).x;
        Shadow += float(ProjectionCoordinates.z - Bias > Fetch);
    }

    Shadow /= 6.0f;
    Shadow = 1.0f - Shadow;

    float Lambertian = max(0.0f, dot(Hit.Normal, -u_SunDirection));
    vec3 Direct = Lambertian * SUN_COLOR * 0.25f * Shadow * Hit.Albedo;
    vec3 Ambient = texture(u_EnvironmentMap, vec3(0.0f, 1.0f, 0.0f)).xyz * 0.2f * Hit.Albedo;
    return Direct + Ambient;
}

void main() {
    vec2 JitteredTexCoords = v_TexCoords;
    JitteredTexCoords += u_Jitter / u_Dimensions;

    HASH2SEED = (JitteredTexCoords.x * JitteredTexCoords.y) * 64.0;
	HASH2SEED += fract(u_Time) * 64.0f;

    float BayerHash = fract(fract(mod(float(u_Frame), 256.0f) * (1.0 / PHI)) + Bayer32(gl_FragCoord.st));
    float Depth = texture(u_Depth, JitteredTexCoords).x;
	vec3 WorldPosition = WorldPosFromDepth(Depth, JitteredTexCoords);
    vec3 Normal = texture(u_Normals, JitteredTexCoords).xyz; 
    vec3 PBR = texture(u_PBR, JitteredTexCoords).xyz;

    vec3 Microfacet = SampleMicrofacet(Normal, PBR.x);
    vec3 Reflected = vec3(0.0f);
    GBufferData Intersection = Raytrace(WorldPosition, Microfacet, Depth, BayerHash, Reflected);
    o_Color = BRDF(Intersection, Reflected);

    o_Transversal = distance(Intersection.Position, WorldPosition);

    if (isnan(o_Color.x) || isnan(o_Color.y) || isnan(o_Color.z) || isinf(o_Color.x) || isinf(o_Color.y) || isinf(o_Color.z)) {
        o_Color = vec3(0.0f);
    }
}