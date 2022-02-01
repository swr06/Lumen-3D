#version 400 core

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

in vec2 v_TexCoords;

uniform float u_Time;

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

vec3 WorldPosFromDepth(float depth, vec2 txc)
{
    float z = depth * 2.0 - 1.0;
    vec4 ClipSpacePosition = vec4(txc * 2.0 - 1.0, z, 1.0);
    vec4 ViewSpacePosition = u_InverseProjection * ClipSpacePosition;
    ViewSpacePosition /= ViewSpacePosition.w;
    vec4 WorldPos = u_InverseView * ViewSpacePosition;
    return WorldPos.xyz;
}

struct GBufferData {
   vec3 Position;
   vec3 Normal;
   vec3 Albedo;
   vec3 Data;
};

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

GBufferData Raytrace(vec3 WorldPosition, vec3 Normal, float Depth) {
    const float Distance = 256.0f;

    // 120 + 8 = 128 total steps 
    const int Steps = 120;
    const int BinarySteps = 8;

    float StepSize = Distance / float(Steps);
    float UnditheredStepSize = StepSize;

    float Hash = Bayer32(gl_FragCoord.xy);
    StepSize *= mix(Hash, 1.0f, 0.8f);

    vec3 ViewDirection = normalize(WorldPosition - u_Incident);
    vec3 ReflectionVector = normalize(reflect(ViewDirection, Normal));
    vec3 RayPosition = (WorldPosition + Normal * 1.0f) + ReflectionVector * StepSize * 2.0f;

    float BestError = 100000.0f;
    vec3 TraceColor = vec3(0.0f);
    vec3 BestPosition = RayPosition;

    vec3 PreviousSampleDirection = ReflectionVector;
    
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
             break;
        }

        RayPosition += ReflectionVector * StepSize;
        
    }

    vec3 IntersectionPosition = vec3(0.0f);

    // Binary refinement ->

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

    return ReturnValue;
}

    const vec3 SUN_COLOR = vec3(6.9f, 6.9f, 10.0f);


vec3 BRDF(GBufferData Hit) {
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
    vec3 Direct = Lambertian * SUN_COLOR * 0.5f * Shadow * Hit.Albedo;
    vec3 Ambient = texture(u_EnvironmentMap, vec3(0.0f, 1.0f, 0.0f)).xyz * 0.2f * Hit.Albedo;
    return Direct + Ambient;
}

void main() {

    float Depth = texture(u_Depth, v_TexCoords).x;
	vec3 WorldPosition = WorldPosFromDepth(Depth, v_TexCoords);
    vec3 Normal = texture(u_Normals, v_TexCoords).xyz; 
    GBufferData Intersection = Raytrace(WorldPosition, Normal, Depth);
    o_Color = BRDF(Intersection);
}