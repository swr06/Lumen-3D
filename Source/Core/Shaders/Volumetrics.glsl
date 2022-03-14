#version 330 core

#define Bayer4(a)   (Bayer2(  0.5 * (a)) * 0.25 + Bayer2(a))
#define Bayer8(a)   (Bayer4(  0.5 * (a)) * 0.25 + Bayer2(a))
#define Bayer16(a)  (Bayer8(  0.5 * (a)) * 0.25 + Bayer2(a))
#define Bayer32(a)  (Bayer16( 0.5 * (a)) * 0.25 + Bayer2(a))
#define Bayer64(a)  (Bayer32( 0.5 * (a)) * 0.25 + Bayer2(a))
#define Bayer128(a) (Bayer64( 0.5 * (a)) * 0.25 + Bayer2(a))
#define Bayer256(a) (Bayer128(0.5 * (a)) * 0.25 + Bayer2(a))
#define PHI 1.6180339
#define PI 3.14159265359 


layout (location = 0) out vec4 o_Volumetrics;

in vec2 v_TexCoords;

uniform sampler2D u_LowResDepth;
uniform sampler2D u_LFNormals;

uniform sampler2D u_HistoryVolumetrics;
uniform sampler2D u_HistoryDepth;

uniform float u_zNear;
uniform float u_zFar;

uniform vec3 u_ViewerPosition;

uniform sampler2D u_Shadowmap;
uniform mat4 u_SunShadowMatrix;
uniform vec3 u_SunDirection;

uniform mat4 u_Projection;
uniform mat4 u_View;
uniform mat4 u_InverseProjection;
uniform mat4 u_InverseView;
uniform mat4 u_InversePrevView;

uniform mat4 u_PrevProjection;
uniform mat4 u_PrevView;

uniform float u_SunVLStrength;
uniform bool u_VolumetricsTemporal;

uniform float u_Time;
uniform int u_Frame;

uniform vec2 u_Jitter;

// Noise functions ->
float noise(in vec3 x); // by iq

// Bayer hash function ->
float Bayer2(vec2 a) 
{
    a = floor(a);
    return fract(dot(a, vec2(0.5, a.y * 0.75)));
}

// GBuffer functions ->
vec3 WorldPosFromDepth(float depth, vec2 txc)
{
    float z = depth * 2.0 - 1.0;
    vec4 ClipSpacePosition = vec4(txc * 2.0 - 1.0, z, 1.0);
    vec4 ViewSpacePosition = u_InverseProjection * ClipSpacePosition;
    ViewSpacePosition /= ViewSpacePosition.w;
    vec4 WorldPos = u_InverseView * ViewSpacePosition;
    return WorldPos.xyz;
}

vec3 Incident(vec2 screenspace)
{
	vec4 clip = vec4(screenspace * 2.0f - 1.0f, -1.0, 1.0);
	vec4 eye = vec4(vec2(u_InverseProjection * clip), -1.0, 0.0);
	return normalize(vec3(u_InverseView * eye));
}

bool IsSky(float NonLinearDepth) {
    if (NonLinearDepth > 0.9999992f || NonLinearDepth == 1.0f) {
        return true;
	}

    return false;
}

float linearizeDepth(float depth)
{
	return (2.0 * u_zNear) / (u_zFar + u_zNear - depth * (u_zFar - u_zNear));
}

vec3 Reprojection(vec3 WorldPos) 
{
	vec4 ProjectedPosition = u_PrevProjection * u_PrevView * vec4(WorldPos, 1.0f);
	ProjectedPosition.xyz /= ProjectedPosition.w;
	ProjectedPosition.xyz = ProjectedPosition.xyz * 0.5f + 0.5f;
	return ProjectedPosition.xyz;
}


// Phase ->

mat3x3 Extinction = mat3x3(vec3(5.8, 13.3, 33.31) * 1e-6, vec3(21.0) * 1e-6 * 1.11, vec3(8.30428e-07, 1.31491e-06, 5.44068e-08));

float CornetteShanks(float cosTheta, float g)
{
    const float cornette = 3.0f / (8.0f * PI);
    float gg = g*g;
    float num   = (1.0f - gg) * (1.0f + pow(cosTheta, 2.0f));
    float denom = (2.0f + gg) * pow(1.0f + gg - 2.0f * g * cosTheta, 1.5f);
    return cornette * (num / denom);
}

float Phase(float cosa, float g) // Normalized HG
{
	float g_sqr = g * g;
	float num = (1 - abs(g));
	float denom = sqrt(max(1 - 2 * g*cosa + g_sqr, 0));
	float frac = num / denom;
	float scale = g_sqr + (1 - g_sqr) / (4 * PI);
	return scale * (frac*frac*frac);
}

float Rayleigh(float CosTheta) 
{
    const float R = 3.0f / (16.0f * PI);
    return R * (1.0f + pow(CosTheta, 2.0f));
}

float IsotropicPhase() {
    return 0.250f / PI;
}

// Samples direct lighting for a point 
float SimpleDirect(vec3 P) 
{
    vec4 ProjectionCoordinates = u_SunShadowMatrix * vec4(P, 1.0f);
	ProjectionCoordinates.xyz = ProjectionCoordinates.xyz / ProjectionCoordinates.w;
    ProjectionCoordinates.xyz = ProjectionCoordinates.xyz * 0.5f + 0.5f;

    if (ProjectionCoordinates.xy == clamp(ProjectionCoordinates.xy, 0., 1.)) {
        float Depth = ProjectionCoordinates.z;
	    float Bias = 0.002f;  
        float SimpleFetch = texture(u_Shadowmap, ProjectionCoordinates.xy).x;
        return 1.0f - float(ProjectionCoordinates.z - Bias > SimpleFetch);
    }

    else {
        return 0.;
    }
}

// Density function 
float SampleDensity(vec3 Point) {
    
    return 1.; 

    float Noise = noise(Point * 0.05f);

    return clamp(Noise * 5.0f, 0., 1.);
    
    // return 1.0f;
}

// Samples transmittance at a ray
float GetTransmittance(vec3 Point, vec3 Direction, int Steps) {

    float Sum = 0.0f;
    float StepSize = 10.0f;

    Point += Direction * 0.5f * StepSize;

    for (int x = 0 ; x < Steps ; x++) {

        Sum += SampleDensity(Point);
        
        Point += Direction * StepSize;
    }

    return Sum;
}

vec4 Reproject(vec3 WorldPosition, vec4 Data, float Transversal, vec3 Direction) {

    bool DoReprojection = true;

    if (!DoReprojection || !u_VolumetricsTemporal) {
        return Data;
    }

    vec3 Reprojected = Reprojection(WorldPosition);

    float bias = 0.015f;

    if (Reprojected.x > bias && Reprojected.x < 1.0f-bias &&
		Reprojected.y > bias && Reprojected.y < 1.0f-bias) {
        bool CameraMoved = distance(u_InverseView[3].xyz, u_InversePrevView[3].xyz) > 0.0001f;
        vec4 History = texture(u_HistoryVolumetrics, Reprojected.xy);

        float PreviousDepth = texture(u_HistoryDepth, Reprojected.xy).x;
		float LinearPrevDepth = linearizeDepth(PreviousDepth);
		float LinearExpectedDepth = linearizeDepth(Reprojected.z);

        float BlendFactor = CameraMoved ? 0.0f : 0.75f;

        BlendFactor *= pow(exp(-abs(LinearExpectedDepth-LinearPrevDepth)), CameraMoved ? 96.0f * 1. : 24.0f);

        return mix(Data, History, BlendFactor); //CameraMoved ? 0.1f : 0.8);
    }

    return Data;
}

vec2 g_TexCoords;

void main() {

    g_TexCoords = v_TexCoords + (u_Jitter * (1.0f / textureSize(u_HistoryVolumetrics, 0).xy));
    
    float HashAnimated = fract(fract(mod(float(u_Frame) + float(0.) * 2., 384.0f) * (1.0 / PHI)) + Bayer16(gl_FragCoord.xy));
    
    bool Checkerstep = int(gl_FragCoord.x + gl_FragCoord.y) % 2 == (u_Frame % 2);
    int Steps = 32;//int(mix(32.,40.,float(Checkerstep)));
    const float MaxDistance = 384.0f;

    float Hash = Bayer64(gl_FragCoord.xy);

    float Depth = texture(u_LowResDepth, g_TexCoords).x;
    vec3 Normals = texture(u_LFNormals, g_TexCoords).xyz;

    vec3 WorldPosition = WorldPosFromDepth(Depth, g_TexCoords);  
    vec3 RawWorldPosition = WorldPosition;

    float Distance = IsSky(Depth) ? MaxDistance : distance(WorldPosition, u_ViewerPosition);
    Distance = clamp(Distance, 0.000001f, MaxDistance);

    vec3 Direction = Incident(g_TexCoords);
    vec3 EndPosition = u_ViewerPosition + Direction * Distance;

    float StepSize = Distance / float(Steps);

    vec3 RayPosition = u_ViewerPosition + Direction * HashAnimated;

    const float G = 0.7f;

    float CosTheta = dot(-Direction, normalize(u_SunDirection));
    float DirectPhase = IsotropicPhase();//clamp(CornetteShanks(CosTheta, G), 0.0f, IsotropicPhase());

    vec3 Transmittance = vec3(1.0f);

    vec3 DirectScattering = vec3(0.0f);

    float Extinction = 0.0f; // extinction for air is (very close to) zero, this will change if the volume differs 

    vec3 SunColor = (vec3(253.,184.,100.)/255.0f) * 0.12f * u_SunVLStrength * 0.3333f;

    for (int Step = 0 ; Step < Steps ; Step++) {

        float Density = SampleDensity(RayPosition);

        if (Density <= 0.0f) {
            continue;
        }

        float DirectVisibility = SimpleDirect(RayPosition);
        vec3 Direct = DirectVisibility * DirectPhase * SunColor;
        //vec3 Indirect = mix(0.0f, 1.0f, 1.0f - clamp(RayPosition.y / 300.0f, 0.0f, 1.0f)) * vec3(0.5f, 0.5f, 1.0f) * 0.0004f;
        vec3 S = (Direct) * StepSize * Density * Transmittance;

        DirectScattering += S;
        Transmittance *= exp(-(StepSize * Density) * Extinction);
        RayPosition += Direction * StepSize * mix(1.0f, Hash, 0.8f);
    }

    vec4 Data = vec4(vec3(DirectScattering), Transmittance);

    o_Volumetrics = Reproject(IsSky(Depth) ? u_ViewerPosition + Direction * 128.0f : RawWorldPosition, Data, 1.0f, Direction);
}



// 3D Noise ->

float hash1(float n) // by iq
{
    return fract( n*17.0*fract( n*0.3183099 ) );
}

float noise(in vec3 x) // by iq
{
    vec3 p = floor(x);
    vec3 w = fract(x);
    vec3 u = w*w*w*(w*(w*6.0-15.0)+10.0);
    float n = p.x + 317.0*p.y + 157.0*p.z;
    float a = hash1(n+0.0);
    float b = hash1(n+1.0);
    float c = hash1(n+317.0);
    float d = hash1(n+318.0);
    float e = hash1(n+157.0);
	float f = hash1(n+158.0);
    float g = hash1(n+474.0);
    float h = hash1(n+475.0);
    float k0 =   a;
    float k1 =   b - a;
    float k2 =   c - a;
    float k3 =   e - a;
    float k4 =   a - b - c + d;
    float k5 =   a - c - e + g;
    float k6 =   a - b - e + f;
    float k7 = - a + b + c - d + e - f - g + h;
    return -1.0+2.0*(k0 + k1*u.x + k2*u.y + k3*u.z + k4*u.x*u.y + k5*u.y*u.z + k6*u.z*u.x + k7*u.x*u.y*u.z);
}