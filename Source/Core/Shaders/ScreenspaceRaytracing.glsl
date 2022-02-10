#version 330 core 

#define PI 3.14159265359
#define PHI 1.6180339

#define Bayer4(a)   (Bayer2(  0.5 * (a)) * 0.25 + Bayer2(a))
#define Bayer8(a)   (Bayer4(  0.5 * (a)) * 0.25 + Bayer2(a))
#define Bayer16(a)  (Bayer8(  0.5 * (a)) * 0.25 + Bayer2(a))
#define Bayer32(a)  (Bayer16( 0.5 * (a)) * 0.25 + Bayer2(a))
#define Bayer64(a)  (Bayer32( 0.5 * (a)) * 0.25 + Bayer2(a))
#define Bayer128(a) (Bayer64( 0.5 * (a)) * 0.25 + Bayer2(a))
#define Bayer256(a) (Bayer128(0.5 * (a)) * 0.25 + Bayer2(a))

const float TAU = radians(360.0f);
const float PHI2 = sqrt(5.0f) * 0.5f + 0.5f;
const float GOLDEN_ANGLE = TAU / PHI2 / PHI2;

float Bayer2(vec2 a) 
{
    a = floor(a);
    return fract(dot(a, vec2(0.5, a.y * 0.75)));
}

layout (location = 0) out float o_AO;
layout (location = 1) out float o_Direct;

in vec2 v_TexCoords;

uniform float u_Time;
uniform int u_Frame;

uniform mat4 u_Projection;
uniform mat4 u_View;
uniform mat4 u_ViewProjection;
uniform mat4 u_InverseProjection;
uniform mat4 u_InverseView;

uniform sampler2D u_Depth;
uniform sampler2D u_Normals;
uniform sampler2D u_PBR;

uniform vec3 u_SunDirection;

uniform float u_zNear;
uniform float u_zFar;

uniform bool u_Shadow;
uniform bool u_AO;

uniform vec2 u_Jitter;
uniform vec2 u_Dimensions;

uniform bool u_Checkerboard;


float HASH2SEED = 0.0f;
vec2 hash2() 
{
	return fract(sin(vec2(HASH2SEED += 0.1, HASH2SEED += 0.1)) * vec2(43758.5453123, 22578.1459123));
}

float LinearizeDepth(float depth)
{
	return (2.0 * u_zNear) / (u_zFar + u_zNear - depth * (u_zFar - u_zNear));
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

vec3 ProjectToViewSpace(vec3 WorldPos) 
{
	vec4 ProjectedPosition = u_View * vec4(WorldPos, 1.0f);
	return ProjectedPosition.xyz;
}


bool RayValid(vec2 x) {
	if (x == clamp(x, 0.0003f, 0.9997)) {
		return true;
	}

	return false;
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

// x, y = transversal, occlusion amount 
vec2 Raytrace(vec3 Origin, vec3 Direction, float RayDistance, int Steps, float Threshold, float Hash)
{
	vec3 StepVector = (Direction * RayDistance) / float(Steps); 
	vec3 RayPosition = Origin + StepVector * Hash; 
	vec2 FinalUV = vec2(-1.0f);

	for(int CurrentStep = 0; CurrentStep < Steps; CurrentStep++) 
	{
		vec3 ProjectedRayScreenspace = ProjectToClipSpace(RayPosition); 
		
		if(abs(ProjectedRayScreenspace.x) > 1.0f || abs(ProjectedRayScreenspace.y) > 1.0f || abs(ProjectedRayScreenspace.z) > 1.0f) 
		{
			return vec2(-1.0f, 1.0f); 
		}
		
		ProjectedRayScreenspace.xyz = ProjectedRayScreenspace.xyz * 0.5f + 0.5f; 

		if (ProjectedRayScreenspace.xy != clamp(ProjectedRayScreenspace.xy, 0.0f, 1.0f))
		{
			return vec2(-1.0f, 1.0f);
		}
		
		float DepthAt = texture(u_Depth, ProjectedRayScreenspace.xy).x; 
		float CurrentRayDepth = LinearizeDepth(ProjectedRayScreenspace.z); 
		float Error = abs(LinearizeDepth(DepthAt) - CurrentRayDepth);
		
		if (Error < Threshold && ProjectedRayScreenspace.z > DepthAt) 
		{
			//return float(StepVector) / float(Steps);

			// Binary refinement : 

			bool DoBinaryRefinement = true;

			if (DoBinaryRefinement) {

				vec3 BinaryStepVector = StepVector / 2.0f;

				for (int BinaryStep = 0 ; BinaryStep < 3 ; BinaryStep++) {
					
					vec3 Projected = ProjectToClipSpace(RayPosition); 
					Projected = Projected * 0.5f + 0.5f;
					float BinaryDepthAt = LinearizeDepth(texture(u_Depth, Projected.xy).x); 
					float BinaryRayDepth = LinearizeDepth(Projected.z); 

					if (BinaryDepthAt < BinaryRayDepth) {
						RayPosition -= BinaryStepVector;

					}

					else {
						RayPosition += BinaryStepVector;

					}

					BinaryStepVector /= 2.0f;
				}
			}

			// Calculate world space transversal  
			vec2 FinalPosition = ProjectToClipSpace(RayPosition).xy * 0.5f + 0.5f;
			float Transversal = distance(WorldPosFromDepth(texture(u_Depth, FinalPosition).x, FinalPosition), Origin.xyz);

			// Invalid hit
			if (Transversal > RayDistance) {
				return vec2(-1.0f, 1.0f);
			}

			float Occlusion = Transversal / RayDistance;
			return vec2(Transversal, Occlusion);
		}


		RayPosition += StepVector; 

	}

	return vec2(-1.0f, 1.0f);

}

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

ivec2 UpscaleOffsets2x2[4] = ivec2[](
	ivec2(1, 1),
	ivec2(1, 0),
	ivec2(0, 0),
	ivec2(0, 1));


bool IsSky(float NonLinearDepth) {
    if (NonLinearDepth > 0.99998f) {
        return true;
	}

    return false;
}

void main() {

	vec2 JitteredTexCoords = v_TexCoords;
	JitteredTexCoords += u_Jitter * (1.0f / u_Dimensions);

    HASH2SEED = (JitteredTexCoords.x * JitteredTexCoords.y) * 64.0;
	HASH2SEED += fract(u_Time) * 64.0f;

	ivec2 Pixel = ivec2(gl_FragCoord.xy);

    if (u_Checkerboard) {
        Pixel.x *= 2;
	    bool IsCheckerStep = Pixel.x % 2 == int(Pixel.y % 2 == (u_Frame % 2));
        Pixel.x += int(IsCheckerStep);
    }

	// Temporal upscale
    Pixel += UpscaleOffsets2x2[u_Frame % 4];

    ivec2 HighResPixel = Pixel * 2;
    vec2 HighResUV = vec2(HighResPixel) / textureSize(u_Depth, 0).xy;

    // GBuffer fetches 
    float Depth = texelFetch(u_Depth, HighResPixel, 0).x;

	// Sky check
	if (IsSky(Depth)) {
		o_AO = 1.0f;
		o_Direct = 1.0f;
		return;
	}

	vec3 WorldPosition = WorldPosFromDepth(Depth, HighResUV);
    vec3 Normal = texelFetch(u_Normals, HighResPixel, 0).xyz; 
	float LinearizedDepth = LinearizeDepth(Depth);

	vec3 ViewDirection = normalize(WorldPosition - u_InverseView[3].xyz);

    float BayerHash = fract(fract(mod(float(u_Frame), 384.0f) * (1.0 / PHI)) + Bayer32(gl_FragCoord.xy));

	const float AOScale = 0.8f;

	vec3 AODirection = mix(CosWeightedHemisphere(Normal), Normal, 1.0f - AOScale);
	AODirection = normalize(AODirection);

	vec3 Hash3D = vec3(hash2(), hash2().x);
	vec3 ContactShadowDirection = -u_SunDirection;
	ContactShadowDirection = normalize(ContactShadowDirection);

	float BiasDirect = mix(0.375f, 0.505f, float(LinearizedDepth > 32.0f));
	float BiasAO = mix(0.38f, 0.505f, float(LinearizedDepth > 24.0f));
	float IndirectAO = u_AO ? Raytrace(WorldPosition + Normal * BiasAO, AODirection, 20.0f, 52, 0.02f, BayerHash).y : 1.0f;
	float DirectContactShadow = u_Shadow ? Raytrace(WorldPosition + Normal * BiasDirect, ContactShadowDirection, 18.0f, 128, 0.0025f, BayerHash).y : 1.0f;
	
	o_AO = IndirectAO;
	o_Direct = DirectContactShadow;
}