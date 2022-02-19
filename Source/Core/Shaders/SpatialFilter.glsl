#version 330 core

#define PI 3.14159265359

#define saturate(x) (clamp(x, 0.0f, 1.0f))

layout (location = 0) out vec4 o_Specular;

in vec2 v_TexCoords;

uniform sampler2D u_Specular;
uniform sampler2D u_Depth;
uniform sampler2D u_Normals;
uniform sampler2D u_PBR;

uniform mat4 u_Projection;
uniform mat4 u_View;
uniform mat4 u_InverseView;
uniform mat4 u_InverseProjection;

uniform mat4 u_PrevView;
uniform mat4 u_PrevProjection;
uniform mat4 u_PrevInverseView;
uniform mat4 u_PrevInverseProjection;

uniform vec3 u_ViewerPosition;

uniform float u_zNear;
uniform float u_zFar;

uniform int u_StepSize;
uniform bool u_Checkered;

// Spherical Gaussian 
struct SG {
	vec3 Axis;
	float Sharpness;
	float Amplitude;
};

vec3 SampleIncidentRayDirection(vec2 screenspace)
{
	vec4 clip = vec4(screenspace * 2.0f - 1.0f, -1.0, 1.0);
	vec4 eye = vec4(vec2(u_InverseProjection * clip), -1.0, 0.0);
	return vec3(u_InverseView * eye);
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

// Calculates a spherical gaussian whose sharpness and amplitude are similar to that of a specular lobe for a roughness value r
// Inspired by UE4
SG RoughnessLobe(float Roughness, vec3 Normal, vec3 Incident) {
	
	Roughness = max(Roughness, 0.001f);
	float a = Roughness * Roughness;
	float a2 = a * a;

	float NDotV = clamp(abs(dot(Incident, Normal)) + 0.00001f, 0.0f, 1.0f);
	vec3 SGAxis = 2.0f * NDotV * Normal - Incident;

	SG ReturnValue;
	ReturnValue.Axis = SGAxis;
	ReturnValue.Sharpness = 0.5f / (a2 * max(NDotV, 0.1f));
	ReturnValue.Amplitude = 1.0f / (PI * a2);

	return ReturnValue;
}

float GetLobeWeight(float CenterRoughness, float SampleRoughness, vec3 CenterNormal, vec3 SampleNormal, const vec3 Incident) {
	
	// Lobe similarity weight strength
	const float Beta = 32.0f;

	float LobeSimilarity = 1.0f;
	float AxisSimilarity = 1.0f;

	SG CenterLobe = RoughnessLobe(CenterRoughness, CenterNormal, Incident);
	SG SampleLobe = RoughnessLobe(SampleRoughness, SampleNormal, Incident);

	float OneOverSharpnessSum = 1.0f / (CenterLobe.Sharpness + SampleLobe.Sharpness);

	LobeSimilarity = pow(2.0f * sqrt(CenterLobe.Sharpness * SampleLobe.Sharpness) * OneOverSharpnessSum, Beta);
	AxisSimilarity = exp(-(Beta * (CenterLobe.Sharpness * SampleLobe.Sharpness) * OneOverSharpnessSum) * clamp(1.0f - dot(CenterLobe.Axis, SampleLobe.Axis), 0.0f, 1.0f));

	return LobeSimilarity * AxisSimilarity;
}

float SpecularWeight(float CenterDepth, float SampleDepth, float CenterTransversal, float SampleTransversal, float CenterRoughness, float SampleRoughness, vec3 CenterNormal, vec3 SampleNormal, float Kernel, const vec3 Incident) {

	float DepthWeight = pow(exp(-abs(CenterDepth - SampleDepth)), 128.0f);
	float LobeWeight = GetLobeWeight(CenterRoughness, SampleRoughness, CenterNormal, SampleNormal, Incident);
	float TransversalWeight = pow(exp(-abs(CenterTransversal - SampleTransversal)), 6.0f);
	return TransversalWeight * LobeWeight * DepthWeight;
}

float DiffuseWeight(float CenterDepth, float SampleDepth, vec3 CenterNormal, vec3 SampleNormal, float CenterLuma, float SampleLuma, float Variance, float SqrtVariance, float Kernel) {

	float LumaDistance = abs(CenterLuma - SampleLuma);

	float DepthWeight = clamp(pow(exp(-abs(CenterDepth - SampleDepth)), 128.0f), 0.0f, 1.0f);
	float NormalWeight = clamp(pow(max(dot(CenterNormal, SampleNormal), 0.0f), 8.0f), 0.0f, 1.0f);

	float PhiL = SqrtVariance;
	PhiL /= 2.0f;

	float DetailWeight = clamp(exp(-(abs(CenterLuma - SampleLuma) / PhiL)), 0.0f, 1.0f);

	float TotalWeight = DetailWeight * DepthWeight * NormalWeight * clamp(Kernel, 0.0f, 1.0f);

	return TotalWeight;
}

float DiffuseWeightBasic(float CenterDepth, float SampleDepth, vec3 CenterNormal, vec3 SampleNormal, float CenterLuma, float Kernel) {

	float DepthWeight = clamp(pow(exp(-abs(CenterDepth - SampleDepth)), 128.0f), 0.0f, 1.0f);
	float NormalWeight = clamp(pow(max(dot(CenterNormal, SampleNormal), 0.0f), 8.0f), 0.0f, 1.0f);
	float TotalWeight = DepthWeight * NormalWeight * clamp(Kernel, 0.0f, 1.0f);
	return TotalWeight;
}

void main() {

	const float Atrous[3] = float[3]( 1.0f, 2.0f / 3.0f, 1.0f / 6.0f );

	int KernelX = u_Checkered ? 1 : 2;
	int KernelY = u_Checkered ? 2 : 2;

	ivec2 Pixel = ivec2(gl_FragCoord.xy);
    ivec2 HighResPixel = Pixel * 2;

    float Depth = texelFetch(u_Depth, HighResPixel, 0).x;
	vec3 Normal = normalize(texelFetch(u_Normals, HighResPixel, 0).xyz); 
    vec3 PBR = texelFetch(u_PBR, HighResPixel, 0).xyz;
	float LinearDepth = LinearizeDepth(Depth);
	vec3 WorldPosition = WorldPosFromDepth(Depth, v_TexCoords);

	vec4 CenterSpecular = texelFetch(u_Specular, Pixel, 0);

	vec3 Incident = normalize(u_ViewerPosition - WorldPosition);

	vec4 SpecularSum = CenterSpecular;
	float TotalSpecularWeight = 1.0f;


	for (int x = -KernelX ; x <= KernelX ; x++)
	{
		for (int y = -KernelY ; y <= KernelY ; y++) 
		{
			if (x == 0 && y == 0) {
				continue; 
			}

			ivec2 SamplePixel = Pixel + ivec2(x, y) * u_StepSize;
			ivec2 SamplePixelHighRes = (Pixel + ivec2(x, y) * u_StepSize) * 2;

			float SampleDepth = LinearizeDepth(texelFetch(u_Depth, SamplePixelHighRes, 0).x);
			vec3 SampleNormal = normalize(texelFetch(u_Normals, SamplePixelHighRes, 0).xyz); 
			vec3 SamplePBR = texelFetch(u_PBR, SamplePixelHighRes, 0).xyz;

			vec4 SpecularSample = texelFetch(u_Specular, SamplePixel, 0).xyzw;

			float KernelWeight = Atrous[abs(x)] * Atrous[abs(y)];
			float SpecularWeight = SpecularWeight(LinearDepth, SampleDepth, CenterSpecular.w, SpecularSample.w, PBR.x, SamplePBR.x, Normal, SampleNormal, KernelWeight, Incident);

			SpecularSum += SpecularSample * SpecularWeight;
			TotalSpecularWeight += SpecularWeight;
		}
	}

	SpecularSum /= max(TotalSpecularWeight, 0.0000001f);

	o_Specular = SpecularSum;
}