// Spatio-temporal denoiser based on SVGF and Atrous wavelet filtering 

#version 450 core

#define PI 3.14159265359

#define saturate(x) (clamp(x, 0.0f, 1.0f))
#define square(x) (x * x)
#define cuberoot(x) (pow(x,0.3333333f))


layout (location = 0) out vec4 o_Specular; // indirect speclar 
layout (location = 1) out vec4 o_Diffuse; // xyz contains indirect diffuse, w component has large scale VXAO
layout (location = 2) out float o_Variance; // Variance of the diffuse buffer 


//layout (location = 3) out float o_Direct; 


in vec2 v_TexCoords;

uniform sampler2D u_Specular;

// Downsampled gbuffers 
uniform sampler2D u_Depth;
uniform sampler2D u_Normals;
uniform sampler2D u_PBR;

// Jitter 
uniform sampler2D u_BlueNoise;

// Frame counters 
uniform sampler2D u_SpecularFrames;

// Indirect lighting
uniform sampler2D u_Diffuse;
uniform sampler2D u_Variance;
uniform sampler2D u_TemporalUtility;

// Matrices
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
uniform int u_Pass;
uniform int u_TotalPasses;

uniform bool u_SVGF;

uniform bool u_SpatialDiffuse;
uniform bool u_SpatialSpecular;
//

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

bool InScreenSpace(in vec2 v) 
{
	return v.x > 0.0f && v.x < 1.0f && v.y > 0.0f && v.y < 1.0f;
}

float GaussianVariance(ivec2 Pixel, out float BaseVariance)
{
	vec2 TexelSize = 1.0f / textureSize(u_Variance, 0);
	float VarianceSum = 0.0f;

	const float Gaussian[2] = float[2](0.60283f, 0.198585f); // gaussian kernel
	float TotalKernel = 0.0f;

	for (int x = -1 ; x <= 1 ; x++)
	{
		for (int y = -1 ; y <= 1 ; y++)
		{
			ivec2 Coord = Pixel + ivec2(x, y);
			float Kernel = Gaussian[abs(x)] * Gaussian[abs(y)];

			float Var = texelFetch(u_Variance, Coord, 0).x;

			if (x == 0 && y == 0) { BaseVariance = Var; }

			VarianceSum += Var * Kernel;
			TotalKernel += Kernel;
		}
	}

	return VarianceSum / max(TotalKernel, 0.0000001f);
}

// Calculates a spherical gaussian whose sharpness and amplitude are similar to that of a specular lobe for a roughness value r
// Based on UE4's implementation 
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

float GetLobeWeight(in SG CenterLobe, float SampleRoughness, vec3 SampleNormal, const vec3 Incident) {
	
	const float Beta = 128.0f;

	float LobeSimilarity = 1.0f;
	float AxisSimilarity = 1.0f;

	SG SampleLobe = RoughnessLobe(SampleRoughness, SampleNormal, Incident);

	float OneOverSharpnessSum = 1.0f / (CenterLobe.Sharpness + SampleLobe.Sharpness);

	LobeSimilarity = pow(2.0f * sqrt(CenterLobe.Sharpness * SampleLobe.Sharpness) * OneOverSharpnessSum, Beta);
	AxisSimilarity = exp(-(Beta * (CenterLobe.Sharpness * SampleLobe.Sharpness) * OneOverSharpnessSum) * clamp(1.0f - dot(CenterLobe.Axis, SampleLobe.Axis), 0.0f, 1.0f));

	return LobeSimilarity * AxisSimilarity;
}

float SpecularWeight(in float CenterDepth, in float SampleDepth, in float CenterTransversal, 
					 in float SampleTransversal, in float CenterRoughness, in float SampleRoughness,
					 in vec3 CenterNormal, in vec3 SampleNormal, const vec3 Incident, 
					 in vec2 Luminances, in SG CenterSG, in float Radius, in float Frames, in float Kernel) 
{
	if (!u_SpatialSpecular) {
		return 0.0f;
	}

	// Linearize transversals 
	SampleTransversal *= 64.0f;
	CenterTransversal *= 64.0f;


	// Specular lobe weight 

	// Handles the roughness transversal and normal weight!
	float LobeWeight = GetLobeWeight(CenterRoughness, SampleRoughness, CenterNormal, SampleNormal, Incident);
	float RawLobeWeight = clamp(LobeWeight, 0.0f, 1.0f);
	LobeWeight = pow(LobeWeight, 1.5f);
	LobeWeight = clamp(LobeWeight, 0.0f, 1.0f);

	float RawTraversalWeight = pow(exp(-(abs(SampleTransversal-CenterTransversal))), 2.0f); 

	// Transversal-Luminance Weight ->

	float LuminanceWeight = 1.0f;

	bool DoLuminanceWeight = true;

	if (DoLuminanceWeight) {
		
		// If LMax is higher, it results in sharper reflections 
		float LMax = mix(24.0f, 12.0f, CenterRoughness * CenterRoughness) / 1.3f; // Weight luminance factor with roughness

		float TransversalRadius = pow(Radius, 1.0f);
		float Exponent = mix(1.0f, LMax, TransversalRadius);
		Exponent = 1.0f + Exponent;
		float LDiff = abs(Luminances.x - Luminances.y);
		float Lw = exp(-LDiff);
		LuminanceWeight = clamp(pow(Lw, Exponent), 0.0f, 1.0f);

	}

	// Frame bias (Filter more in case of disocclusions)
	float Framebias = clamp(float(Frames) / 24.0f, 0.0f, 1.0f);

	// Combine and account for the framebias 
	if (CenterRoughness > 0.25f) {
		LuminanceWeight = mix(1.0f, LuminanceWeight, Framebias);
	}

	else {
		LuminanceWeight = mix(clamp(LuminanceWeight * 3.5f, 0.142f, 1.0f), LuminanceWeight, Framebias);
	}

	float CombinedWeight = clamp(LobeWeight  * Kernel * LuminanceWeight, 0.0f, 1.0f);

	return CombinedWeight;
}

float DiffuseWeightSVGF(float CenterDepth, float SampleDepth, vec3 CenterNormal, vec3 SampleNormal, float CenterLuma, float SampleLuma, float Variance, float PhiL, float PhiD, vec2 xy, float framebias, vec2 DepthGradient, float Kernel) 
{
	float LumaDistance = abs(CenterLuma - SampleLuma);

	float DepthWeight =1.;// clamp(pow(exp(-abs(CenterDepth - SampleDepth)), 32.0f), 0.0f, 1.0f);
	//float DepthWeight = exp((-abs(CenterDepth - SampleDepth)) / (abs(dot(DepthGradient, xy)) + 0.00001f));
	float NormalWeight = clamp(pow(max(dot(CenterNormal, SampleNormal), 0.0f), 8.0f), 0.0f, 1.0f);

	float LumaError = abs(CenterLuma - SampleLuma);
	float DetailWeight = clamp(exp(-(abs(CenterLuma - SampleLuma) / PhiL)), 0.0f, 1.0f);

	DetailWeight = mix(1.0f, DetailWeight, pow(framebias, 3.0f));

	float TotalWeight = DetailWeight * DepthWeight * NormalWeight * clamp(Kernel, 0.0f, 1.0f);

	return TotalWeight;
}

float DiffuseWeightBasic(float CenterDepth, float SampleDepth, vec3 CenterNormal, vec3 SampleNormal, float Kernel) {

	float DepthWeight = 1.; //clamp(pow(exp(-abs(CenterDepth - SampleDepth)), 32.0f), 0.0f, 1.0f);
	float NormalWeight = clamp(pow(max(dot(CenterNormal, SampleNormal), 0.0f), 8.0f), 0.0f, 1.0f);
	float TotalWeight = DepthWeight * NormalWeight * clamp(Kernel, 0.0f, 1.0f);
	return TotalWeight;
}

float SVGFGetPhiL(float GaussianVariance, float RawVariance, float Frames) {

	GaussianVariance = clamp(GaussianVariance, 0.00000001f, 256.0f);
	RawVariance = clamp(RawVariance, 0.00000001f, 256.0f);

	//return sqrt(GaussianVariance);

	float f = pow(clamp(Frames / 180.0f, 0.0f, 1.0f), 4.0f);
	
	float t = exp(-pow(max(0.00000001f, GaussianVariance), 1.0f / 3.0f) * 800.0f * mix(1.0f, 1.004f, f));

	return t * mix(sqrt(GaussianVariance), 1.0f, 0.9f);
}

float Luminance(vec3 x) 
{
	return dot(x, vec3(0.2126, 0.7152,0.0722)); 
}

void main() {

	ivec2 Pixel = ivec2(gl_FragCoord.xy);
	ivec2 Dimensions = ivec2(textureSize(u_Specular, 0).xy);

	const float Atrous[3] = float[3]( 1.0f, 2.0f / 3.0f, 1.0f / 6.0f );

	const bool DoSpatial = true;

	bool CheckerStep = (Pixel.x + Pixel.y) % 2 == 0;

	int KernelX = 1;
	int KernelY = 1; //CheckerStep ? 2 : 1;

	vec4 CenterSpecular = texelFetch(u_Specular, Pixel, 0);

	if (!DoSpatial) {

		o_Specular = CenterSpecular;
		return;
	}

	// Sample GBuffers 
    float Depth = texelFetch(u_Depth, Pixel, 0).x;
	float CenterDepth = LinearizeDepth(Depth);
	vec3 Normal = normalize(texelFetch(u_Normals, Pixel, 0).xyz); 
    vec3 PBR = texelFetch(u_PBR, Pixel, 0).xyz;

	// Screenspace derivatives to get depth gradient
	float Derivative = max(dFdxFine(abs(CenterDepth)), dFdxFine(abs(CenterDepth)));

	// Calculate positions 
	float LinearDepth = LinearizeDepth(Depth);
	vec3 WorldPosition = WorldPosFromDepth(Depth, v_TexCoords);
	vec3 ViewSpaceBase = vec3(u_View * vec4(WorldPosition.xyz, 1.0f));
	float ViewLength = length(ViewSpaceBase);

	// Variance 
	float RawVar = 0.0f;
	vec4 UtilityFetch = texelFetch(u_TemporalUtility, Pixel, 0).xyzw;
	float GaussianVar = u_SVGF ? GaussianVariance(Pixel, RawVar) : 0.0f;
	float PhiL = 3. * 1. * sqrt(max(0.0f, 0.00000001f + GaussianVar));
	float PhiDepth = 0.005f * float(u_StepSize);
	float FramebiasDiffuse = clamp(UtilityFetch.z / 16.0f, 0.0f, 1.0f);

	// Specular radius 
	float ViewLengthWeight = 0.001f + (ViewLength /  64.0f);
	float SpecularHitDistance = CenterSpecular.w * 64.0f;
	float TransversalContrib = SpecularHitDistance / max((SpecularHitDistance + ViewLengthWeight), 0.00001f);
	float Rr = PBR.x;
	float SpecularRadius = clamp(pow(mix(1.0f * Rr, 1.0f, TransversalContrib), pow((1.0f-Rr),1.0/1.4f)*5.0f), 0.0f, 1.0f);
	float SpecularFrames = texture(u_SpecularFrames, v_TexCoords).x * 64.0f;

	// Eye vector 
	vec3 Incident = normalize(u_ViewerPosition - WorldPosition);

	// Calculate center roughness lobe 
	SG CenterLobe = RoughnessLobe(PBR.x, Normal, Incident);
	
	// Jitter sample pixel with blue noise to reduce aliasing patterns with the wavelet filter
	vec3 BlueNoiseSample = texelFetch(u_BlueNoise, Pixel % ivec2(256), 0).xyz;
	ivec2 JitterValue = ivec2(0);
	JitterValue = ivec2((vec2(BlueNoiseSample.xy) - 0.5f) * float(u_StepSize));

	// Total  
	vec4 SpecularSum = CenterSpecular;
	float TotalSpecularWeight = 1.0f;

	// Diffuse 
	vec4 CenterDiffuse = texelFetch(u_Diffuse, Pixel, 0).xyzw;

	float CenterDiffLuma = Luminance(CenterDiffuse.xyz);

	float CenterAOL = dot(CenterDiffuse.w, 0.333f);
	vec4 DiffuseSum = CenterDiffuse;
	float TotalDiffuseWeight = 1.0f;
	float TotalAOWeight = 1.0f;

	// Variance Sum
	float VarianceSum = 0.0f;
	
	for (int x = -KernelX ; x <= KernelX ; x++)
	{
		for (int y = -KernelY ; y <= KernelY ; y++) 
		{
			if (x == 0 && y == 0) {
				continue; 
			}

			ivec2 SamplePixel = Pixel + (ivec2(x, y) * u_StepSize) + JitterValue;

			if (SamplePixel.x < 0 || SamplePixel.x > Dimensions.x || SamplePixel.y < 0 || SamplePixel.y > Dimensions.y) {
				continue;
			}

			float KernelWeight = Atrous[abs(x)] * Atrous[abs(y)];

			// Sample gbuffers 
			float SampleDepth = LinearizeDepth(texelFetch(u_Depth, SamplePixel, 0).x);

			float DepthError = abs(SampleDepth - CenterDepth);

			float Threshold = clamp((Derivative * 6.5), 0.00175f, 0.0175f);

			if (DepthError > Threshold) {
				continue;
			}

			vec3 SampleNormal = normalize(texelFetch(u_Normals, SamplePixel, 0).xyz); 
			vec3 SamplePBR = texelFetch(u_PBR, SamplePixel, 0).xyz;

			// Sample buffers 
			vec4 SpecularSample = u_SpatialSpecular ? texelFetch(u_Specular, SamplePixel, 0).xyzw : vec4(0.);
			vec4 DiffuseSample = u_SpatialDiffuse ? texelFetch(u_Diffuse, SamplePixel, 0) : vec4(0.);
			float VarianceSample = u_SVGF ? texelFetch(u_Variance, SamplePixel, 0).x : 0.0f;

			// Calculate weights 
			float SpecularWeight = SpecularWeight(LinearDepth, SampleDepth, CenterSpecular.w, SpecularSample.w, PBR.x, SamplePBR.x, Normal, SampleNormal, Incident, vec2(Luminance(CenterSpecular.xyz), Luminance(SpecularSample.xyz)), CenterLobe, SpecularRadius, SpecularFrames, KernelWeight);
			float GBufferWeight = DiffuseWeightBasic(LinearDepth, SampleDepth, Normal, SampleNormal, KernelWeight);
			float DiffuseWeight = 0.0f; 
			
			if (u_SpatialDiffuse) {
				
				DiffuseWeight = GBufferWeight;
				float LumaError = abs(CenterDiffLuma - Luminance(DiffuseSample.xyz));
				float LumaWeight = pow(clamp(exp(-LumaError / (1.0f * PhiL + 0.0000001f)), 0.0f, 1.0f), 1.0f);
			}


			DiffuseSum.xyz += DiffuseSample.xyz * DiffuseWeight;

			// AO 
			float SampleAO = DiffuseSample.w;
			float AOWeight = 0.0f;
			
			if (u_SpatialDiffuse) {
				AOWeight = GBufferWeight * mix(1.0f, clamp(exp( -(abs(CenterAOL - dot(SampleAO, 0.333f)) * 1.0f)), 0.0f, 1.0f), clamp(FramebiasDiffuse * 1.6f, 0.0f, 1.0f));
			}

			DiffuseSum.w += SampleAO * AOWeight;
			TotalAOWeight += AOWeight;

			TotalDiffuseWeight += DiffuseWeight;

			VarianceSum += (DiffuseWeight * DiffuseWeight) * VarianceSample;

			SpecularSum += SpecularSample * SpecularWeight;
			TotalSpecularWeight += SpecularWeight;

		}
	}

	SpecularSum /= max(TotalSpecularWeight, 0.000001f);
	DiffuseSum.xyz /= max(TotalDiffuseWeight, 0.000001f);
	DiffuseSum.w /= max(TotalAOWeight, 0.000001f);
	VarianceSum /= max(TotalDiffuseWeight * TotalDiffuseWeight, 0.000001f);

	o_Specular = SpecularSum;
	o_Variance = VarianceSum;

	o_Diffuse = DiffuseSum;

	//if (GaussianVar < 0.0001f) {
	//	o_Diffuse = vec4(1.0f, 0.0f, 0.0f, 1.0f);
	//}



	
}