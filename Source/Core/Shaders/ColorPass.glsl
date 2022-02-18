#version 400 core
#define PI 3.14159265359

layout (location = 0) out vec3 o_Color;

in vec2 v_TexCoords;

uniform sampler2D u_AlbedoTexture;
uniform sampler2D u_NormalTexture;
uniform sampler2D u_PBRTexture;
uniform sampler2D u_DepthTexture;
uniform sampler2D u_ShadowTexture;
uniform sampler2D u_BlueNoise;
uniform sampler2D u_ResolvedSpecular;
uniform sampler2D u_RTAO;
uniform sampler2D u_ScreenspaceShadows;
uniform samplerCube u_Skymap;
uniform samplerCube u_Probe;

uniform vec3 u_ViewerPosition;
uniform vec3 u_LightDirection;
uniform mat4 u_InverseView;
uniform mat4 u_InverseProjection;
uniform mat4 u_Projection;
uniform mat4 u_View;
uniform mat4 u_LightVP;

uniform float u_zNear;
uniform float u_zFar;

uniform int u_CurrentFrame;

uniform float u_RTAOStrength;

uniform vec2 u_Dims;

uniform vec3 u_ProbeCapturePoints[6];


const vec3 SUN_COLOR = vec3(6.9f, 6.9f, 10.0f);
vec3 CookTorranceBRDF(vec3 world_pos, vec3 light_dir, vec3 radiance, vec3 albedo, vec3 normal, vec2 rm, float shadow);
float FilterShadows(vec3 WorldPosition, vec3 N);
vec3 FresnelSchlickRoughness(float cosTheta, vec3 F0, float roughness);

// Spherical gaussian
struct SG {
	vec3 Axis;
	float Sharpness;
	float Amplitude;
};

// Returns a spherical gaussian approximation of the specular lobe 
// Based on the Blinn-Phong BRDF (Good enough for filtering) 
// https://www.desmos.com/calculator/rvtqpze0g7
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

float SpecularWeight(float CenterDepth, float SampleDepth, float CenterTransversal, float SampleTransversal, float CenterRoughness, float SampleRoughness, vec3 CenterNormal, vec3 SampleNormal, const vec3 Incident) {

	float DepthWeight = pow(exp(-abs(CenterDepth - SampleDepth)), 128.0f);
	float LobeWeight = GetLobeWeight(CenterRoughness, SampleRoughness, CenterNormal, SampleNormal, Incident);
	float TransversalWeight = pow(abs(CenterTransversal - SampleTransversal), 16.0f);
	return TransversalWeight * LobeWeight * DepthWeight;
}

// Gets world position from screenspace texcoord 
vec3 WorldPosFromCoord(vec2 txc)
{
	float depth = texture(u_DepthTexture, txc).r;
    float z = depth * 2.0 - 1.0;
    vec4 ClipSpacePosition = vec4(txc * 2.0 - 1.0, z, 1.0);
    vec4 ViewSpacePosition = u_InverseProjection * ClipSpacePosition;
    ViewSpacePosition /= ViewSpacePosition.w;
    vec4 WorldPos = u_InverseView * ViewSpacePosition;
    return WorldPos.xyz;
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

// Environment BRDF 
vec2 KarisEnvironmentBRDF(float NdotV, float roughness)
{
	vec4 c0 = vec4(-1.0, -0.0275, -0.572, 0.022);
	vec4 c1 = vec4(1.0, 0.0425, 1.040, -0.040);
	vec4 r = roughness * c0 + c1;
	float a004 = min(r.x * r.x, exp2(-9.28 * NdotV)) * r.x + r.y;
	return vec2(-1.04, 1.04) * a004 + r.zw;
}

float Luminance(vec3 rgb)
{
    const vec3 W = vec3(0.2125, 0.7154, 0.0721);
    return dot(rgb, W);
}

// Spatial upscaler 
void SpatialUpscale(float Depth, vec3 Normal, float Roughness, out float AO, out float ContactShadow, out vec3 Specular) {

	const float Atrous[3] = float[3]( 1.0f, 2.0f / 3.0f, 1.0f / 6.0f );

	const bool DoSpatialUpscaling = true;

	if (!DoSpatialUpscaling) {

		Specular = texture(u_ResolvedSpecular, v_TexCoords).xyz;
		AO = texture(u_RTAO, v_TexCoords).x;
		ContactShadow = texture(u_ScreenspaceShadows, v_TexCoords).x;

		return;
	}

	vec2 TexelSize = 1.0f / (u_Dims * 0.5f);

	float TotalWeight = 0.0f;
	float TotalSpecularWeight = 0.0f;

	Specular = vec3(0.0f);
	AO = 0.0f;
	ContactShadow = 0.0f;

	vec3 CenterSpecular = texture(u_ResolvedSpecular, v_TexCoords).xyz;

	float Exponent = mix(2.75f, 1.0f, clamp(pow(Roughness, 1.5f), 0.0f, 1.0f));

	for (int x = -1 ; x <= 1 ; x++) {

		for (int y = -1 ; y <= 1 ; y++) {
			
			float KernelWeight = Atrous[abs(x)] * Atrous[abs(y)];
			
			vec2 SampleCoord = v_TexCoords + vec2(x, y) * TexelSize;

			float SampleDepth = LinearizeDepth(texture(u_DepthTexture, SampleCoord).x);
			vec3 SampleNormal = texture(u_NormalTexture, SampleCoord).xyz;

			float DepthWeight = pow(exp(-abs(Depth - SampleDepth)), 128.0f);
			float NormalWeight = pow(max(dot(SampleNormal, Normal), 0.0f), 12.0f);

			float Weight = DepthWeight * NormalWeight * KernelWeight;

			vec3 SpecularSample = ivec2(x,y) == ivec2(0) ? CenterSpecular : texture(u_ResolvedSpecular, SampleCoord).xyz;
			//float LDiff = abs(Luminance(CenterSpecular) - Luminance(SpecularSample));
			//float LWeight = pow(1.0f/(1.0f+LDiff*4.5f), Exponent);
			float SpecWeight = DepthWeight * NormalWeight * KernelWeight; // * LWeight


			Specular += SpecularSample * SpecWeight;
			AO += texture(u_RTAO, SampleCoord).x * Weight;
			ContactShadow += texture(u_ScreenspaceShadows, SampleCoord).x * Weight;

			TotalWeight += Weight;
			TotalSpecularWeight += SpecWeight;
		}

	}

	TotalWeight = max(TotalWeight, 0.00001f);

	Specular /= TotalSpecularWeight;
	AO /= TotalWeight;
	//ContactShadow = texture(u_ScreenspaceShadows, v_TexCoords).x;
	ContactShadow /= TotalWeight;
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
    return u_ProbeCapturePoints[clamp(GetFaceID(Direction) ,0,5)];
}

const bool ProbeDebug = false;

void main() 
{	
	float Depth = texture(u_DepthTexture, v_TexCoords).r;
	vec3 rD = normalize(SampleIncidentRayDirection(v_TexCoords).xyz);

	// Sky check
	if (Depth > 0.999995f && (!ProbeDebug)) {
		vec3 Sample = texture(u_Skymap, rD).xyz;
		o_Color = Sample*Sample;
		return;
	}

	// Common 
	vec3 NormalizedSunDir = normalize(u_LightDirection);

	// GBuffer 
	vec3 WorldPosition = WorldPosFromDepth(Depth, v_TexCoords);
	vec3 Normal = normalize(texture(u_NormalTexture, v_TexCoords).xyz);
	vec3 Albedo = texture(u_AlbedoTexture, v_TexCoords).xyz;
	vec4 PBR = texture(u_PBRTexture, v_TexCoords).xyzw;

	// Model emission
	vec3 Emission = PBR.w * 1.5f * Albedo;

	// View vector 
    vec3 Lo = normalize(u_ViewerPosition - WorldPosition);

	// Spatial upscale ->
	float AO;
	float ScreenspaceShadow;
	vec3 SpecularIndirect;

	SpatialUpscale(LinearizeDepth(Depth), Normal, PBR.x, AO, ScreenspaceShadow, SpecularIndirect);
	SpecularIndirect = SpecularIndirect * 1.25f;

	// Direct lighting ->
	float Shadowmap = clamp(FilterShadows(WorldPosition, Normal), 0.0f, 1.0f);
	ScreenspaceShadow = ScreenspaceShadow * ScreenspaceShadow * ScreenspaceShadow * ScreenspaceShadow * ScreenspaceShadow;
	float DirectionalShadow = clamp((1.-ScreenspaceShadow) + Shadowmap, 0.0f, 1.0f); //max(pow((1.-ScreenspaceShadow), 3.0f), Shadowmap); 
	vec3 DirectLighting = CookTorranceBRDF(WorldPosition, normalize(u_LightDirection), SUN_COLOR * 0.15f, Albedo, Normal, PBR.xy, DirectionalShadow).xyz;
	
	// Indirect lighting ->
	vec3 AmbientTerm = (texture(u_Skymap, vec3(0.0f, 1.0f, 0.0f)).xyz * 0.18f) * Albedo;

	// Combine indirect diffuse and specular using environment BRDF
	
	const float TEXTURE_AO_STRENGTH = 24.0f;
	float TextureAO = pow(1.0f - PBR.z, TEXTURE_AO_STRENGTH);

	AO =  clamp(pow(AO, u_RTAOStrength * 4.0f), 0.001f, 1.0f);
	vec3 IndirectDiffuse = (AmbientTerm * AO * TextureAO);

	float Roughness = PBR.x;
	
	// Metalness is binary, so use a simple thresholding function
	float Metalness = PBR.y > 0.04f ? 1.0f : 0.0f;

	SpecularIndirect = mix(SpecularIndirect, SpecularIndirect * 1.8f, Metalness);

	// F0 
	vec3 F0 = mix(vec3(0.04f), Albedo, Metalness);

	// Fresnel roughness 
	vec3 FresnelTerm = FresnelSchlickRoughness(max(dot(Lo, Normal.xyz), 0.000001f), vec3(F0), Roughness); 
    FresnelTerm = clamp(FresnelTerm, 0.0f, 1.0f);

    vec3 kS = FresnelTerm;
    vec3 kD = 1.0f - kS;
    kD *= 1.0f - Metalness;
				
	// Samples Karis environment brdf 
    vec2 EnvironmentBRDFSampleLocation = vec2(max(dot(Lo, Normal.xyz), 0.000001f), Roughness);
    EnvironmentBRDFSampleLocation = clamp(EnvironmentBRDFSampleLocation, 0.0f, 1.0f);
    vec2 EnvironmentBRDF = KarisEnvironmentBRDF(EnvironmentBRDFSampleLocation.x, EnvironmentBRDFSampleLocation.y);

	// Integrate indirect specular lighting 
	vec3 IndirectSpecularFinal = SpecularIndirect * (FresnelTerm * EnvironmentBRDF.x + EnvironmentBRDF.y);

	// Combine indirect diffuse and indirect specular 
    vec3 IndirectLighting = (kD * IndirectDiffuse) + IndirectSpecularFinal;

	o_Color = DirectLighting + Emission + IndirectLighting;

	if (ProbeDebug) {
		vec3 sLo = -normalize(GetCapturePoint(Lo) - WorldPosition);
		int FaceID = clamp(GetFaceID(sLo),0,5);
		o_Color = vec3(pow(float(FaceID) / 5.0f, 2.0f));//texture(u_Probe, sLo).xyz;
	}

	// Nan/inf check
	if (isnan(o_Color.x) || isnan(o_Color.y) || isnan(o_Color.z) || isinf(o_Color.x) || isinf(o_Color.y) || isinf(o_Color.z)) {
        o_Color = vec3(0.0f);
    }
}

// Percentage closer filtering 
// Todo : Implement PCSS for contact hardening shadows 
float FilterShadows(vec3 WorldPosition, vec3 N)
{
	const vec2 PoissonDisk[32] = vec2[] ( vec2(-0.613392, 0.617481), vec2(0.751946, 0.453352), vec2(0.170019, -0.040254), vec2(0.078707, -0.715323), vec2(-0.299417, 0.791925), vec2(-0.075838, -0.529344), vec2(0.645680, 0.493210), vec2(0.724479, -0.580798), vec2(-0.651784, 0.717887), vec2(0.222999, -0.215125), vec2(0.421003, 0.027070), vec2(-0.467574, -0.405438), vec2(-0.817194, -0.271096), vec2(-0.248268, -0.814753), vec2(-0.705374, -0.668203), vec2(0.354411, -0.887570), vec2(0.977050, -0.108615), vec2(0.175817, 0.382366), vec2(0.063326, 0.142369), vec2(0.487472, -0.063082), vec2(0.203528, 0.214331), vec2(-0.084078, 0.898312), vec2(-0.667531, 0.326090), vec2(0.488876, -0.783441), vec2(-0.098422, -0.295755), vec2(0.470016, 0.217933), vec2(-0.885922, 0.215369), vec2(-0.696890, -0.549791), vec2(0.566637, 0.605213), vec2(-0.149693, 0.605762), vec2(0.039766, -0.396100), vec2(0.034211, 0.979980) );
	vec4 ProjectionCoordinates = u_LightVP * vec4(WorldPosition, 1.0f);
	ProjectionCoordinates.xyz = ProjectionCoordinates.xyz / ProjectionCoordinates.w; // Perspective division is not really needed for orthagonal projection but whatever
    ProjectionCoordinates.xyz = ProjectionCoordinates.xyz * 0.5f + 0.5f;
	float shadow = 0.0;

	if (ProjectionCoordinates.z > 1.0)
	{
		return 0.0f;
	}

    float ClosestDepth = texture(u_ShadowTexture, ProjectionCoordinates.xy).r; 
    float Depth = ProjectionCoordinates.z;
	float Bias = max(0.00025f * (1.0f - dot(N, u_LightDirection)), 0.0008f);  
	vec2 TexelSize = 1.0 / textureSize(u_ShadowTexture, 0);
	float noise = texture(u_BlueNoise, v_TexCoords * (u_Dims / textureSize(u_BlueNoise, 0).xy)).r;
	noise = mod(noise + 1.61803f * mod(float(u_CurrentFrame), 100.0f), 1.0f);
	float scale = 1.0f;

    int Samples = 16;

	for(int x = 0; x < Samples; x++)
	{
		float theta = noise * (2.0f*PI);
		float cosTheta = cos(theta);
		float sinTheta = sin(theta);
		mat2 dither = mat2(vec2(cosTheta, -sinTheta), vec2(sinTheta, cosTheta));
		vec2 jitter_value;
		jitter_value = PoissonDisk[x] * dither;
		float pcf = texture(u_ShadowTexture, ProjectionCoordinates.xy + (jitter_value * scale * TexelSize)).r;  // force hardware bilinear
		shadow += ProjectionCoordinates.z - Bias > pcf ? 1.0f : 0.0f;
		
		scale *= 1.06f;
	}

	shadow /= float(Samples);
    return shadow;
}





// Analytical direct lighting BRDF

float NDF(float cosLh, float roughness)
{
	float alpha   = roughness * roughness;
	float alphaSq = alpha * alpha;

	float denom = (cosLh * cosLh) * (alphaSq - 1.0) + 1.0;
	return alphaSq / (PI * denom * denom);
}

float Schlick(float cosTheta, float k)
{
	return cosTheta / (cosTheta * (1.0 - k) + k);
}

float GGX(float cosLi, float cosLo, float roughness)
{
	float r = roughness + 1.0f;
	float k = (r * r) / 8.0f; 
	return Schlick(cosLi, k) * Schlick(cosLo, k);
}

vec3 FresnelSchlick(vec3 F0, float cosTheta)
{
	return F0 + (vec3(1.0) - F0) * pow(1.0 - cosTheta, 5.0);
}

vec3 FresnelSchlickRoughness(vec3 Eye, vec3 norm, vec3 F0, float roughness) 
{
	return F0 + (max(vec3(pow(1.0f - roughness, 3.0f)) - F0, vec3(0.0f))) * pow(max(1.0 - clamp(dot(Eye, norm), 0.0f, 1.0f), 0.0f), 5.0f);
}

vec3 FresnelSchlickRoughness(float cosTheta, vec3 F0, float roughness)
{
    return F0 + (max(vec3(1.0-roughness), F0) - F0) * pow(1.0 - cosTheta, 5.0);
}

vec3 CookTorranceBRDF(vec3 world_pos, vec3 light_dir, vec3 radiance, vec3 albedo, vec3 N, vec2 rm, float shadow)
{
    const float Epsilon = 0.00001;

	// Material
	float roughness = rm.x;
	float metalness = rm.y;

    vec3 Lo = normalize(u_ViewerPosition - world_pos);
	vec3 Li = -light_dir;
	vec3 Radiance = radiance;

	// Half vector 
	vec3 Lh = normalize(Li + Lo);

	// Compute cosines 
	float CosLo = max(0.0, dot(N, Lo));
	float CosLi = max(0.0, dot(N, Li));
	float CosLh = max(0.0, dot(N, Lh));

	// Fresnel 
	vec3 F0 = mix(vec3(0.04), albedo, metalness);
	vec3 F  = FresnelSchlick(F0, max(0.0, dot(Lh, Lo)));
	
	// Distribution 
	float D = NDF(CosLh, roughness);

	// Geometry 
	float G = GGX(CosLi, CosLo, roughness);

	// Direct diffuse 
	vec3 kd = mix(vec3(1.0) - F, vec3(0.0), metalness);
	vec3 DiffuseBRDF = kd * albedo;

	// Direct specular 
	vec3 SpecularBRDF = (F * D * G) / max(Epsilon, 4.0 * CosLi * CosLo);

	// Combine 
	vec3 Combined = (DiffuseBRDF + SpecularBRDF) * Radiance * CosLi;

	// Multiply by visibility and return
	return Combined * (1.0f - shadow);
}