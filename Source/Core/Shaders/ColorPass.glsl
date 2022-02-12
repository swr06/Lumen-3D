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

uniform int u_CurrentFrame;

uniform float u_RTAOStrength;

uniform vec2 u_Dims;

const vec3 SUN_COLOR = vec3(6.9f, 6.9f, 10.0f);

const vec2 PoissonDisk[32] = vec2[]
(
    vec2(-0.613392, 0.617481),  vec2(0.751946, 0.453352),
    vec2(0.170019, -0.040254),  vec2(0.078707, -0.715323),
    vec2(-0.299417, 0.791925),  vec2(-0.075838, -0.529344),
    vec2(0.645680, 0.493210),   vec2(0.724479, -0.580798),
    vec2(-0.651784, 0.717887),  vec2(0.222999, -0.215125),
    vec2(0.421003, 0.027070),   vec2(-0.467574, -0.405438),
    vec2(-0.817194, -0.271096), vec2(-0.248268, -0.814753),
    vec2(-0.705374, -0.668203), vec2(0.354411, -0.887570),
    vec2(0.977050, -0.108615),  vec2(0.175817, 0.382366),
    vec2(0.063326, 0.142369),   vec2(0.487472, -0.063082),
    vec2(0.203528, 0.214331),   vec2(-0.084078, 0.898312),
    vec2(-0.667531, 0.326090),  vec2(0.488876, -0.783441),
    vec2(-0.098422, -0.295755), vec2(0.470016, 0.217933),
    vec2(-0.885922, 0.215369),  vec2(-0.696890, -0.549791),
    vec2(0.566637, 0.605213),   vec2(-0.149693, 0.605762),
    vec2(0.039766, -0.396100),  vec2(0.034211, 0.979980)
);

vec3 CookTorranceBRDF(vec3 world_pos, vec3 light_dir, vec3 radiance, vec3 albedo, vec3 normal, vec2 rm, float shadow);
float CalculateSunShadow(vec3 WorldPosition, vec3 N);
vec3 FresnelSchlickRoughness(float cosTheta, vec3 F0, float roughness);

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

vec3 GetRayDirectionAt(vec2 screenspace)
{
	vec4 clip = vec4(screenspace * 2.0f - 1.0f, -1.0, 1.0);
	vec4 eye = vec4(vec2(u_InverseProjection * clip), -1.0, 0.0);
	return vec3(u_InverseView * eye);
}

vec2 KarisEnvironmentBRDF(float NdotV, float roughness)
{
	vec4 c0 = vec4(-1.0, -0.0275, -0.572, 0.022);
	vec4 c1 = vec4(1.0, 0.0425, 1.040, -0.040);
	vec4 r = roughness * c0 + c1;
	float a004 = min(r.x * r.x, exp2(-9.28 * NdotV)) * r.x + r.y;
	return vec2(-1.04, 1.04) * a004 + r.zw;
}


void main() 
{	
	float Depth = texture(u_DepthTexture, v_TexCoords).r;
	vec3 rD = GetRayDirectionAt(v_TexCoords).xyz;

	// Sky check
	if (Depth > 0.99998f) {
		vec3 Sample = texture(u_Skymap, normalize(rD)).xyz;
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

	// Direct lighting ->
	float ScreenspaceShadow = texture(u_ScreenspaceShadows, v_TexCoords).x;
	float Shadowmap = CalculateSunShadow(WorldPosition, Normal);
	float DirectionalShadow = clamp(max(Shadowmap, 1.-clamp(pow(clamp(ScreenspaceShadow + 0.05f, 0.0f, 1.0f), 42.0f),0.0f,1.0f)), 0.0f, 1.0f);
	vec3 DirectLighting = CookTorranceBRDF(WorldPosition, normalize(u_LightDirection), SUN_COLOR * 0.15f, Albedo, Normal, PBR.xy, DirectionalShadow).xyz;
	
	// Indirect lighting ->
	float AO = texture(u_RTAO, v_TexCoords).x;
	vec3 SpecularIndirect = texture(u_ResolvedSpecular, v_TexCoords).xyz;
	vec3 AmbientTerm = (texture(u_Skymap, vec3(0.0f, 1.0f, 0.0f)).xyz * 0.18f) * Albedo;

	// Combine indirect diffuse and specular using environment BRDF
	
	const float TEXTURE_AO_STRENGTH = 24.0f;
	float TextureAO = pow(1.0f - PBR.z, TEXTURE_AO_STRENGTH);

	AO =  clamp(pow(AO, u_RTAOStrength * 4.0f), 0.001f, 1.0f);
	vec3 IndirectDiffuse = (AmbientTerm * AO * TextureAO);
	vec3 IndirectSpecular = SpecularIndirect * 3.0f;

	float Roughness = PBR.x;
	
	// Metalness is binary, so use a simple thresholding function
	float Metalness = PBR.y > 0.04f ? 1.0f : 0.0f;

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

	// Nan/inf check
	if (isnan(o_Color.x) || isnan(o_Color.y) || isnan(o_Color.z) || isinf(o_Color.x) || isinf(o_Color.y) || isinf(o_Color.z)) {
        o_Color = vec3(0.0f);
    }
}


float CalculateSunShadow(vec3 WorldPosition, vec3 N)
{
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

// -- // 

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
	vec3 Lradiance = radiance;

	// Half vector 
	vec3 Lh = normalize(Li + Lo);

	// Compute dots (Since they are used multiple times)
	float cosLo = max(0.0, dot(N, Lo));
	float cosLi = max(0.0, dot(N, Li));
	float cosLh = max(0.0, dot(N, Lh));

	// Fresnel 
	vec3 F0 = mix(vec3(0.04), albedo, metalness);
	vec3 F  = FresnelSchlick(F0, max(0.0, dot(Lh, Lo)));
	
	// Distribution 
	float D = NDF(cosLh, roughness);

	// Geometry 
	float G = GGX(cosLi, cosLo, roughness);

	// Direct diffuse 
	vec3 kd = mix(vec3(1.0) - F, vec3(0.0), metalness);
	vec3 diffuseBRDF = kd * albedo;

	// Direct specular 
	vec3 specularBRDF = (F * D * G) / max(Epsilon, 4.0 * cosLi * cosLo);

	// Combine 
	vec3 final = (diffuseBRDF + specularBRDF) * Lradiance * cosLi;

	// Multiply by visibility and return
	return final * (1.0f - shadow);
}