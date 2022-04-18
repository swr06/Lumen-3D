#version 330 core

in vec2 v_TexCoords;
layout(location = 0) out vec3 o_Color;

uniform sampler2D u_MainTexture;
uniform sampler2D u_DOFTexture;
uniform sampler2D u_DepthTexture;
uniform float u_FilmGrainStrength;
uniform float u_Time;
uniform float u_ExposureMultiplier;

uniform bool u_DOF;
uniform float u_FocalDepthTemporal;
uniform float u_COCScale;



vec4 textureBicubic(sampler2D sampler, vec2 texCoords);




mat3 ACESInputMat = mat3(
    0.59719, 0.07600, 0.02840,
    0.35458, 0.90834, 0.13383,
    0.04823, 0.01566, 0.83777
);

// ODT_SAT => XYZ => D60_2_D65 => sRGB
mat3 ACESOutputMat = mat3(
    1.60475, -0.10208, -0.00327,
    -0.53108, 1.10813, -0.07276,
    -0.07367, -0.00605, 1.07602
);

vec3 RRTAndODTFit(vec3 v)
{
    vec3 a = v * (v + 0.0245786f) - 0.000090537f;
    vec3 b = v * (0.983729f * v + 0.4329510f) + 0.238081f;
    return a / b;
}

vec4 ACESFitted(vec4 Color, float Exposure)
{
    Color.rgb *= Exposure;
    
    Color.rgb = ACESInputMat * Color.rgb;
    Color.rgb = RRTAndODTFit(Color.rgb);
    Color.rgb = ACESOutputMat * Color.rgb;

    return Color;
}

float Vignette() {
	float dist = distance(vec2(0.5f), v_TexCoords) * sqrt(2.0f);
	float vig = clamp((1.0f - dist) / (1.0f - 0.05f), 0.0, 1.0);
	return clamp(vig + 0.2f, 0.0f, 1.0f);
}

float HASH2SEED = 0.0f;
vec2 hash2() 
{
	return fract(sin(vec2(HASH2SEED += 0.1, HASH2SEED += 0.1)) * vec2(43758.5453123, 22578.1459123));
}

void FilmGrain(inout vec3 oc) 
{
    if (u_FilmGrainStrength < 0.01f) {
        return;
    }

	float Strength = 0.08;
	vec3 NoiseColor = vec3(0.2001f, 0.804f, 1.02348f);
    vec3 Noise = vec3(hash2().x, hash2().y, hash2().x);
    vec3 NoiseC = Noise * exp(-oc) * NoiseColor * 0.01f;
	//oc += clamp(NoiseC, 0.0f, 1.0f);
	oc += mix(clamp(NoiseC, 0.0f, 1.0f), vec3(0.0f), 1.-u_FilmGrainStrength);
    oc *= mix(vec3(1.0f), Noise, u_FilmGrainStrength * u_FilmGrainStrength);
}

float GetCircleOfConfusion(float Depth, float CenterDepth, float Scale) 
{
	float CircleOfConfusion = abs(Depth - CenterDepth) / 0.6f;
	return CircleOfConfusion / (1.0f / Scale + CircleOfConfusion);
}

void main()
{   
    HASH2SEED = (v_TexCoords.x * v_TexCoords.y) * 489.0 * 20.0f;
	HASH2SEED += fract(u_Time) * 100.0f;

    float Depth = texture(u_DepthTexture, v_TexCoords).x;

    vec3 BaseColor = (texture(u_MainTexture, v_TexCoords).xyz);

    o_Color = BaseColor;

    if (u_DOF) {

        vec3 DOFFetch = textureBicubic(u_DOFTexture, v_TexCoords).xyz;

	    float CoC = GetCircleOfConfusion(Depth, u_FocalDepthTemporal, u_COCScale);
        float MixFactor = CoC * 2050.0f;

        MixFactor = clamp(MixFactor, 0.0f, 1.0f);

        o_Color = mix(BaseColor, DOFFetch, MixFactor);
    }

    o_Color = ACESFitted(vec4(o_Color, 1.0f), 2.0f*u_ExposureMultiplier).xyz;
	o_Color *= clamp(Vignette(), 0.0f, 1.0f);
    FilmGrain(o_Color);
    o_Color = clamp(o_Color, 0.0f, 1.0f);
}


vec4 cubic(float v){
    vec4 n = vec4(1.0, 2.0, 3.0, 4.0) - v;
    vec4 s = n * n * n;
    float x = s.x;
    float y = s.y - 4.0 * s.x;
    float z = s.z - 4.0 * s.y + 6.0 * s.x;
    float w = 6.0 - x - y - z;
    return vec4(x, y, z, w) * (1.0/6.0);
}

vec4 textureBicubic(sampler2D sampler, vec2 texCoords)
{

   vec2 texSize = textureSize(sampler, 0);
   vec2 invTexSize = 1.0 / texSize;

   texCoords = texCoords * texSize - 0.5;


    vec2 fxy = fract(texCoords);
    texCoords -= fxy;

    vec4 xcubic = cubic(fxy.x);
    vec4 ycubic = cubic(fxy.y);

    vec4 c = texCoords.xxyy + vec2 (-0.5, +1.5).xyxy;

    vec4 s = vec4(xcubic.xz + xcubic.yw, ycubic.xz + ycubic.yw);
    vec4 offset = c + vec4 (xcubic.yw, ycubic.yw) / s;

    offset *= invTexSize.xxyy;

    vec4 sample0 = texture(sampler, offset.xz);
    vec4 sample1 = texture(sampler, offset.yz);
    vec4 sample2 = texture(sampler, offset.xw);
    vec4 sample3 = texture(sampler, offset.yw);

    float sx = s.x / (s.x + s.y);
    float sy = s.z / (s.z + s.w);

    return mix(
       mix(sample3, sample2, sx), mix(sample1, sample0, sx)
    , sy);
}