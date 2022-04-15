#version 400 core

const float PI = 3.1415926f;

layout (location = 0) out vec4 o_Color;

in vec2 v_TexCoords;

uniform sampler2D u_CurrentColorTexture;
uniform sampler2D u_DepthTexture;
uniform sampler2D u_PreviousColorTexture;
uniform sampler2D u_PreviousDepthTexture;

uniform vec2 u_CurrentJitter;
uniform bool u_TAAU;
uniform float u_TAAUConfidenceExponent;

uniform mat4 u_InverseView;
uniform mat4 u_InverseProjection;
uniform mat4 u_PrevProjection;
uniform mat4 u_PrevView;
uniform mat4 u_InversePrevProjection;
uniform mat4 u_InversePrevView;

uniform bool u_Enabled;

uniform float u_zNear;
uniform float u_zFar;

bool g_TAAU;

float Luminance(vec3 rgb)
{
    const vec3 W = vec3(0.2125, 0.7154, 0.0721);
    return dot(rgb, W);
}

const float Multiplier = 1.0f; // idk could be anything
vec3 Reinhard(vec3 RGB)
{
    RGB *= Multiplier;
    return vec3(RGB) / (vec3(1.0f) + Luminance(RGB));
}

vec3 InverseReinhard(vec3 RGB)
{
    return (RGB / (vec3(1.0f) - Luminance(RGB))) / Multiplier;
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

float FastDistance(in vec3 p1, in vec3 p2)
{
	return abs(p1.x - p2.x) + abs(p1.y - p2.y) + abs(p1.z - p2.z);
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

vec3 WorldPosFromDepthPrev(float depth, vec2 txc)
{
    float z = depth * 2.0 - 1.0;
    vec4 ClipSpacePosition = vec4(txc * 2.0 - 1.0, z, 1.0);
    vec4 ViewSpacePosition = u_InversePrevProjection * ClipSpacePosition;
    ViewSpacePosition /= ViewSpacePosition.w;
    vec4 WorldPos = u_InversePrevView * ViewSpacePosition;
    return WorldPos.xyz;
}

vec3 clipAABB(vec3 prevColor, vec3 minColor, vec3 maxColor)
{
    vec3 pClip = 0.5 * (maxColor + minColor); 
    vec3 eClip = 0.5 * (maxColor - minColor); 
    vec3 vClip = prevColor - pClip;
    vec3 vUnit = vClip / eClip;
    vec3 aUnit = abs(vUnit);
    float denom = max(aUnit.x, max(aUnit.y, aUnit.z));
    return denom > 1.0 ? pClip + vClip / denom : prevColor;
}

vec3 rgb2ycocg(in vec3 rgb)
{
    float co = rgb.r - rgb.b;
    float t = rgb.b + co / 2.0;
    float cg = rgb.g - t;
    float y = t + cg / 2.0;
    return vec3(y, co, cg);
}

vec3 ycocg2rgb(in vec3 ycocg)
{
    float t = ycocg.r - ycocg.b / 2.0;
    float g = ycocg.b + t;
    float b = t - ycocg.g / 2.0;
    float r = ycocg.g + b;
    return vec3(r, g, b);
}

vec4 CatmullRom(sampler2D tex, in vec2 uv)
{
    vec2 texSize = textureSize(tex, 0).xy;
    vec2 samplePos = uv * texSize;
    vec2 texPos1 = floor(samplePos - 0.5f) + 0.5f;
    vec2 f = samplePos - texPos1;
    vec2 w0 = f * (-0.5f + f * (1.0f - 0.5f * f));
    vec2 w1 = 1.0f + f * f * (-2.5f + 1.5f * f);
    vec2 w2 = f * (0.5f + f * (2.0f - 1.5f * f));
    vec2 w3 = f * f * (-0.5f + 0.5f * f);
    
    vec2 w12 = w1 + w2;
    vec2 offset12 = w2 / (w1 + w2);

    vec2 texPos0 = texPos1 - 1;
    vec2 texPos3 = texPos1 + 2;
    vec2 texPos12 = texPos1 + offset12;

    texPos0 /= texSize;
    texPos3 /= texSize;
    texPos12 /= texSize;

    vec4 result = vec4(0.0f);

    result += texture(tex, vec2(texPos0.x, texPos0.y), 0.0f) * w0.x * w0.y;
    result += texture(tex, vec2(texPos12.x, texPos0.y), 0.0f) * w12.x * w0.y;
    result += texture(tex, vec2(texPos3.x, texPos0.y), 0.0f) * w3.x * w0.y;
    result += texture(tex, vec2(texPos0.x, texPos12.y), 0.0f) * w0.x * w12.y;
    result += texture(tex, vec2(texPos12.x, texPos12.y), 0.0f) * w12.x * w12.y;
    result += texture(tex, vec2(texPos3.x, texPos12.y), 0.0f) * w3.x * w12.y;
    result += texture(tex, vec2(texPos0.x, texPos3.y), 0.0f) * w0.x * w3.y;
    result += texture(tex, vec2(texPos12.x, texPos3.y), 0.0f) * w12.x * w3.y;
    result += texture(tex, vec2(texPos3.x, texPos3.y), 0.0f) * w3.x * w3.y;

    return result;
}

float LanczosWeight(float x, const float a) 
{
	float w = a * sin(PI * x) * sin(PI * x / a) / (PI * PI * x * x);
	return x == 0.0 ? 1.0 : w;
}

float Gaussian(float DistanceSqr)
{
    return exp(-2.29f * DistanceSqr);
}

float GaussianConfidence(vec2 Difference, float OneOverStdDev2, float Scale)
{
    float ScaleSqr = Scale * Scale;
    return ScaleSqr * exp2(-0.5f * dot(Difference, Difference) * ScaleSqr * OneOverStdDev2);
}

vec3 GetNearestFragment(vec2 Txc) {
    
    ivec2 Kernel[5] = ivec2[5](ivec2(0,0), ivec2(1,0), ivec2(-1,0), ivec2(0,1), ivec2(0,-1));

    vec3 BestTexel = vec3(vec2(-1.0f), 1000.0f);

    vec2 Dimensions = textureSize(u_DepthTexture,0).xy;
    ivec2 Texel = ivec2(Txc * Dimensions - 0.5f);

    for (int x = 0 ; x < 5 ; x++) {

        ivec2 SampleTexel = Texel + Kernel[x];

        float Depth = texelFetch(u_DepthTexture, SampleTexel, 0).x;

        if (Depth < BestTexel.z) {
            BestTexel.xy = vec2(SampleTexel);
            BestTexel.z = Depth;
        }
    }

    BestTexel.xy /= Dimensions;
    return BestTexel;
}

vec4 Resampling(sampler2D tex, vec2 txc) {

    const float Alpha = 2.0f;

    vec3 Total = vec3(0.0f);
    float TotalWeight = 0.0f;
    vec2 Dimensions = textureSize(tex,0).xy;
    vec2 TexelF = txc * Dimensions - 0.5;
    ivec2 Texel = ivec2(txc * Dimensions - 0.5f);

    float MaximumWeight = -100.0f;

    const int Kernel = 2;

    for (int x = -Kernel ; x <= Kernel ; x++) {

        for (int y = -Kernel ; y <= Kernel ; y++) {

            ivec2 SampleTexel = Texel + ivec2(x, y);

            float WeightX = LanczosWeight(TexelF.x - float(SampleTexel.x), Alpha);
			float WeightY = LanczosWeight(TexelF.y - float(SampleTexel.y), Alpha);
			float Weight = WeightX * WeightY;

            vec3 Sample = texelFetch(u_CurrentColorTexture, SampleTexel, 0).xyz;

            Total += Sample * Weight;
            TotalWeight += Weight;

            MaximumWeight = max(MaximumWeight, Weight);
        }

    }

    Total /= TotalWeight;

    return vec4(Total, MaximumWeight);
}

vec4 SampleCurrent(sampler2D tex, vec2 txc) {

    return vec4(CatmullRom(u_CurrentColorTexture, v_TexCoords).xyz, 1.0f);
}

vec3 Clip(vec3 History, float ClipBoxSize, out float BlendContrastWeight) {
    
    vec3 Moments[2];
    vec3 Mean;
    float TotalWeight = 0.0f;

    int Kernel = 1;

    vec2 Dimensions = textureSize(u_CurrentColorTexture,0).xy;
    ivec2 Texel = ivec2(v_TexCoords * Dimensions - 0.5f);

    vec3 Center = texelFetch(u_CurrentColorTexture, Texel, 0).xyz;

    for (int x = -Kernel ; x <= Kernel ; x++) {

        for (int y = -Kernel ; y <= Kernel ; y++) {

            ivec2 SampleTexel = Texel + ivec2(x, y);
            vec3 Sample = (x == 0 && y == 0) ? Center : texelFetch(u_CurrentColorTexture, SampleTexel, 0).xyz;

            Sample = Reinhard(Sample);

            Moments[0] += Sample;
            Moments[1] += Sample * Sample;
            Mean += Sample;
            TotalWeight += 1.0f;
        }

    }

    Mean /= TotalWeight;
    Moments[0] /= TotalWeight;
    Moments[1] /= TotalWeight;

    vec3 StandardDeviation = sqrt(abs(Moments[1] - (Moments[0] * Moments[0])));

    float CurrentLuma = Luminance(Center);
    float HistoryLuma = Luminance(History);

    vec3 AABBMin = Moments[0] - ClipBoxSize * StandardDeviation;
    vec3 AABBMax = Moments[0] + ClipBoxSize * StandardDeviation;

    BlendContrastWeight = 1.0f / (1.0f + (max(AABBMax.x - AABBMin.y, 0) / HistoryLuma));

    return clipAABB(History, AABBMin, AABBMax);
}

float GetAABBMultiplier(float MoveLength, in vec3 Current, in vec3 History) {

    return 1.414f;

    // TODO -> Use motion vectors to detect motion and reduce the size of the aabb
    float CurrentLuma = Luminance(Current);
    float PreviousLuma = Luminance(History);
    float Contrast = clamp(abs(CurrentLuma - PreviousLuma) / max(0.2f, max(CurrentLuma, PreviousLuma)), 0.0f, 1.0f);
    float M = mix(0.95f, 1.4f, Contrast);
    return M;
}

// Todo : check this (I yolo'd it)
vec4 NearestFiltering(sampler2D Texture, vec2 UV)
{
    vec2 Dimensions = textureSize(Texture, 0).xy;
    UV *= Dimensions;
    UV = floor(UV) + 0.5f;
    UV /= Dimensions;
    return textureLod(Texture, UV, 0.0f);
}

void main()
{
    ivec2 Texel = ivec2(gl_FragCoord.xy);
    vec2 iDimensions = textureSize(u_CurrentColorTexture, 0).xy;

    g_TAAU = u_TAAU && u_Enabled;

    // Upscale 
    float JitterStrength = 0.5f;
    vec2 UpscaleCoord = g_TAAU ? (v_TexCoords + ((u_CurrentJitter * JitterStrength) / iDimensions)) : v_TexCoords;
	vec4 Resampled = g_TAAU ? Resampling(u_CurrentColorTexture, UpscaleCoord).rgba : SampleCurrent(u_CurrentColorTexture, UpscaleCoord);

    Resampled.w = g_TAAU ? Resampled.w : 1.0f;

    vec3 CurrentColor = Resampled.xyz;

	if (!u_Enabled) {
		o_Color.xyz = CurrentColor;
        o_Color.w = 0.0f;
		return;
	}

	float CurrentDepth = texture(u_DepthTexture, v_TexCoords).x;
	vec3 WorldPosition = WorldPosFromDepth(CurrentDepth, v_TexCoords).xyz;
	vec3 PreviousCoord = Reprojection(WorldPosition.xyz); 

    // Inside TAA -> Use velocity vector from the nearest fragment in a region around the current pixel 
    vec3 ClosestFragment = GetNearestFragment(v_TexCoords);
    vec3 ReprojectedClosestFrag = Reprojection(WorldPosFromDepth(ClosestFragment.z, ClosestFragment.xy));
    vec3 VelocityVector = ClosestFragment - ReprojectedClosestFrag.xyz;
    //vec2 NewPreviousCoord = v_TexCoords.xy - VelocityVector.xy;
    PreviousCoord.xy = v_TexCoords.xy - VelocityVector.xy;

	float bias = (u_InverseView != u_InversePrevView) ? 0.02f : 0.0f;

    // Screen space check 
	if (PreviousCoord.x > bias && PreviousCoord.x < 1.0f-bias &&
		PreviousCoord.y > bias && PreviousCoord.y < 1.0f-bias && 
		v_TexCoords.x > bias && v_TexCoords.x < 1.0f-bias &&
		v_TexCoords.y > bias && v_TexCoords.y < 1.0f-bias)
	{
        // Depths 
		float PreviousDepth = texture(u_PreviousDepthTexture, PreviousCoord.xy).x;
		float LinearPrevDepth = linearizeDepth(PreviousDepth);
		float LinearExpectedDepth = linearizeDepth(PreviousCoord.z);

        // Sample history 
		vec3 PrevColor = CatmullRom(u_PreviousColorTexture, PreviousCoord.xy).rgb;

        // Tonemap colors to reduce artifacts when neighbourhood clamping 
        CurrentColor = Reinhard(CurrentColor);
        PrevColor = Reinhard(PrevColor);

        // Get frame count 
        float Frames = NearestFiltering(u_PreviousColorTexture, PreviousCoord.xy).w; 
        Frames += 1.0f;

        // AABB Size 
        float ContrastBlendWeight = 1.0f;
        float AABBSize = GetAABBMultiplier(0.0f, CurrentColor, PrevColor);
        PrevColor = Clip(PrevColor, AABBSize, ContrastBlendWeight);
        ContrastBlendWeight = clamp(ContrastBlendWeight, 0.0f, 1.0f);

        // Velocity weight 
		vec2 DimensionsFull = textureSize(u_PreviousColorTexture, 0).xy;
		vec2 velocity = (v_TexCoords - PreviousCoord.xy) * DimensionsFull;
		
        // Apply velocity weight 
        float BlendFactor = exp(-length(velocity)) * 0.9f + 0.8f;
		BlendFactor = clamp(BlendFactor, 0.0f, 0.98f);
        
        // Depth weight 
		const float DepthRejectionStrength = 64.0f; // 4.0f
		float DepthRejectionWeight = pow(exp(-abs(LinearExpectedDepth-LinearPrevDepth)), DepthRejectionStrength);

        // Calculate frame blend based on accumulated frames 
        float FrameBlend = 1.0f - clamp(1.0f / Frames, 0.0f, 1.0f);
        o_Color.w = Frames * BlendFactor;

        // Confidence weight (To help with TAA-U)
        // http://behindthepixels.io/assets/files/TemporalAA.pdf
        float ConfidenceWeight = u_TAAU ? clamp(pow(Resampled.w * 1.0f, u_TAAUConfidenceExponent), 0.0f, 1.0f) : 1.0f;

        if (g_TAAU) {
            BlendFactor *= ConfidenceWeight;
        }
        
        BlendFactor *= clamp(FrameBlend, 0.0f, 1.0f);

        // Mix in tonemapped space (reduces fireflies)
		o_Color.xyz = max(InverseReinhard(mix((CurrentColor.xyz), (PrevColor.xyz), clamp(BlendFactor, 0.01f, 0.98f))), 0.0f);
        
        // Debug ->
        //o_Color.xyz = vec3(DepthRejectionWeight);
        //o_Color.xyz = vec3(distance(PreviousCoord.xy, NewPreviousCoord.xy) * 10.0f);
	}

	else 
	{
		o_Color.xyz = CurrentColor;
		o_Color.w = 0.0f;
	}
}

