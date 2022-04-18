#version 330 core 


layout (location = 0) out vec4 o_Color;
layout (location = 1) out float o_Frames;

in vec2 v_TexCoords;

uniform sampler2D u_Depth;
uniform sampler2D u_PreviousDepth;
uniform sampler2D u_Specular;
uniform sampler2D u_HistorySpecular;
uniform sampler2D u_Normals;
uniform sampler2D u_PrevTransversals;
uniform sampler2D u_PBR;
uniform sampler2D u_HitMask;


uniform sampler2D u_Frames;

uniform mat4 u_Projection;
uniform mat4 u_View;
uniform mat4 u_InverseView;
uniform mat4 u_InverseProjection;

uniform vec2 u_HaltonJitter;

uniform mat4 u_PrevView;
uniform mat4 u_PrevProjection;
uniform mat4 u_PrevInverseView;
uniform mat4 u_PrevInverseProjection;

uniform vec3 u_ViewerPosition;

uniform float u_zNear;
uniform float u_zFar;

uniform bool u_RoughSpecular;
uniform bool u_SpecularChecker;
uniform bool u_SpecularTemporal;

float remap(float x, float a, float b, float c, float d)
{
    return (((x - a) / (b - a)) * (d - c)) + c;
}

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

vec3 WorldPosFromDepth(float depth, vec2 txc)
{
    float z = depth * 2.0 - 1.0;
    vec4 ClipSpacePosition = vec4(txc * 2.0 - 1.0, z, 1.0);
    vec4 ViewSpacePosition = u_InverseProjection * ClipSpacePosition;
    ViewSpacePosition /= ViewSpacePosition.w;
    vec4 WorldPos = u_InverseView * ViewSpacePosition;
    return WorldPos.xyz;
}

vec3 Reprojection(vec3 WorldPosition) 
{
	vec4 ProjectedPosition = u_PrevProjection * u_PrevView * vec4(WorldPosition, 1.0f);
	ProjectedPosition.xyz /= ProjectedPosition.w;
	ProjectedPosition.xyz = ProjectedPosition.xyz * 0.5f + 0.5f;
	return ProjectedPosition.xyz;
}

float DistanceSqr(vec3 A, vec3 B)
{
    vec3 C = A - B;
    return dot(C, C);
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

vec4 Resampler(sampler2D tex, in vec2 uv)
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

bool GetDisocclusion(vec3 History, vec3 Specular, vec2 PBR) {

    ivec2 RawPixel = ivec2(gl_FragCoord.xy);
    ivec2 Downscaled = RawPixel / (u_SpecularChecker ? ivec2(2, 1) : ivec2(1));
    
    // Mask is 1 when there was a hit and 0 when there wasnt
    float Mask = texelFetch(u_HitMask,Downscaled,0).x; 
    bool MaskCheck = Mask < 0.1f; 

    return MaskCheck;
}

vec2 GetMinTransversals(vec2 Reprojected) {

    ivec2 Kernel[5] = ivec2[5](ivec2(0,0), ivec2(1,0), ivec2(-1,0), ivec2(0,1), ivec2(0,-1));

    vec2 Dimensions = textureSize(u_Specular, 0).xy;
    ivec2 Pixel = ivec2(gl_FragCoord.xy);
    ivec2 ReprojectedPixel = ivec2(Reprojected * Dimensions - 0.5f);

    vec2 MinTransversals = vec2(10000.0f);
    
    for (int x = 0 ; x < 5 ; x++) {

        ivec2 Offset = Kernel[x];
        float Transversal = texelFetch(u_Specular, Pixel + Offset, 0).w;
        float PrevTransversal = texelFetch(u_PrevTransversals, (ReprojectedPixel + Offset) / (u_SpecularChecker ? ivec2(2, 1) : ivec2(1)), 0).w;
        Transversal *= 64.0f;
        PrevTransversal *= 64.0f;
        MinTransversals.x = min(MinTransversals.x, Transversal);
        MinTransversals.y = min(MinTransversals.y, PrevTransversal);
    }

    return MinTransversals;
}

// Clips specular (reduces ghosting)
vec3 Clip(ivec2 Pixel, bool LesserConservative, vec3 History, vec3 Specular, float Roughness, float Metalness, bool DisocclusionBias, float TransversalError) {
    
    float ClipBoxSize = Roughness < 0.26f ? mix(1.5f, 2.2f, clamp(remap(Roughness, 0.0f, 0.26, 0.0f, 1.0f), 0.0f, 1.0f))
                        : mix(3.0f, 8.0f, clamp(pow(Roughness, 0.5f), 0.0f, 1.0f));
    
    ClipBoxSize *= mix(0.95f, 2.0f, float(TransversalError < mix(3.5f, 6.5f, Roughness)));
    bool Metallic = Metalness > 0.05f;

    if (Metallic) {

        if (Roughness < 0.51f) {
            ClipBoxSize *= 0.8f;
            ClipBoxSize = max(ClipBoxSize, 1.4f);
        }

        else {
            ClipBoxSize *= 0.95f;
            ClipBoxSize = max(ClipBoxSize, 1.4f);
        }
    }

    // Increase size of the AABB when no hit was found..
    ClipBoxSize *= mix(1.0f, 1.55f, float(DisocclusionBias && (Roughness > 0.2f)));

    // Calculate variance using kernel 
    vec3 Moments[2];
    vec3 Mean;
    float TotalWeight = 0.0f;

    int Kernel = Roughness < 0.25f ? 1 : 2;

    vec2 Dimensions = textureSize(u_Specular,0).xy;
    ivec2 Texel = ivec2(v_TexCoords * Dimensions - 0.5f);

    vec3 Center = texelFetch(u_Specular, Texel, 0).xyz;

    for (int x = -Kernel ; x <= Kernel ; x++) {

        for (int y = -Kernel ; y <= Kernel ; y++) {

            ivec2 SampleTexel = Texel + ivec2(x, y);
            vec3 Sample = (x == 0 && y == 0) ? Center : texelFetch(u_Specular, SampleTexel, 0).xyz;

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

    // Calculate std deviation 
    vec3 StandardDeviation = sqrt(abs(Moments[1] - (Moments[0] * Moments[0])));

    float CurrentLuma = Luminance(Center);
    float HistoryLuma = Luminance(History);

    vec3 AABBMin = Moments[0] - ClipBoxSize * StandardDeviation;
    vec3 AABBMax = Moments[0] + ClipBoxSize * StandardDeviation;

    vec3 Clipped = clipAABB(History, AABBMin, AABBMax);
    return Clipped;
}


vec3 ReprojectReflection(vec3 P, vec3 Incident, float Transversal) {

    vec3 ReflectedPosition = (P.xyz) - Incident * Transversal;
    vec3 Reprojected = Reprojection(ReflectedPosition).xyz;

    return Reprojected;
}

void main() {

    const bool BE_USELESS = false;

    ivec2 Pixel = ivec2(gl_FragCoord.xy);

    vec2 NormalizedJitter = u_HaltonJitter / textureSize(u_Specular, 0).xy;

    vec3 PBRFetch = texture(u_PBR, v_TexCoords).xyz;
    float Roughness = !u_RoughSpecular ? 0.05f : PBRFetch.x;

    float Depth = texture(u_Depth, v_TexCoords).x;
    float LinearizedDepth = linearizeDepth(Depth);
	vec3 WorldPosition = WorldPosFromDepth(Depth, v_TexCoords);
    vec3 Normal = texture(u_Normals, v_TexCoords).xyz;

    // Current sample 
    vec4 Current = Resampler(u_Specular, v_TexCoords).xyzw;

    if (!u_SpecularTemporal) {
        o_Color = Current;
        return;
    }

    // Traversal
    float Transversal = Current.w;
    Transversal *= 64.0f;

    // Reproject 
    vec3 I = normalize((u_ViewerPosition) - WorldPosition.xyz);
	vec3 Reprojected = ReprojectReflection(WorldPosition, I, Transversal);

    // Reprojected surface 
    vec3 ReprojectedSurface = Reprojection(WorldPosition);

    bool ChangedViewMatrix = u_PrevInverseView != u_InverseView;

    // Screenspace check
    float Cutoff = ChangedViewMatrix ? 0.02f : 0.0f;
    if (Reprojected.x > Cutoff && Reprojected.x < 1.0f - Cutoff && Reprojected.y > Cutoff && Reprojected.y < 1.0f - Cutoff &&
        ReprojectedSurface.x > Cutoff && ReprojectedSurface.x < 1.0f - Cutoff && ReprojectedSurface.y > Cutoff && ReprojectedSurface.y < 1.0f - Cutoff && (!BE_USELESS)) 
    {
        // Sample Prev Transversal 
        vec2 MinTransversals = GetMinTransversals(Reprojected.xy);
        float TransversalError = abs(MinTransversals.x - MinTransversals.y);
    
        // Depth rejection
        float ReprojectedDepth = linearizeDepth(texture(u_PreviousDepth, Reprojected.xy).x);
        float ReprojectedSurfaceDepth = linearizeDepth(texture(u_PreviousDepth, ReprojectedSurface.xy).x);

        // Calculate depth error 
        float Error = abs(ReprojectedDepth - linearizeDepth(Reprojected.z));
        float ErrorSurface = abs(ReprojectedSurfaceDepth - linearizeDepth(ReprojectedSurface.z));

        // Velocity rejection
        vec2 Dimensions = textureSize(u_HistorySpecular, 0).xy;
		vec2 Velocity = (v_TexCoords - Reprojected.xy) * Dimensions;
      
        // Sample History 
        vec3 History = Resampler(u_HistorySpecular, Reprojected.xy).xyz;

        // Reinhard tonemapping (massively reduces fireflies, increases temporal stability)
        History = Reinhard(History);
        Current.xyz = Reinhard(Current.xyz);

        bool DisocclusionDetected = GetDisocclusion(History.xyz, Current.xyz, PBRFetch.xy);

        // Clip specular 
        vec3 ClippedSpecular = Clip(Pixel, true, History, Current.xyz, Roughness, PBRFetch.y, DisocclusionDetected, TransversalError);
        History = ClippedSpecular; 

        // Calculate temporal blur
        float MovedBlurFactor = clamp(exp(-length(Velocity)) * 0.8f + 0.875f, 0.0f, 0.975f);
        
        const float DepthWeightStrength = 1.6f;
        float DepthRejection = pow(exp(-ErrorSurface), 64.0f); 

        // Calculate temporal blur and frame increment
        float Frames = texture(u_Frames, Reprojected.xy).x * 64.0f;
        Frames += 1.0f;

        float Framerejection = DepthRejection;

        float TemporalBlur = (1.0f / max(Frames, 1.0f)) * Framerejection;
        TemporalBlur = clamp(TemporalBlur, 0.02f, 0.98f);

        o_Frames = (Frames * Framerejection);

        // Blur
        o_Color.xyz = InverseReinhard(mix(History.xyz, Current.xyz, TemporalBlur));

        o_Color.xyz = max(o_Color.xyz, 0.0f);

        vec4 RawHistory = texture(u_HistorySpecular, Reprojected.xy).xyzw;
        o_Color.w = mix(Transversal / 64.0f, RawHistory.w, TemporalBlur);
    }

    else {
        o_Color.xyz = Current.xyz;
        o_Color.w = Transversal / 64.0f;
        o_Frames = 1.0f;
    }

    if (isnan(o_Color.x) || isnan(o_Color.y) || isnan(o_Color.z) || isinf(o_Color.x) || isinf(o_Color.y) || isinf(o_Color.z)) {
        o_Color.xyz = vec3(0.0f);
        o_Color.w = vec3(0.0f).x;
    }

    o_Frames = clamp(o_Frames / 64.0f, 0.0f, 1.0f);
}