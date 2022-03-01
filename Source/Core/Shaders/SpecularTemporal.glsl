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

vec3 ClipAABBMinMax(vec3 prevColor, vec3 minColor, vec3 maxColor)
{
    vec3 pClip = 0.5 * (maxColor + minColor); 
    vec3 eClip = 0.5 * (maxColor - minColor); 
    vec3 vClip = prevColor - pClip;
    vec3 vUnit = vClip / eClip;
    vec3 aUnit = abs(vUnit);
    float denom = max(aUnit.x, max(aUnit.y, aUnit.z));
    return denom > 1.0 ? pClip + vClip / denom : prevColor;
}

vec3 ClipAABB(in vec3 cOld, in vec3 cNew, in vec3 centre, in vec3 halfSize)
{
    if (all(lessThanEqual(abs(cOld - centre), halfSize))) 
    {
        return cOld;
    }
    
    vec3 dir = (cNew - cOld);
    vec3 near = centre - sign(dir) * halfSize;
    vec3 tAll = (near - cOld) / dir;

    float t = 1e20;

    for (int i = 0; i < 3; i++) 
    {
        if (tAll[i] >= 0.0 && tAll[i] < t) 
        {
            t = tAll[i];
        }
    }
    
    if (t >= 1e20) 
    {
		return cOld;
    }

    return cOld + dir * t;
}

void CalculateStatistics(ivec2 Pixel, vec4 Specular, float Roughness, out vec3 Mean, out vec3 StandardDeviation, out float AverageTraversal, out float MinTransversal, out vec3 Min, out vec3 Max) {

    Mean = Specular.xyz;
    StandardDeviation = Mean * Mean;

    int KernelX = Roughness > 0.3f ? 2 : 2; 
    int KernelY = Roughness > 0.3f ? 4 : 2; 
    float TotalWeight = 1.0f;

    AverageTraversal = Specular.w;
    MinTransversal = Specular.w;

    Min = vec3(1000.0f);
    Max = vec3(-1000.0f);

    for (int x = -KernelX; x <= KernelX; x++) 
    {
        for (int y = -KernelY ; y <= KernelY ; y++) 
        {
            vec4 Sample = texelFetch(u_Specular, Pixel + ivec2(x, y), 0);
            Min = min(Sample.xyz, Min);
            Max = max(Sample.xyz, Max);
            Mean += Sample.xyz;
            StandardDeviation += Sample.xyz * Sample.xyz;
            AverageTraversal += Sample.w;
            MinTransversal = min(MinTransversal, Sample.w);
            TotalWeight += 1.0f;
        }
    }

    Mean /= TotalWeight;
    AverageTraversal /= TotalWeight;

    StandardDeviation = sqrt(StandardDeviation / TotalWeight - Mean * Mean);


}

vec3 Clip(ivec2 Pixel, bool LesserConservative, vec3 History, vec3 Specular, float Roughness, vec3 Mean, vec3 StandardDeviation, vec3 Min, vec3 Max) {

    return History;

    bool UseNewClipping = true;

    if (UseNewClipping) {

        float RoughnessWeight = Roughness < 0.1f ? 0.0f : mix(0.0f, 0.4f, clamp(pow(Roughness, 1.1f) * 1.45f, 0.0f, 1.0f)) * mix(1.0f, 0.85f, float(LesserConservative));

        vec3 VarianceClipped = ClipAABB(History, Specular, Mean, StandardDeviation + RoughnessWeight);
        vec3 StrictClipped = ClipAABBMinMax(History, Min + 0.015f, Max + 0.015f);
        vec3 StrictishClipped = clamp(VarianceClipped, Min - 0.075f, Max + 0.075f);
        vec3 Unclipped = History;

        vec3 Resolved = VarianceClipped;

        if (Roughness < 0.1f) {
            Resolved = StrictishClipped;
            if (Roughness < 0.05f) {
                Resolved = StrictClipped;
            }
        }

        return Resolved;
    }

    else {

        float Bias = mix(0.0f, 0.5f, pow(Roughness, 1.25f));
        vec3 Clipped = ClipAABBMinMax(History, Min - Bias, Max + Bias);
        return Clipped;
    }
}

void main() {

    const bool BE_USELESS = false;

    ivec2 Pixel = ivec2(gl_FragCoord.xy);

    vec2 NormalizedJitter = u_HaltonJitter / textureSize(u_Specular, 0).xy;

    float Roughness = !u_RoughSpecular ? 0.05f : texture(u_PBR, v_TexCoords).x;

    float Depth = texture(u_Depth, v_TexCoords).x;
    float LinearizedDepth = linearizeDepth(Depth);
	vec3 WorldPosition = WorldPosFromDepth(Depth, v_TexCoords);
    vec3 Normal = texture(u_Normals, v_TexCoords).xyz;

    // Current sample 
    vec4 Current = texture(u_Specular, v_TexCoords).xyzw;

    // Calculate statistics for the current pixel
    vec3 Mean, StandardDeviation, Min, Max;
    float AverageTraversal;
    float MinTransversal;
    CalculateStatistics(Pixel, Current, Roughness, Mean, StandardDeviation, AverageTraversal, MinTransversal, Min, Max);

    // Traversal
    float Transversal = MinTransversal; // Roughness < 0.1f ? Current.w : MinTransversal;
    Transversal *= 64.0f;

    // Reproject 
    vec3 I = normalize((u_ViewerPosition) - WorldPosition.xyz);
	vec3 ReflectedPosition = (WorldPosition.xyz) - I * Transversal;
    vec3 Reprojected = Reprojection(ReflectedPosition).xyz;

    // Reprojected surface 
    vec3 ReprojectedSurface = Reprojection(WorldPosition);

    bool ChangedViewMatrix = u_PrevInverseView != u_InverseView;

    // Screenspace check
    float Cutoff = ChangedViewMatrix ? 0.01f : 0.0f;
    if (Reprojected.x > Cutoff && Reprojected.x < 1.0f - Cutoff && Reprojected.y > Cutoff && Reprojected.y < 1.0f - Cutoff &&
        ReprojectedSurface.x > Cutoff && ReprojectedSurface.x < 1.0f - Cutoff && ReprojectedSurface.y > Cutoff && ReprojectedSurface.y < 1.0f - Cutoff && (!BE_USELESS)) 
    {
        

        // Depth rejection
        float ReprojectedDepth = linearizeDepth(texture(u_PreviousDepth, Reprojected.xy).x);
        float ReprojectedSurfaceDepth = linearizeDepth(texture(u_PreviousDepth, ReprojectedSurface.xy).x);

        // Calculate depth error 
        float Error = abs(ReprojectedDepth - linearizeDepth(Reprojected.z));
        float ErrorSurface = abs(ReprojectedSurfaceDepth - linearizeDepth(ReprojectedSurface.z));

        // Velocity rejection
        vec2 Dimensions = textureSize(u_HistorySpecular, 0).xy;
		vec2 Velocity = (v_TexCoords - Reprojected.xy) * Dimensions;
        bool MovedCamera = distance(u_InverseView[3].xyz, u_PrevInverseView[3].xyz) > 0.0005f;

      
        // Sample History 
        vec3 History = texture(u_HistorySpecular, Reprojected.xy).xyz;

        // Clip specular 
        vec3 ClippedSpecular = Clip(Pixel, MovedCamera, History, Current.xyz, Roughness, Mean, StandardDeviation, Min, Max);
        vec4 RawHistory = texture(u_HistorySpecular, Reprojected.xy).xyzw;
        History = MovedCamera ? ClippedSpecular : RawHistory.xyz; 

        // Calculate temporal blur
        float MovedBlurFactor = clamp(exp(-length(Velocity)) * 0.8f + 0.86f, 0.0f, 0.975f);
        
        float Frames = texture(u_Frames, Reprojected.xy).x * 64.0f;

        const float DepthWeightStrength = 1.6f;
        float DepthRejection = pow(exp(-ErrorSurface), 32.0f); 

        // Calculate temporal blur and frame increment
        float Framerejection = (DepthRejection < 0.25f ? 0.0f : (DepthRejection <= 0.5 ? 0.5f : (DepthRejection <= 0.8f ? 0.75f : 1.0f))) * mix(1.0f, MovedBlurFactor, float(MovedCamera));
        o_Frames = (Frames * Framerejection) + 1;
        float TemporalBlur = (1.0f / max(o_Frames, 1.0f));

        // Blur
        o_Color.xyz = mix(History.xyz, Current.xyz, TemporalBlur);
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