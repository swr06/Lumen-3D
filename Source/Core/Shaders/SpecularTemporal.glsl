#version 330 core 

layout (location = 0) out vec3 o_Color;

in vec2 v_TexCoords;

uniform sampler2D u_Depth;
uniform sampler2D u_PreviousDepth;
uniform sampler2D u_Specular;
uniform sampler2D u_HistorySpecular;
uniform sampler2D u_Normals;
uniform sampler2D u_Transversals;
uniform sampler2D u_PrevTransversals;
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

void CalculateStatistics(ivec2 Pixel, vec4 Specular, float Roughness, out vec3 Mean, out vec3 StandardDeviation, out float AverageTraversal, out vec3 Min, out vec3 Max) {

    Mean = Specular.xyz;
    StandardDeviation = Mean * Mean;

    int KernelX = 1; //Roughness < 0.1f ? 1 : 2;
    int KernelY = 2; //Roughness < 0.1f ? 1 : 2;
    float TotalWeight = 1.0f;

    AverageTraversal = Specular.w;

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
            TotalWeight += 1.0f;
        }
    }

    Mean /= TotalWeight;
    AverageTraversal /= TotalWeight;

    StandardDeviation = sqrt(StandardDeviation / TotalWeight - Mean * Mean);
}

vec3 Clip(ivec2 Pixel, vec3 History, vec3 Specular, float Roughness, vec3 Mean, vec3 StandardDeviation, vec3 Min, vec3 Max) {

    float RoughnessWeight = Roughness < 0.1f ? 0.0f : mix(0.0f, 0.45f, clamp(pow(Roughness, 1.1f) * 1.55f, 0.0f, 1.0f));
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

void main() {

    const bool BE_USELESS = false;

    ivec2 Pixel = ivec2(gl_FragCoord.xy);

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
    CalculateStatistics(Pixel, Current, Roughness, Mean, StandardDeviation, AverageTraversal, Min, Max);

    // Traversal
    float SampleTransversal = texture(u_Transversals, v_TexCoords).x;
    float Transversal = Roughness < 0.1f ? SampleTransversal : AverageTraversal;
    Transversal *= 64.0f;

    // Reproject 
    vec3 I = normalize((u_ViewerPosition) - WorldPosition.xyz);
	vec3 ReflectedPosition = (WorldPosition.xyz) - I * Transversal;
    vec3 Reprojected = Reprojection(ReflectedPosition).xyz;

    // Reprojected surface 
    vec3 ReprojectedSurface = Reprojection(WorldPosition);

    // Screenspace check
    float Cutoff = 0.009f;
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

        // Calculate temporal blur
        float TemporalBlur = MovedCamera ? clamp(exp(-length(Velocity)) * 0.825f + 0.785f, 0.0f, 0.975f) : 0.975f;

        // Sample History 
        vec3 History = texture(u_HistorySpecular, Reprojected.xy).xyz;

        // Clip specular 
        vec3 ClippedSpecular = Clip(Pixel, History, Current.xyz, Roughness, Mean, StandardDeviation, Min, Max);
        History = MovedCamera ? ClippedSpecular : texture(u_HistorySpecular, Reprojected.xy).xyz; 

        // Apply depth weight
        const float DepthWeightStrength = 1.6f;
        TemporalBlur *= pow(exp(-ErrorSurface), 128.0f);

        // Blur
        TemporalBlur = clamp(TemporalBlur, 0.0f, 0.96f);
        o_Color = mix(Current.xyz, History, TemporalBlur);
    }

    else {
        o_Color = Current.xyz;
    }

    if (isnan(o_Color.x) || isnan(o_Color.y) || isnan(o_Color.z) || isinf(o_Color.x) || isinf(o_Color.y) || isinf(o_Color.z)) {
        o_Color = vec3(0.0f);
    }

}