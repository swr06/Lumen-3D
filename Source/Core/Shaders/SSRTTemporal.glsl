#version 330 core 

layout (location = 0) out float o_AO;
layout (location = 1) out float o_DirectShadows;

in vec2 v_TexCoords;

uniform sampler2D u_Depth;
uniform sampler2D u_PreviousDepth;
uniform sampler2D u_Current;
uniform sampler2D u_History;
uniform sampler2D u_CurrentDirect;
uniform sampler2D u_HistoryDirect;
uniform sampler2D u_Normals;

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


float WaveletFilter(float Depth, vec3 Normal) {

    const float AtrousWeights[3] = float[3]( 1.0f, 2.0f / 3.0f, 1.0f / 6.0f );

    vec2 TexelSize = 1.0f / textureSize(u_Current, 0).xy;

    float Total = 0.0f;
    float TotalWeight = 0.0f;

    for (int x = -2 ; x <= 2 ; x++) {

        for (int y = -2 ; y <= 2 ; y++) {
            
            vec2 SampleCoord = v_TexCoords + vec2(x, y) * TexelSize * 1.41f;
            float KernelWeight = AtrousWeights[abs(x)] * AtrousWeights[abs(y)];

            float DepthSample = linearizeDepth(texture(u_Depth, SampleCoord).x);
            vec3 NormalSample = texture(u_Normals, SampleCoord).xyz;

            float Sample = texture(u_Current, SampleCoord).x;

            float Weight = clamp(pow(exp(-abs(DepthSample - Depth)), 32.0f), 0.0f, 1.0f) * clamp(pow(max(dot(Normal, NormalSample), 0.0f), 12.0f), 0.0f, 1.0f) * clamp(KernelWeight, 0.0f, 1.0f);

            Total += Sample * Weight;
            TotalWeight += Weight;
        }
        
    }

    Total /= max(TotalWeight, 0.000001f);

    return Total;
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

float ClampDirect(float Current, vec2 Projection, bool DoClamp) {

    float History = texture(u_HistoryDirect, Projection.xy).x;
    
    if (!DoClamp) {

        return History;

    }

    ivec2 Offsets[4] = ivec2[4](ivec2(0.0f, 1.0f), ivec2(0.0f, -1.0f), ivec2(1.0f, 0.0f), ivec2(-1.0f, 0.0f));

    float Min = 1000.0f, Max = -1000.0f;

    for (int x = 0 ; x < 4 ; x++) {

        float Sample = texelFetch(u_CurrentDirect, ivec2(gl_FragCoord.xy) + Offsets[x], 0).x;
        Min = min(Min, Sample);
        Max = max(Max, Sample);
    }

    Min = min(Min, Current) - 0.03f;
    Max = min(Max, Current) + 0.03f;


    return ClipAABBMinMax(vec3(History), vec3(Min), vec3(Max)).x;
}

void main() {

    float Depth = texture(u_Depth, v_TexCoords).x;
    float LinearizedDepth = linearizeDepth(Depth);

    vec3 Normal = texture(u_Normals, v_TexCoords).xyz;

	vec3 WorldPosition = WorldPosFromDepth(Depth, v_TexCoords);
    vec3 Reprojected = Reprojection(WorldPosition);

    float Current = WaveletFilter(LinearizedDepth, Normal); //texture(u_Current, v_TexCoords).x;
    float CurrentDirect = texture(u_CurrentDirect, v_TexCoords).x;

    float Cutoff = 0.02f;

    bool MovedCamera = distance(u_InverseView[3].xyz, u_PrevInverseView[3].xyz) > 0.002f;

    if (Reprojected.x > Cutoff && Reprojected.x < 1.0f - Cutoff && Reprojected.y > Cutoff && Reprojected.y < 1.0f - Cutoff) 
    {

        float ReprojectedDepth = linearizeDepth(texture(u_PreviousDepth, Reprojected.xy).x);
        float Error = abs(ReprojectedDepth - linearizeDepth(Reprojected.z));

        vec2 Dimensions = textureSize(u_History, 0).xy;
		vec2 Velocity = (v_TexCoords - Reprojected.xy) * Dimensions;

        float TemporalBlur = MovedCamera ? clamp(exp(-length(Velocity)) * 0.95f + 0.7725f, 0.0f, 0.95f) : 0.95f;

        float History = texture(u_History, Reprojected.xy).x;
        float HistoryDirect = ClampDirect(CurrentDirect, Reprojected.xy, MovedCamera);//texture(u_HistoryDirect, Reprojected.xy).x;

        TemporalBlur *= pow(exp(-Error), 70.0f * 1.0f);
        TemporalBlur = clamp(TemporalBlur, 0.0f, 0.95f);

        o_AO = mix(Current, History, TemporalBlur);
        o_DirectShadows = mix(CurrentDirect, HistoryDirect, clamp(TemporalBlur * 1.05f, 0.0f, 0.95f));
    }

    else {
        o_AO = Current;
        o_DirectShadows = CurrentDirect;
    }

}