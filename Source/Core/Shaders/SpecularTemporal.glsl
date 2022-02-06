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

vec2 Reprojection(vec3 WorldPosition) 
{
	vec4 ProjectedPosition = u_PrevProjection * u_PrevView * vec4(WorldPosition, 1.0f);
	ProjectedPosition.xyz /= ProjectedPosition.w;
	ProjectedPosition.xy = ProjectedPosition.xy * 0.5f + 0.5f;
	return ProjectedPosition.xy;
}

float DistanceSqr(vec3 A, vec3 B)
{
    vec3 C = A - B;
    return dot(C, C);
}

void main() {

    float Depth = texture(u_Depth, v_TexCoords).x;
    float LinearizedDepth = linearizeDepth(Depth);

    vec3 Normal = texture(u_Normals, v_TexCoords).xyz;

	vec3 WorldPosition = WorldPosFromDepth(Depth, v_TexCoords);
    float Transversal = texture(u_Transversals, v_TexCoords).x * 64.0f;
    vec3 I = normalize((u_ViewerPosition) - WorldPosition.xyz);
	vec3 ReflectedPosition = (WorldPosition.xyz) - I * Transversal;
    vec2 Reprojected = Reprojection(ReflectedPosition);
    vec2 ReprojectedSurface = Reprojection(WorldPosition);

    vec3 Current = texture(u_Specular, v_TexCoords).xyz;

    if (Reprojected == clamp(Reprojected, 0.002f, 0.998f) && ReprojectedSurface == clamp(ReprojectedSurface, 0.002f, 0.998f)) {

        float ReprojectedDepth = linearizeDepth(texture(u_PreviousDepth, Reprojected.xy).x);
        float ReprojectedSurfaceDepth = linearizeDepth(texture(u_PreviousDepth, ReprojectedSurface.xy).x);
        float Error = abs(ReprojectedDepth - LinearizedDepth);
        float ErrorSurface = abs(ReprojectedSurfaceDepth - LinearizedDepth);

        float PreviousT = texture(u_PrevTransversals, Reprojected).x * 64.0f;

        vec2 Dimensions = textureSize(u_HistorySpecular, 0).xy;
		vec2 Velocity = (v_TexCoords - Reprojected.xy) * Dimensions;

        float TemporalBlur = clamp(exp(-length(Velocity)) * 0.95f + 0.88f, 0.0f, 0.97f);

        vec3 History = texture(u_HistorySpecular, Reprojected.xy).xyz;

        const float DepthWeightStrength = 1.6f;
        TemporalBlur *= pow(exp(-Error), 24.0f * DepthWeightStrength) * pow(exp(-ErrorSurface), 32.0f * DepthWeightStrength);

        TemporalBlur = clamp(TemporalBlur, 0.0f, 0.96f);
        o_Color = mix(Current, History, TemporalBlur);
    }

    else {
        o_Color = Current;
    }

    if (isnan(o_Color.x) || isnan(o_Color.y) || isnan(o_Color.z) || isinf(o_Color.x) || isinf(o_Color.y) || isinf(o_Color.z)) {
        o_Color = vec3(0.0f);
    }
}