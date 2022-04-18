#version 330 core 

layout (location = 0) out vec4 o_SurfaceMotionVectors;

in vec2 v_TexCoords;

uniform sampler2D u_DepthTexture;

uniform vec2 u_CurrentJitter;

uniform mat4 u_InverseView;
uniform mat4 u_InverseProjection;
uniform mat4 u_PrevProjection;
uniform mat4 u_PrevView;
uniform mat4 u_InversePrevProjection;
uniform mat4 u_InversePrevView;

vec3 Reprojection(vec3 WorldPos) 
{
	vec4 ProjectedPosition = u_PrevProjection * u_PrevView * vec4(WorldPos, 1.0f);
	ProjectedPosition.xyz /= ProjectedPosition.w;
	ProjectedPosition.xyz = ProjectedPosition.xyz * 0.5f + 0.5f;
	return ProjectedPosition.xyz;
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

void main() {
	
    vec2 Dimensions = textureSize(u_DepthTexture, 0).xy;
    vec2 Coord = v_TexCoords + u_CurrentJitter * (1.0f / Dimensions);
    float Depth = texture(u_DepthTexture, Coord).x;
    vec3 WorldPosition = WorldPosFromDepth(Depth, v_TexCoords);
    vec3 Reprojected = Reprojection(WorldPosition);
    o_SurfaceMotionVectors = vec4(vec2(v_TexCoords.xy - Reprojected.xy), vec2(0.0f));
}