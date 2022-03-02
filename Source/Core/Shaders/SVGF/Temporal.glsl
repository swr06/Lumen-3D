#version 330 core 

layout (location = 0) out vec4 o_Diffuse;
layout (location = 1) out float o_Frames;

in vec2 v_TexCoords;

uniform sampler2D u_Depth;
uniform sampler2D u_PreviousDepth;
uniform sampler2D u_Normals;

uniform sampler2D u_History;
uniform sampler2D u_Current;

uniform sampler2D u_Frames;

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

void main() {

	vec4 Sample = texture(u_Current, v_TexCoords);

    float Depth = texture(u_Depth, v_TexCoords).x;
    vec3 WorldPosition = WorldPosFromDepth(Depth, v_TexCoords);

    vec3 Reprojected = Reprojection(WorldPosition);
	 
	vec2 PreviousCoord = Reprojected.xy;

	float bias = (u_InverseView != u_PrevInverseView) ? 0.001f : 0.0f;

	o_Frames = 1.0f;

	if (PreviousCoord.x > bias && PreviousCoord.x < 1.0f-bias &&
		PreviousCoord.y > bias && PreviousCoord.y < 1.0f-bias && 
		v_TexCoords.x > bias && v_TexCoords.x < 1.0f-bias &&
		v_TexCoords.y > bias && v_TexCoords.y < 1.0f-bias)

	{
		float PreviousDepth = texture(u_PreviousDepth, PreviousCoord.xy).x;
		float LinearPrevDepth = linearizeDepth(PreviousDepth);
		float LinearExpectedDepth = linearizeDepth(Reprojected.z);

		float Weights = 1.0f;

		float D = distance(WorldPosition, u_ViewerPosition);

		float DistanceExp = D < 128.0f ? 128.0f : 32.0f;

		Weights *= pow(exp(-abs(LinearPrevDepth - LinearExpectedDepth)), DistanceExp);

		float Frames = texture(u_Frames, Reprojected.xy).x * 32.0f;

		float SumFramesIncremented = (Frames * Weights) + 1.0f;

		float BlendFactor = 1.0f - clamp(1.0f / SumFramesIncremented, 0.0f, 1.0f);

		vec4 History = texture(u_History, Reprojected.xy);

		BlendFactor = clamp(BlendFactor, 0.01f, 0.99f);
		o_Diffuse = mix(Sample, History, BlendFactor);
		o_Frames = SumFramesIncremented;
	}

	else {

		o_Diffuse = Sample;
		o_Frames = 0.0f;
	}

	o_Frames /= 32.0f;
	o_Frames = clamp(o_Frames, 0.0f, 256.0f);
}