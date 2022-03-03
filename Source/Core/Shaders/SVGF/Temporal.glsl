#version 330 core 

layout (location = 0) out vec4 o_Diffuse;
layout (location = 1) out float o_Frames;

in vec2 v_TexCoords;

uniform sampler2D u_Depth;
uniform sampler2D u_PreviousDepth;
uniform sampler2D u_Normals;
uniform sampler2D u_PreviousNormals;

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

ivec2 UpscaleCoordinate(ivec2 x) {
	return x * 2;
}

ivec2 BilinearOffsets[4] = ivec2[4](ivec2(0, 0), ivec2(1, 0), ivec2(0, 1), ivec2(1, 1));

void main() {

	ivec2 Pixel = ivec2(gl_FragCoord.xy);
	ivec2 PixelHighRes = ivec2(gl_FragCoord.xy) * 2;

	vec4 CurrentSample = texture(u_Current, v_TexCoords);

    float Depth = texelFetch(u_Depth, PixelHighRes, 0).x;

	float CurrentDepth = linearizeDepth(Depth);

	vec3 Normals = texelFetch(u_Normals, PixelHighRes, 0).xyz;

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
		// Bilateral/bilinear filter ->

		vec2 Dims = textureSize(u_History, 0).xy;

		vec2 PrevPixel = (Reprojected.xy * Dims) - 0.5f;
		vec2 Interp = fract(PrevPixel);
		vec2 UVInteger = ivec2(floor(PrevPixel));

		// Raw bilinear weights 
		float BilinearWeights[4] = {
			(1.0 - Interp.x) * (1.0 - Interp.y),
			(Interp.x) * (1.0 - Interp.y),
			(1.0 - Interp.x) * (Interp.y),
			(Interp.x) * (Interp.y)
		};

		// Depth weight 
		float LinearExpectedDepth = linearizeDepth(Reprojected.z);

		// Sum 
		vec4 TemporalDiffuseSum = vec4(0.0f);
		float FrameCount = 0.0f;
		float TotalWeights = 0.0f;

		float DepthDistanceBias = abs(CurrentDepth);

		for (int Sample = 0 ; Sample < 4 ; Sample++) {

			ivec2 SampleCoord = ivec2(UVInteger + vec2(BilinearOffsets[Sample]));
			ivec2 SampleCoordHighRes = SampleCoord * 2;

			float PreviousDepth = texelFetch(u_PreviousDepth, SampleCoordHighRes, 0).x;
			float LinearPrevDepth = linearizeDepth(PreviousDepth);
			float DepthError = abs(LinearPrevDepth - LinearExpectedDepth) / DepthDistanceBias;

			vec3 PreviousNormal = texelFetch(u_PreviousNormals, SampleCoordHighRes, 0).xyz;
			float NormalDot = dot(Normals, PreviousNormal);

			if (DepthError < 0.025f && NormalDot > 0.5f) {
				float Weight = BilinearWeights[Sample];

				TemporalDiffuseSum += texelFetch(u_History, SampleCoord, 0).xyzw * Weight;
				TotalWeights += Weight;

				float Frames = texelFetch(u_Frames, SampleCoord, 0).x * 32.0f;
				FrameCount += Frames * Weight;
			}
		}

		if (TotalWeights > 0.00001f) {
			TemporalDiffuseSum /= TotalWeights;
			FrameCount /= TotalWeights;
		}

		else {
			FrameCount = 0.0f;
		}

		float SumFramesIncremented = FrameCount + 1.0f;

		float BlendFactor = 1.0f - (clamp(max(1.0f / max(SumFramesIncremented, 1.0f), 0.02f), 0.01f, 0.99f));

		o_Diffuse = mix(CurrentSample, TemporalDiffuseSum, BlendFactor);
		o_Frames = SumFramesIncremented;
	}

	else {

		o_Diffuse = CurrentSample;
		o_Frames = 0.0f;
	}

	o_Frames /= 32.0f;
	o_Frames = clamp(o_Frames, 0.0f, 256.0f);
}