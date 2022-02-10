#version 330 core 

layout (location = 0) out float o_AO;
layout (location = 1) out float o_Direct;

uniform sampler2D u_Input;
uniform sampler2D u_InputDirect;
uniform sampler2D u_Depth;
uniform sampler2D u_Normals;

uniform float u_zNear;
uniform float u_zFar;
uniform int u_Frame;

const ivec2 UpscaleOffsets[4] = ivec2[](ivec2(1, 0), ivec2(-1, 0), ivec2(0, 1), ivec2(0, -1)); 

float linearizeDepth(float depth)
{
	return (2.0 * u_zNear) / (u_zFar + u_zNear - depth * (u_zFar - u_zNear));
}

void main() {

	ivec2 Pixel = ivec2(gl_FragCoord.st);

	ivec2 HighResPixel = Pixel * 2;

	// Frame checkering for temporal integration 
	bool IsCheckerStep = Pixel.x % 2 == int(Pixel.y % 2 == (u_Frame % 2));

	if (IsCheckerStep) {

		ivec2 PixelHalvedX = Pixel;
		PixelHalvedX.x /= 2;

		// Spatial upscale 

		float BaseDepth = linearizeDepth(texelFetch(u_Depth, HighResPixel, 0).x);
		vec3 BaseNormal = texelFetch(u_Normals, HighResPixel, 0).xyz;

		float TotalWeight = 0.0f;
		float TotalAO = 0.0f;
		float TotalDirect = 0.0f;

		for (int i = 0 ; i < 4 ; i++) {
			
			ivec2 Offset = UpscaleOffsets[i];
			ivec2 Coord = PixelHalvedX + Offset;
			ivec2 HighResCoord = HighResPixel + Offset;

			float SampleDepth = linearizeDepth(texelFetch(u_Depth, HighResCoord, 0).x);
			vec3 SampleNormal = texelFetch(u_Normals, HighResCoord, 0).xyz;

			float SampleAO = texelFetch(u_Input, Coord, 0).x;
			float SampleDirect = texelFetch(u_InputDirect, Coord, 0).x;

			float CurrentWeight = pow(exp(-(abs(SampleDepth - BaseDepth))), 52.0f) * pow(max(dot(SampleNormal, BaseNormal), 0.0f), 12.0f);
			CurrentWeight = clamp(CurrentWeight, 0.0000000001f, 1.0f);

			TotalAO += SampleAO * CurrentWeight;
			TotalDirect += SampleDirect * CurrentWeight;
			TotalWeight += CurrentWeight;
		}

		TotalAO /= max(TotalWeight, 0.000001f);
		TotalDirect /= max(TotalWeight, 0.000001f);
		o_AO = TotalAO;
		o_Direct = TotalDirect;
	}

	else {
		ivec2 PixelHalvedX = Pixel;
		PixelHalvedX.x /= 2;
		o_AO = texelFetch(u_Input, PixelHalvedX, 0).x;
		o_Direct = texelFetch(u_InputDirect, PixelHalvedX, 0).x;
	}

}