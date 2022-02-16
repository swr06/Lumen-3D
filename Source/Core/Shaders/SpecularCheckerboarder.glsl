#version 330 core 

layout (location = 0) out vec4 o_Data;
layout (location = 1) out float o_Transversal;

uniform sampler2D u_Input;
uniform sampler2D u_InputTransversal;

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

		float BaseDepth = linearizeDepth(texelFetch(u_Depth, HighResPixel, 0).x);
		vec3 BaseNormal = texelFetch(u_Normals, HighResPixel, 0).xyz;

		float TotalWeight = 0.0f;
		vec4 Total = vec4(0.0f);
		float TotalTraversal = 0.0f;

		for (int i = 0 ; i < 4 ; i++) {
			
			ivec2 Offset = UpscaleOffsets[i];
			ivec2 Coord = PixelHalvedX + Offset;
			ivec2 HighResCoord = HighResPixel + Offset;

			float SampleDepth = linearizeDepth(texelFetch(u_Depth, HighResCoord, 0).x);
			vec3 SampleNormal = texelFetch(u_Normals, HighResCoord, 0).xyz;

			vec4 SampleData = texelFetch(u_Input, Coord, 0).xyzw;
			float SampleTransversal = texelFetch(u_InputTransversal, Coord, 0).x;

			float CurrentWeight = pow(exp(-(abs(SampleDepth - BaseDepth))), 54.0f) * pow(max(dot(SampleNormal, BaseNormal), 0.0f), 12.0f);
			CurrentWeight = clamp(CurrentWeight, 0.0000000001f, 1.0f);

			Total += SampleData * CurrentWeight;
			TotalTraversal += SampleTransversal * CurrentWeight;
			TotalWeight += CurrentWeight;
		}

		Total /= max(TotalWeight, 0.000001f);
		TotalTraversal /= max(TotalWeight, 0.000001f);
		o_Data = Total;
		o_Transversal = TotalTraversal;
	}

	else {
		ivec2 PixelHalvedX = Pixel;
		PixelHalvedX.x /= 2;
		o_Data = texelFetch(u_Input, PixelHalvedX, 0).xyzw;
		o_Transversal = texelFetch(u_InputTransversal, PixelHalvedX, 0).x;
	}

	o_Data.w =o_Transversal;
}