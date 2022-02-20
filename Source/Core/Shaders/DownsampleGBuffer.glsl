// 4x downsample (2x downsample on each axis)

#version 330 core 

layout (location = 0) out float o_Depth;
layout (location = 1) out vec4 o_Normals;
layout (location = 2) out vec4 o_PBR;

uniform sampler2D u_Depth;
uniform sampler2D u_Normals;
uniform sampler2D u_PBR;

void main() {
	
	ivec2 Offsets[4] = ivec2[4](ivec2(-1, -1), ivec2(-1, 1), ivec2(1, -1), ivec2(1, 1));

	float Min = 1000.0f;

	ivec2 Pixel = ivec2(gl_FragCoord.xy);
	ivec2 HighResPixel = Pixel * 2;

	ivec2 BestPixel = ivec2(0);

	for (int x = 0 ; x < 4 ; x++) {
		float Depth = texelFetch(u_Depth, HighResPixel + Offsets[x], 0).x;

		if (Depth < Min) {
			
			Min = min(Min, texelFetch(u_Depth, HighResPixel + Offsets[x], 0).x);
			BestPixel = Offsets[x];
		}

	}

	o_Depth = Min;

	o_Normals = texelFetch(u_Normals, HighResPixel + BestPixel, 0).xyzw;
	o_PBR = texelFetch(u_PBR, HighResPixel + BestPixel, 0).xyzw;
}