#version 330 core 

layout (location = 0) out float o_Depth;

uniform sampler2D u_Depth;

void main() {
	
	ivec2 Offsets[4] = ivec2[4](ivec2(-1, -1), ivec2(-1, 1), ivec2(1, -1), ivec2(1, 1));

	float Min = 1000.0f;

	ivec2 Pixel = ivec2(gl_FragCoord.xy);
	ivec2 HighResPixel = Pixel * 2;

	for (int x = 0 ; x < 4 ; x++) {
		Min = min(Min, texelFetch(u_Depth, HighResPixel + Offsets[x], 0).x);
	}

	o_Depth = Min;
}