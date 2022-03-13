#version 330 core

layout (location = 0) out vec4 o_Data;

uniform sampler2D u_Input;
uniform sampler2D u_LowResDepth;
uniform float u_zNear;
uniform float u_zFar;

uniform float u_SigmaB; // Higher -> Fewer details are preserved, Lower -> More details are preserved

float linearizeDepth(float depth)
{
	return (2.0 * u_zNear) / (u_zFar + u_zNear - depth * (u_zFar - u_zNear));
}

float Luminance(vec3 rgb) {
	const vec3 W = vec3(0.2125, 0.7154, 0.0721);
    return dot(rgb, W);
}

// http://people.csail.mit.edu/sparis/bf_course/course_notes.pdf
float BilateralWeight(in vec3 v, in float SigmaB)
{
	return 0.39894f * exp(-0.5f * dot(v, v) / (SigmaB * SigmaB)) / SigmaB;
}

void main() {
	const float Atrous[3] = float[3]( 1.0f, 2.0f / 3.0f, 1.0f / 6.0f );

	ivec2 Pixel = ivec2(gl_FragCoord.xy);

	float Depth = linearizeDepth(texelFetch(u_LowResDepth, Pixel, 0).x);

	vec4 Total = texelFetch(u_Input, Pixel, 0).xyzw;
	vec4 BaseData = Total;
	float Luma = Luminance(Total.xyz);

	float TotalWeight = 1.0f;

	const int Kernel = 2;
	const float PhiL = 0.5f;
	float SigmaB = u_SigmaB; //0.8f;

	for (int x = -Kernel ; x <= Kernel ; x++) {
		for (int y = -Kernel ; y <= Kernel ; y++) {

			if (x == 0 && y == 0) { continue; }
			
			float KernelWeight = 1.0f;// Atrous[abs(x)] * Atrous[abs(y)];
			
			float DepthSample = linearizeDepth(texelFetch(u_LowResDepth, Pixel + ivec2(x,y), 0).x);
			vec4 Data = texelFetch(u_Input, Pixel + ivec2(x,y), 0).xyzw;

			float LumaAt = Luminance(Data.xyz);

			//float Lw = clamp(pow(exp(-abs(LumaAt-Luma)), PhiL),0.,1.);

			float BilateralWeight = BilateralWeight(abs(Data.xyz - BaseData.xyz), SigmaB);

			float Weight = KernelWeight * clamp((pow(exp(-(abs(DepthSample-Depth))), 16.0f)),0.,1.) * BilateralWeight;

			Total += Data * Weight;
			TotalWeight += Weight;

		}

	}

	Total /= TotalWeight;

	o_Data = Total;

}