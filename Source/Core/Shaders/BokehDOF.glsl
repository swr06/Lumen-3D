#version 330 core

#define PI 3.14159265359

#define IGN GradientNoise

layout (location = 0) out vec3 o_Color;

in vec2 v_TexCoords;

uniform sampler2D u_InputTexture;
uniform sampler2D u_DepthTexture;

uniform vec2 u_Dimensions;
uniform float u_TexelScale;
uniform float u_FocalDepthTemporal;
uniform vec2 u_CameraPlanes;

uniform float u_COCScale;
uniform float u_BlurScale;
uniform float u_CAScale;
uniform bool u_CAEnabled;
uniform bool u_DOF_HQ;

// Pentagonal Bokeh Offsets
const vec2 BokehOffsets30[30] = vec2[30] ( vec2(2.0022898675353416e-17,0.32699875081197516), vec2(-0.3109942927801041,0.10104817114027898), vec2(-0.19220504324534038,-0.2645475465462665), vec2(0.1922050432453403,-0.26454754654626655), vec2(0.31099429278010404,0.10104817114027885), vec2(2.002289867535342e-16,0.6539975016239503), vec2(-0.3451202185434629,0.47501722919228284), vec2(-0.6219885855602081,0.20209634228055814), vec2(-0.558416243808115,-0.18144043630965306), vec2(-0.3844100864906809,-0.5290950930925328), vec2(-2.516695157953604e-16,-0.587153585765259), vec2(0.3844100864906805,-0.5290950930925331), vec2(0.5584162438081148,-0.18144043630965356), vec2(0.6219885855602082,0.2020963422805576), vec2(0.3451202185434625,0.47501722919228323), vec2(2.2832176666751127e-15,0.9809962524359253), vec2(-0.3667545543301486,0.8237442160223655), vec2(-0.6545102235329249,0.5893236523033282), vec2(-0.9329828783403116,0.30314451342083726), vec2(-0.8967606944609388,-0.0942533470076339), vec2(-0.7627348817937598,-0.44036518932394336), vec2(-0.5766151297360201,-0.7936426396388001), vec2(-0.18747403462167145,-0.8819959880265213), vec2(0.18311414217924682,-0.8614843067678086), vec2(0.5766151297360221,-0.7936426396387992), vec2(0.7808953690566715,-0.45085015153380154), vec2(0.8759056454813151,-0.09206139303320351), vec2(0.9329828783403129,0.3031445134208344), vec2(0.6700939143560876,0.6033552705455913), vec2(0.35822531766612176,0.804587236821626) );
const vec2 BokehOffsets50[50] = vec2[50] ( vec2(1.5017174006515064e-17, 0.24524906310898137), vec2(-0.2332457195850781, 0.07578612835520925), vec2(-0.1441537824340053, -0.1984106599096999), vec2(0.14415378243400523, -0.19841065990969994), vec2(0.23324571958507806, 0.07578612835520915), vec2(1.5017174006515063e-16, 0.49049812621796274), vec2(-0.25884016390759723, 0.35626292189421216), vec2(-0.4664914391701561, 0.15157225671041863), vec2(-0.41881218285608623, -0.13608032723223978), vec2(-0.2883075648680107, -0.39682131981939966), vec2(-1.887521368465203e-16, -0.4403651893239443), vec2(0.28830756486801035, -0.3968213198193999), vec2(0.4188121828560862, -0.13608032723224017), vec2(0.46649143917015623, 0.15157225671041819), vec2(0.25884016390759684, 0.35626292189421244), vec2(1.7124132500063345e-15, 0.7357471893269439), vec2(-0.27506591574761147, 0.6178081620167741), vec2(-0.4908826676496936, 0.4419927392274961), vec2(-0.6997371587552336, 0.22735838506562794), vec2(-0.672570520845704, -0.07069001025572542), vec2(-0.5720511613453199, -0.33027389199295754), vec2(-0.4324613473020151, -0.5952319797291), vec2(-0.1406055259662536, -0.661496991019891), vec2(0.13733560663443514, -0.6461132300758565), vec2(0.43246134730201663, -0.5952319797290995), vec2(0.5856715267925036, -0.33813761365035117), vec2(0.6569292341109864, -0.06904604477490263), vec2(0.6997371587552346, 0.22735838506562583), vec2(0.5025704357670656, 0.45251645290919346), vec2(0.26866898824959135, 0.6034404276162195), vec2(-9.617063541017877e-16, 0.9809962524359257), vec2(-0.29009284723614176, 0.8928139801909716), vec2(-0.5176803278151952, 0.7125258437884238), vec2(-0.7154325033838056, 0.5197921396256749), vec2(-0.9329828783403125, 0.303144513420836), vec2(-0.9387601734426196, 8.047547491014627e-16), vec2(-0.8376243657121732, -0.2721606544644777), vec2(-0.7154325033838066, -0.5197921396256738), vec2(-0.5766151297360206, -0.7936426396388002), vec2(-0.29009284723614365, -0.8928139801909724), vec2(7.555537024365194e-16, -0.8807303786478886), vec2(0.27327089963618906, -0.8410413489993727), vec2(0.5766151297360188, -0.793642639638801), vec2(0.7594729339574509, -0.5517893853890959), vec2(0.8376243657121717, -0.27216065446448223), vec2(0.8843232074952327, -8.663868683766059e-16), vec2(0.9329828783403126, 0.3031445134208341), vec2(0.7594729339574525, 0.5517893853890949), vec2(0.5176803278151966, 0.7125258437884228), vec2(0.27327089963619067, 0.8410413489993722) );

// https://en.wikipedia.org/wiki/Circle_of_confusion
float GetCircleOfConfusion(float Depth, float CenterDepth, float Scale) 
{
	float CircleOfConfusion = abs(Depth - CenterDepth) / 0.6f;
	return CircleOfConfusion / (1.0f / Scale + CircleOfConfusion);
}

// Interleaved gradient noise, used for dithering
float GradientNoise()
{
	vec2 coord = gl_FragCoord.xy + mod(6.0f * 100.493850275f, 500.0f);
	float noise = fract(52.9829189f * fract(0.06711056f * coord.x + 0.00583715f * coord.y));
	return noise;
}

void main() {

	vec2 PixelSizeNormalized = 1. / textureSize(u_InputTexture, 0).xy;

	float Depth = texture(u_DepthTexture, v_TexCoords).x;

	vec3 TotalColor = vec3(0.0f);
	float TotalWeight = 1.0f;

	// Rotation matrix from hash ->
	float Hash = GradientNoise();

	float Theta = Hash * 2.0f * PI;
    float CosTheta = cos(Theta);
    float SinTheta = sin(Theta);
    mat2 RotationMatrix = mat2(vec2(CosTheta, -SinTheta), vec2(SinTheta, CosTheta));

	float CoC = GetCircleOfConfusion(Depth, u_FocalDepthTemporal, u_COCScale);

	vec2 AspectCorrect = 1.0f / vec2(u_Dimensions.x / u_Dimensions.y, 1.0f);

	float DistanceToCenter = pow(distance(v_TexCoords, vec2(0.5)), 2.0f);

	bool CAEnabled = u_CAScale > 0.05f;

    vec2 ChromaticOffset = vec2(14.0f * DistanceToCenter * u_CAScale) * CoC;

	int SAMPLES = u_DOF_HQ ? 50 : 30;

	for (int Sample = 0 ; Sample < SAMPLES; Sample++) {
		
		float Weight = 1.0f;

		vec2 BokehOffset = u_DOF_HQ ? BokehOffsets50[Sample] : BokehOffsets30[Sample];
		vec2 BokehSampleOffset = RotationMatrix * BokehOffset;

		BokehSampleOffset *= AspectCorrect;

		vec2 SampleOffset = v_TexCoords + BokehSampleOffset * CoC * u_TexelScale * 8.25f * u_BlurScale;
		
		if (SampleOffset != clamp(SampleOffset, 0.0001f, 0.9999f)) {
			continue;
		}

		//TotalColor += clamp(texture(u_InputTexture, SampleOffset).xyz, 0.0f, 1.0f);

		if (u_CAEnabled) {
		
		
			vec2 ROffset = SampleOffset + ChromaticOffset;
			if (ROffset != clamp(ROffset, 0.0001f, 0.9999f)) { ROffset = SampleOffset; }
			
			vec2 BOffset = SampleOffset - ChromaticOffset;
			if (BOffset != clamp(BOffset, 0.0001f, 0.9999f)) { BOffset = SampleOffset; }
			
			vec3 CurrentSampleColor = vec3(texture(u_InputTexture, ROffset).r, texture(u_InputTexture, SampleOffset).g, texture(u_InputTexture, BOffset).b);
			TotalColor += max(CurrentSampleColor, 0.0f);
		}

		else {
			TotalColor += max(texture(u_InputTexture, SampleOffset).xyz, 0.0f);
		}

		TotalWeight += Weight;
		
	}

	TotalWeight = max(TotalWeight, 0.00001f);

	o_Color = TotalColor / TotalWeight;
}


