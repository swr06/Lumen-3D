#version 330 core

in vec2 v_TexCoords;
layout(location = 0) out vec3 o_Color;

uniform sampler2D u_MainTexture;
uniform sampler2D u_Vol;

mat3 ACESInputMat = mat3(
    0.59719, 0.07600, 0.02840,
    0.35458, 0.90834, 0.13383,
    0.04823, 0.01566, 0.83777
);

// ODT_SAT => XYZ => D60_2_D65 => sRGB
mat3 ACESOutputMat = mat3(
    1.60475, -0.10208, -0.00327,
    -0.53108, 1.10813, -0.07276,
    -0.07367, -0.00605, 1.07602
);

vec3 RRTAndODTFit(vec3 v)
{
    vec3 a = v * (v + 0.0245786f) - 0.000090537f;
    vec3 b = v * (0.983729f * v + 0.4329510f) + 0.238081f;
    return a / b;
}

vec4 ACESFitted(vec4 Color, float Exposure)
{
    Color.rgb *= Exposure;
    
    Color.rgb = ACESInputMat * Color.rgb;
    Color.rgb = RRTAndODTFit(Color.rgb);
    Color.rgb = ACESOutputMat * Color.rgb;

    return Color;
}

float Vignette() {
	float dist = distance(vec2(0.5f), v_TexCoords) * sqrt(2.0f);
	float vig = clamp((1.0f - dist) / (1.0f - 0.05f), 0.0, 1.0);
	return clamp(vig + 0.2f, 0.0f, 1.0f);
}


void main()
{
    o_Color.xyz = (texture(u_MainTexture, v_TexCoords).xyz);
    o_Color = ACESFitted(vec4(o_Color, 1.0f), 2.0f).xyz;
	o_Color *= clamp(Vignette(), 0.0f, 1.0f);
    o_Color = clamp(o_Color, 0.0f, 1.0f);
}


