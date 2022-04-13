#version 330 core
#define EPSILON 0.01f

layout (location = 0) out vec3 o_Color;

uniform sampler2D u_Texture;
uniform float u_SharpenAmount;

// Bayer dither 
float bayer2(vec2 a){
    a = floor(a);
    return fract(dot(a, vec2(0.5, a.y * 0.75)));
}
#define bayer4(a)   (bayer2(  0.5 * (a)) * 0.25 + bayer2(a))
#define bayer8(a)   (bayer4(  0.5 * (a)) * 0.25 + bayer2(a))
#define bayer16(a)  (bayer8(  0.5 * (a)) * 0.25 + bayer2(a))
#define bayer32(a)  (bayer16( 0.5 * (a)) * 0.25 + bayer2(a))
#define bayer64(a)  (bayer32( 0.5 * (a)) * 0.25 + bayer2(a))
#define bayer128(a) (bayer64( 0.5 * (a)) * 0.25 + bayer2(a))
#define bayer256(a) (bayer128(0.5 * (a)) * 0.25 + bayer2(a))

float linearToSrgb(float linear){
    float SRGBLo = linear * 12.92;
    float SRGBHi = (pow(abs(linear), 1.0/2.4) * 1.055) - 0.055;
    float SRGB = mix(SRGBHi, SRGBLo, step(linear, 0.0031308));
    return SRGB;
}

float srgbToLinear(float color) {
    float linearRGBLo = color / 12.92;
    float linearRGBHi = pow((color + 0.055) / 1.055, 2.4);
    float linearRGB = mix(linearRGBHi, linearRGBLo, step(color, 0.04045));
    return linearRGB;
}

vec3 linearToSrgb(vec3 linear) {
    vec3 SRGBLo = linear * 12.92;
    vec3 SRGBHi = (pow(abs(linear), vec3(1.0/2.4)) * 1.055) - 0.055;
    vec3 SRGB = mix(SRGBHi, SRGBLo, step(linear, vec3(0.0031308)));
    return SRGB;
}

vec3 srgbToLinear(vec3 color) {
    vec3 linearRGBLo = color / 12.92;
    vec3 linearRGBHi = pow((color + 0.055) / 1.055, vec3(2.4));
    vec3 linearRGB = mix(linearRGBHi, linearRGBLo, step(color, vec3(0.04045)));
    return linearRGB;
}

float Luminance(vec3 x) {
    return dot(x, vec3(0.2722287168, 0.6740817658, 0.0536895174));
}

float GetSat(vec3 x) { 
    return length(x); 
}

float CASWeight(vec3 x) {
    //return GetSat(x);
    //return Luminance(x);
    return max(x.g, EPSILON);
}

vec3 ContrastAdaptiveSharpening(sampler2D Texture, ivec2 Pixel, float SharpeningAmount)
{
    // Samples 
    vec3 a = texelFetch(Texture, Pixel + ivec2(0, -1), 0).rgb;
    vec3 b = texelFetch(Texture, Pixel + ivec2(-1, 0), 0).rgb;
    vec3 c = texelFetch(Texture, Pixel + ivec2(0, 0), 0).rgb;
    vec3 d = texelFetch(Texture, Pixel + ivec2(1, 0), 0).rgb;
    vec3 e = texelFetch(Texture, Pixel + ivec2(0, 1), 0).rgb;

    // Weight by luminance 
    float WeightA = CASWeight(a.xyz);
    float WeightB = CASWeight(b.xyz);
    float WeightC = CASWeight(c.xyz);
    float WeightD = CASWeight(d.xyz);
    float WeightE = CASWeight(e.xyz);

    // Calculate bounds :
    float MinWeighter = min(WeightA, min(WeightB, min(WeightC, min(WeightD, WeightE))));
    float MaxWeighter = max(WeightA, max(WeightB, max(WeightC, max(WeightD, WeightE))));

    // Apply weights :
    float FinalSharpenAmount = sqrt(min(1.0f - MaxWeighter, MinWeighter) / MaxWeighter);
    float w = FinalSharpenAmount * mix(-0.125f, -0.2f, SharpeningAmount);
    return (w * (a + b + d + e) + c) / (4.0f * w + 1.0f);
}

void BasicColorDither(inout vec3 color)
{
	const vec2 LestynCoefficients = vec2(171.0f, 231.0f);
    vec3 Lestyn = vec3(dot(LestynCoefficients, gl_FragCoord.xy));
    Lestyn = fract(Lestyn.rgb / vec3(103.0f, 71.0f, 97.0f));
    color += Lestyn.rgb / 255.0f;
}

void main() {
    ivec2 Pixel = ivec2(gl_FragCoord.xy);
    vec3 OriginalColor = texelFetch(u_Texture, Pixel, 0).xyz;

    if (false) {
        o_Color = OriginalColor;
        return;
    }

    float SharpeningAmount = u_SharpenAmount;
    vec3 SharpenedColor = SharpeningAmount < 0.001f ? OriginalColor : ContrastAdaptiveSharpening(u_Texture, Pixel, SharpeningAmount+0.02f);
    o_Color = pow(SharpenedColor, vec3(1.0f / 2.2f));
    o_Color = clamp(o_Color,0.0f,1.0f);
	BasicColorDither(o_Color);
}
