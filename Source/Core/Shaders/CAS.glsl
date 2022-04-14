#version 330 core
#define EPSILON 0.01f

layout (location = 0) out vec3 o_Color;

uniform sampler2D u_Texture;
uniform float u_SharpenAmount;
uniform bool u_Upscaled;

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
    int Radius = u_Upscaled ? 2 : 1;
    
    // Samples 
    vec3 a = texelFetch(Texture, Pixel + ivec2(0, -1) * Radius, 0).rgb;
    vec3 b = texelFetch(Texture, Pixel + ivec2(-1, 0) * Radius, 0).rgb;
    vec3 c = texelFetch(Texture, Pixel + ivec2(0, 0) * Radius, 0).rgb;
    vec3 d = texelFetch(Texture, Pixel + ivec2(1, 0) * Radius, 0).rgb;
    vec3 e = texelFetch(Texture, Pixel + ivec2(0, 1) * Radius, 0).rgb;

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


#define FSR_RCAS_LIMIT (0.25-(1.0/16.0))

vec4 FsrRcasLoadF(vec2 p) {
    return texture(u_Texture, p / textureSize(u_Texture, 0).xy);
}

void FsrRcasCon(
    out float con,
    float sharpness
){
    con = exp2(-sharpness);
}

vec3 FsrRcasF(vec2 ip, float con)
{
    vec2 sp = vec2(ip);
    vec3 b = FsrRcasLoadF(sp + vec2( 0,-1)).rgb;
    vec3 d = FsrRcasLoadF(sp + vec2(-1, 0)).rgb;
    vec3 e = FsrRcasLoadF(sp).rgb;
    vec3 f = FsrRcasLoadF(sp+vec2( 1, 0)).rgb;
    vec3 h = FsrRcasLoadF(sp+vec2( 0, 1)).rgb;

    float bL = b.g + .5 * (b.b + b.r);
    float dL = d.g + .5 * (d.b + d.r);
    float eL = e.g + .5 * (e.b + e.r);
    float fL = f.g + .5 * (f.b + f.r);
    float hL = h.g + .5 * (h.b + h.r);

    float nz = .25 * (bL + dL + fL + hL) - eL;
    nz=clamp(
        abs(nz)
        /(
            max(max(bL,dL),max(eL,max(fL,hL)))
            -min(min(bL,dL),min(eL,min(fL,hL)))
        ),
        0., 1.
    );

    nz=1.-.5*nz;

    vec3 mn4 = min(b, min(f, h));
    vec3 mx4 = max(b, max(f, h));
    vec2 peakC = vec2(1., -4.);
    vec3 hitMin = mn4 / (4. * mx4);
    vec3 hitMax = (peakC.x - mx4) / (4.* mn4 + peakC.y);
    vec3 lobeRGB = max(-hitMin, hitMax);
    float lobe = max(
        -FSR_RCAS_LIMIT,
        min(max(lobeRGB.r, max(lobeRGB.g, lobeRGB.b)), 0.)
    )*con;

    #ifdef FSR_RCAS_DENOISE
    lobe *= nz;
    #endif

    return (lobe * (b + d + h + f) + e) / (4. * lobe + 1.);
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

    float SharpeningAmount = u_SharpenAmount + 0.005f;
    float FCon;
    FsrRcasCon(FCon, SharpeningAmount);

    vec3 SharpenedColor = SharpeningAmount < 0.001f ? OriginalColor : 
                          ContrastAdaptiveSharpening(u_Texture, Pixel, SharpeningAmount+0.02f);
    o_Color = pow(SharpenedColor, vec3(1.0f / 2.2f));
    o_Color = clamp(o_Color,0.0f,1.0f);
	BasicColorDither(o_Color);
}
