#version 330 core 

layout (location = 0) out vec4 o_Diffuse;
layout (location = 1) out float o_Variance;

in vec2 v_TexCoords;

// 4x downsampled gbuffers (faster to sample)
uniform sampler2D u_LowResDepth;
uniform sampler2D u_LowResNormals;

uniform sampler2D u_Diffuse;
uniform sampler2D u_Utility;

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

float CalcVariance(float x, float x2) {
	return abs(x2 - x * x);
}

float Luminance(vec3 rgb)
{
    const vec3 W = vec3(0.2125, 0.7154, 0.0721);
    return dot(rgb, W);
}

const float FRAME_BIAS = 24.0f;
const int VAR_KERNEL_X = 3;
const int VAR_KERNEL_Y = 3;
const float PhiL = 2.0f;

void main() {

    ivec2 Pixel = ivec2(gl_FragCoord.xy);
    ivec2 HighResPixel = Pixel * 2;

    vec4 CenterUtilityFetch = texelFetch(u_Utility, Pixel, 0).xyzw;
    vec4 CenterDiffuse = texelFetch(u_Diffuse, Pixel, 0);

    float Frames = CenterUtilityFetch.z;

    if (Frames < FRAME_BIAS) {

        float Depth = texelFetch(u_LowResDepth, Pixel, 0).x; 
        float CenterLinDepth = linearizeDepth(Depth);
	    vec3 Normals = texelFetch(u_LowResNormals, Pixel, 0).xyz;

        float CenterL = Luminance(CenterDiffuse.xyz);

        vec4 TotalDiffuse = CenterDiffuse;
        vec4 TotalUtil = CenterUtilityFetch;
        float TotalLuma = CenterL;

        float TotalWeight = 1.0f;

        for (int x = -VAR_KERNEL_X ; x <= VAR_KERNEL_X ; x++) {

            for (int y = -VAR_KERNEL_Y ; y <= VAR_KERNEL_Y; y++) {

                if (x == 0 && y == 0) {
                    continue;
                }

                ivec2 SamplePixel = Pixel + ivec2(x, y) * 2;

                float SampleDepth = linearizeDepth(texelFetch(u_LowResDepth, SamplePixel, 0).x);
                vec3 SampleNormals = texelFetch(u_LowResNormals, SamplePixel, 0).xyz;

                vec4 SampleDiffuse = texelFetch(u_Diffuse, SamplePixel, 0);
                vec4 SampleUtility = texelFetch(u_Utility, SamplePixel, 0);

                float L = Luminance(SampleDiffuse.xyz);

                float DepthWeight = clamp(pow(exp(-abs(SampleDepth - CenterLinDepth)), 32.0f), 0.0f, 1.0f);
                float NormalWeight = clamp(pow(max(dot(SampleNormals, Normals), 0.0f), 8.0f), 0.0f, 1.0f);
                float DetailWeight = 1.0f;//clamp(exp(-(abs(L-CenterL) / PhiL)), 0.0f, 1.0f);

                const float KernelWeight = 1.0f;

                float Weight = KernelWeight * DepthWeight * NormalWeight * DetailWeight;

                TotalDiffuse += SampleDiffuse * Weight;
                TotalUtil += SampleUtility * Weight;
                TotalLuma += L * Weight;
                TotalWeight += Weight;
            }
        }

        TotalWeight = max(TotalWeight, 0.0000001f);
        TotalDiffuse /= TotalWeight;
        TotalUtil /= TotalWeight;
        TotalLuma /= TotalWeight;

        o_Diffuse = TotalDiffuse;

        float Variance = 0.0f;

        Variance = TotalUtil.y - (TotalLuma * TotalLuma);
        
        Variance = Variance * clamp((float(FRAME_BIAS) / Frames) * 1.0f, 0.0f, 10.0f);

        o_Variance = Variance;
    }

    else {
        
        o_Diffuse = CenterDiffuse;
        o_Variance = CenterUtilityFetch.w;

    }
}
