#version 330 core

const float PI = 3.1415926535897932384626433832795;

layout (location = 0) out vec3 o_Color;

in vec2 v_TexCoords;

uniform sampler2D u_CurrentColorTexture;
uniform sampler2D u_DepthTexture;
uniform sampler2D u_PreviousColorTexture;
uniform sampler2D u_PreviousDepthTexture;

uniform vec2 u_CurrentJitter;

uniform mat4 u_InverseView;
uniform mat4 u_InverseProjection;
uniform mat4 u_PrevProjection;
uniform mat4 u_PrevView;
uniform mat4 u_InversePrevProjection;
uniform mat4 u_InversePrevView;

uniform bool u_Enabled;

uniform float u_zNear;
uniform float u_zFar;

vec3 min4(vec3 a, vec3 b, vec3 c, vec3 d)
{
    return min(a, min(b, min(c, d)));
}

vec3 max4(vec3 a, vec3 b, vec3 c, vec3 d)
{
    return max(a, max(b, max(c, d)));
}

float linearizeDepth(float depth)
{
	return (2.0 * u_zNear) / (u_zFar + u_zNear - depth * (u_zFar - u_zNear));
}


vec3 Reprojection(vec3 WorldPos) 
{
	vec4 ProjectedPosition = u_PrevProjection * u_PrevView * vec4(WorldPos, 1.0f);
	ProjectedPosition.xyz /= ProjectedPosition.w;
	ProjectedPosition.xyz = ProjectedPosition.xyz * 0.5f + 0.5f;
	return ProjectedPosition.xyz;
}

float FastDistance(in vec3 p1, in vec3 p2)
{
	return abs(p1.x - p2.x) + abs(p1.y - p2.y) + abs(p1.z - p2.z);
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

vec3 WorldPosFromDepthPrev(float depth, vec2 txc)
{
    float z = depth * 2.0 - 1.0;
    vec4 ClipSpacePosition = vec4(txc * 2.0 - 1.0, z, 1.0);
    vec4 ViewSpacePosition = u_InversePrevProjection * ClipSpacePosition;
    ViewSpacePosition /= ViewSpacePosition.w;
    vec4 WorldPos = u_InversePrevView * ViewSpacePosition;
    return WorldPos.xyz;
}

vec3 clipAABB(vec3 prevColor, vec3 minColor, vec3 maxColor)
{
    vec3 pClip = 0.5 * (maxColor + minColor); 
    vec3 eClip = 0.5 * (maxColor - minColor); 
    vec3 vClip = prevColor - pClip;
    vec3 vUnit = vClip / eClip;
    vec3 aUnit = abs(vUnit);
    float denom = max(aUnit.x, max(aUnit.y, aUnit.z));
    return denom > 1.0 ? pClip + vClip / denom : prevColor;
}

float GetHistoryContrast(float historyLuma, float minNeighbourLuma, float maxNeighbourLuma)
{
    float lumaContrast = max(maxNeighbourLuma - minNeighbourLuma, 0) / historyLuma;
    float blendFactor = 1.0 / 20.0f;
    return clamp(blendFactor * (1.0f / (1.0 + lumaContrast)), 0.0f, 1.0f);
}

vec3 rgb2ycocg(in vec3 rgb)
{
    float co = rgb.r - rgb.b;
    float t = rgb.b + co / 2.0;
    float cg = rgb.g - t;
    float y = t + cg / 2.0;
    return vec3(y, co, cg);
}


vec3 ycocg2rgb(in vec3 ycocg)
{
    float t = ycocg.r - ycocg.b / 2.0;
    float g = ycocg.b + t;
    float b = t - ycocg.g / 2.0;
    float r = ycocg.g + b;
    return vec3(r, g, b);
}


vec3 clipToAABB(in vec3 cOld, in vec3 cNew, in vec3 centre, in vec3 halfSize)
{
    if (all(lessThanEqual(abs(cOld - centre), halfSize))) {
        return cOld;
    }
    
    vec3 dir = (cNew - cOld);
    vec3 near = centre - sign(dir) * halfSize;
    vec3 tAll = (near - cOld) / dir;
    float t = 1e20;
    for (int i = 0; i < 3; i++) {
        if (tAll[i] >= 0.0 && tAll[i] < t) {
            t = tAll[i];
        }
    }
    
    if (t >= 1e20) {
		return cOld;
    }
    return cOld + dir * t;
}

vec3 ClampColor(vec3 Color, bool CameraMoved) 
{

    vec3 MinColor = vec3(100.0);
	vec3 MaxColor = vec3(-100.0); 
	vec2 TexelSize = 1.0f / textureSize(u_CurrentColorTexture,0);

    for(int x = -1; x <= 1; x++) 
	{
        for(int y = -1; y <= 1; y++) 
		{
            vec3 Sample = texture(u_CurrentColorTexture, v_TexCoords + vec2(x, y) * TexelSize).rgb; 
            Sample = rgb2ycocg(Sample);

            MinColor = min(Sample, MinColor); 
			MaxColor = max(Sample, MaxColor); 
        }
    }

    //return ycocg2rgb(clipToAABB(rgb2ycocg(Color), rgb2ycocg(CenterColor.rgb), mean, stddev));

    float Bias = 0.003f;

    return ycocg2rgb(clipAABB(rgb2ycocg(Color), MinColor, MaxColor));
}

vec4 SampleHistoryBicubic(vec2 uv) 
{
    vec2 resolution = textureSize(u_PreviousColorTexture, 0).xy;
    vec2 invResolution = 1.0f / resolution;
    vec2 position = uv * resolution;

    vec2 center = floor(position - 0.5) + 0.5;
    vec2 f = position - center;
    vec2 f2 = f * f;
    vec2 f3 = f2 * f;

    vec2 w0 = f2 - 0.5 * (f3 + f);
    vec2 w1 = 1.5 * f3 - 2.5 * f2 + 1.0;
    vec2 w3 = 0.5 * (f3 - f2);
    vec2 w2 = 1.0 - w0 - w1 - w3;

    vec2 w12 = w1 + w2;

    vec2 tc0 = (center - 1.0) * invResolution;
    vec2 tc12 = (center + w2 / w12) * invResolution;
    vec2 tc3 = (center + 2.0) * invResolution;

    vec2 uv0 = clamp(vec2(tc12.x, tc0.y), vec2(0.0), vec2(1.0));
    vec2 uv1 = clamp(vec2(tc0.x, tc12.y), vec2(0.0), vec2(1.0));
    vec2 uv2 = clamp(vec2(tc12.x, tc12.y), vec2(0.0), vec2(1.0));
    vec2 uv3 = clamp(vec2(tc3.x, tc12.y), vec2(0.0), vec2(1.0));
    vec2 uv4 = clamp(vec2(tc12.x, tc3.y), vec2(0.0), vec2(1.0));

    float weight0 = w12.x * w0.y;
    float weight1 = w0.x * w12.y;
    float weight2 = w12.x * w12.y;
    float weight3 = w3.x * w12.y;
    float weight4 = w12.x * w3.y;

    vec4 sample0 = texture(u_PreviousColorTexture, uv0) * weight0;
    vec4 sample1 = texture(u_PreviousColorTexture, uv1) * weight1;
    vec4 sample2 = texture(u_PreviousColorTexture, uv2) * weight2;
    vec4 sample3 = texture(u_PreviousColorTexture, uv3) * weight3;
    vec4 sample4 = texture(u_PreviousColorTexture, uv4) * weight4;

    float totalWeight = weight0 + weight1 + 
        weight2 + weight3 + weight4;

    vec4 totalSample = sample0 + sample1 +
        sample2 + sample3 + sample4;

    return totalSample / totalWeight;    

}

vec2 v_sel(vec2 f, float val, float eq, vec2 neq)
{
    vec2 result;
    result.x = (f.x == val) ? eq : neq.x;
    result.y = (f.y == val) ? eq : neq.y;
    return result;
}

vec3 Lanczos(sampler2D img, vec2 uv)
{
    const float M_PI = PI;
	ivec2 size = textureSize(img, 0);
    
    // Lanczos 3
    vec2 UV = uv.xy * size;
    vec2 tc = floor(UV - 0.5) + 0.5;
    vec2 f = UV - tc + 2;

    // compute at f, f-1, f-2, f-3, f-4, and f-5 using trig angle addition
    vec2 fpi = f * M_PI, fpi3 = f * (M_PI / 3.0);
    vec2 sinfpi = sin(fpi), sinfpi3 = sin(fpi3), cosfpi3 = cos(fpi3);
    const float r3 = sqrt(3.0);
    vec2 w0 = v_sel(f, 0, M_PI * M_PI * 1.0 / 3.0, (sinfpi *       sinfpi3) / (f       * f));
    vec2 w1 = v_sel(f, 1, M_PI * M_PI * 2.0 / 3.0, (-sinfpi * (sinfpi3 - r3 * cosfpi3)) / ((f - 1.0)*(f - 1.0)));
    vec2 w2 = v_sel(f, 2, M_PI * M_PI * 2.0 / 3.0, (sinfpi * (-sinfpi3 - r3 * cosfpi3)) / ((f - 2.0)*(f - 2.0)));
    vec2 w3 = v_sel(f, 3, M_PI * M_PI * 2.0 / 3.0, (-sinfpi * (-2.0*sinfpi3)) / ((f - 3.0)*(f - 3.0)));
    vec2 w4 = v_sel(f, 4, M_PI * M_PI * 2.0 / 3.0, (sinfpi * (-sinfpi3 + r3 * cosfpi3)) / ((f - 4.0)*(f - 4.0)));
    vec2 w5 = v_sel(f, 5, M_PI * M_PI * 2.0 / 3.0, (-sinfpi * (sinfpi3 + r3 * cosfpi3)) / ((f - 5.0)*(f - 5.0)));

    // use bilinear texture weights to merge center two samples in each dimension
    vec2 Weight[5];
    Weight[0] = w0;
    Weight[1] = w1;
    Weight[2] = w2 + w3;
    Weight[3] = w4;
    Weight[4] = w5;

    vec2 invTextureSize = 1.0 / vec2(size);

    vec2 Sample[5];
    Sample[0] = invTextureSize * (tc - 2);
    Sample[1] = invTextureSize * (tc - 1);
    Sample[2] = invTextureSize * (tc + w3 / Weight[2]);
    Sample[3] = invTextureSize * (tc + 2);
    Sample[4] = invTextureSize * (tc + 3);

    vec4 o_rgba = vec4(0);

    // 5x5 footprint with corners dropped to give 13 texture taps
    o_rgba += vec4(textureLod(img, vec2(Sample[0].x, Sample[2].y), 0).rgb, 1.0) * Weight[0].x * Weight[2].y;
    o_rgba += vec4(textureLod(img, vec2(Sample[1].x, Sample[1].y), 0).rgb, 1.0) * Weight[1].x * Weight[1].y;
    o_rgba += vec4(textureLod(img, vec2(Sample[1].x, Sample[2].y), 0).rgb, 1.0) * Weight[1].x * Weight[2].y;
    o_rgba += vec4(textureLod(img, vec2(Sample[1].x, Sample[3].y), 0).rgb, 1.0) * Weight[1].x * Weight[3].y;
    o_rgba += vec4(textureLod(img, vec2(Sample[2].x, Sample[0].y), 0).rgb, 1.0) * Weight[2].x * Weight[0].y;
    o_rgba += vec4(textureLod(img, vec2(Sample[2].x, Sample[1].y), 0).rgb, 1.0) * Weight[2].x * Weight[1].y;
    o_rgba += vec4(textureLod(img, vec2(Sample[2].x, Sample[2].y), 0).rgb, 1.0) * Weight[2].x * Weight[2].y;
    o_rgba += vec4(textureLod(img, vec2(Sample[2].x, Sample[3].y), 0).rgb, 1.0) * Weight[2].x * Weight[3].y;
    o_rgba += vec4(textureLod(img, vec2(Sample[2].x, Sample[4].y), 0).rgb, 1.0) * Weight[2].x * Weight[4].y;
    o_rgba += vec4(textureLod(img, vec2(Sample[3].x, Sample[1].y), 0).rgb, 1.0) * Weight[3].x * Weight[1].y;
    o_rgba += vec4(textureLod(img, vec2(Sample[3].x, Sample[2].y), 0).rgb, 1.0) * Weight[3].x * Weight[2].y;
    o_rgba += vec4(textureLod(img, vec2(Sample[3].x, Sample[3].y), 0).rgb, 1.0) * Weight[3].x * Weight[3].y;
    o_rgba += vec4(textureLod(img, vec2(Sample[4].x, Sample[2].y), 0).rgb, 1.0) * Weight[4].x * Weight[2].y;

    return o_rgba.rgb / o_rgba.w;
}


void main()
{
    vec2 CSampleCoord = v_TexCoords;
    CSampleCoord -= u_CurrentJitter / (textureSize(u_CurrentColorTexture, 0) * 2.);
	vec3 CurrentColor = Lanczos(u_CurrentColorTexture, CSampleCoord).rgb;

	if (!u_Enabled) {
		o_Color = CurrentColor;
		return;
	}

	float CurrentDepth = texture(u_DepthTexture, v_TexCoords).x;

	vec3 WorldPosition = WorldPosFromDepth(CurrentDepth, v_TexCoords).xyz;
	vec3 PreviousCoord = Reprojection(WorldPosition.xyz); 
	float bias = (u_InverseView != u_InversePrevView) ? 0.01f : 0.0f;

	if (PreviousCoord.x > bias && PreviousCoord.x < 1.0f-bias &&
		PreviousCoord.y > bias && PreviousCoord.y < 1.0f-bias && 
		v_TexCoords.x > bias && v_TexCoords.x < 1.0f-bias &&
		v_TexCoords.y > bias && v_TexCoords.y < 1.0f-bias)
	{
        // TODO : Replace this with a motion vector check (this is actually braindead)
        bool CameraMoved = distance(u_InversePrevView[3].xyz, u_InverseView[3].xyz) > 0.0005f;

		// Depth rejection
		float PreviousDepth = texture(u_PreviousDepthTexture, PreviousCoord.xy).x;
		float LinearPrevDepth = linearizeDepth(PreviousDepth);
		float LinearExpectedDepth = linearizeDepth(PreviousCoord.z);

		vec3 PrevColor = SampleHistoryBicubic(PreviousCoord.xy).rgb;
		
	    PrevColor = ClampColor(PrevColor, CameraMoved);

		vec2 Dimensions = textureSize(u_CurrentColorTexture, 0).xy;
		vec2 velocity = (v_TexCoords - PreviousCoord.xy) * Dimensions;
		float BlendFactor = exp(-length(velocity)) * 0.9f + 0.6f;
		BlendFactor = clamp(BlendFactor, 0.0f, 0.95f);

		const float DepthRejectionStrength = 0.0f; // 4.0f
		BlendFactor *= pow(exp(-abs(LinearExpectedDepth-LinearPrevDepth)), 32.0f * DepthRejectionStrength);
		o_Color = mix(CurrentColor.xyz, PrevColor.xyz, clamp(BlendFactor, 0.01f, 0.95f));
	}

	else 
	{
		o_Color = CurrentColor;
	}
}

