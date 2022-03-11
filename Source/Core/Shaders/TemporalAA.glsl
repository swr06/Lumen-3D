#version 330 core

layout (location = 0) out vec3 o_Color;

in vec2 v_TexCoords;

uniform sampler2D u_CurrentColorTexture;
uniform sampler2D u_DepthTexture;
uniform sampler2D u_PreviousColorTexture;
uniform sampler2D u_PreviousDepthTexture;

uniform mat4 u_InverseView;
uniform mat4 u_InverseProjection;
uniform mat4 u_PrevProjection;
uniform mat4 u_PrevView;
uniform mat4 u_InversePrevProjection;
uniform mat4 u_InversePrevView;

uniform bool u_Enabled;

uniform float u_zNear;
uniform float u_zFar;

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
            MinColor = min(Sample, MinColor); 
			MaxColor = max(Sample, MaxColor); 
        }
    }

    float Bias = CameraMoved ? 0.01f : 0.025f;

    return clipAABB(Color, MinColor - Bias, MaxColor + Bias);
}

vec4 SampleHistory(vec2 uv) 
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

void main()
{
	vec3 CurrentColor = texture(u_CurrentColorTexture, v_TexCoords).rgb;

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
        bool CameraMoved = distance(u_InversePrevView[3].xyz, u_InverseView[3].xyz) > 0.0005f;

		// Depth rejection
		float PreviousDepth = texture(u_PreviousDepthTexture, PreviousCoord.xy).x;
		float LinearPrevDepth = linearizeDepth(PreviousDepth);
		float LinearExpectedDepth = linearizeDepth(PreviousCoord.z);

		vec3 PrevColor = SampleHistory(PreviousCoord.xy).rgb;
		
		if (CameraMoved) {
			PrevColor = ClampColor(PrevColor, CameraMoved);
		}

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

