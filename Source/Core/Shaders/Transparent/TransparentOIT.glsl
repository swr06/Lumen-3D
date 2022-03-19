#version 330 core

layout (location = 0) out vec4 o_Color;

// Vertex shader inputs ->
in vec2 v_TexCoords;
in vec3 v_FragPosition;
in vec3 v_Normal;
in mat3 v_TBNMatrix;

uniform sampler2D u_AlbedoMap; // <- Glass albedo 
uniform sampler2D u_NormalMap; // <- Glass normals

uniform sampler2D u_LowResDepth; // <- 1/4 Resolution depth
uniform sampler2D u_ColorBuffer; // <- Combined lighting of opaque objects 

uniform bool u_UsesNormalMap;
uniform bool u_UsesAlbedoTexture;
uniform vec4 u_ModelColor;

// Camera stuff 
uniform float u_zNear;
uniform float u_zFar;
uniform mat4 u_ViewProjection;
uniform mat4 u_Projection;
uniform mat4 u_View;
uniform mat4 u_InverseProjection;
uniform mat4 u_InverseView;

uniform float u_Time;

struct GBufferData {
   vec3 Position;
   vec2 UV;
   bool ValidMask;
};

vec3 WorldPosFromDepth(float depth, vec2 txc)
{
    float z = depth * 2.0 - 1.0;
    vec4 ClipSpacePosition = vec4(txc * 2.0 - 1.0, z, 1.0);
    vec4 ViewSpacePosition = u_InverseProjection * ClipSpacePosition;
    ViewSpacePosition /= ViewSpacePosition.w;
    vec4 WorldPos = u_InverseView * ViewSpacePosition;
    return WorldPos.xyz;
}

vec3 ProjectToScreenSpace(vec3 WorldPos) 
{
	vec4 ProjectedPosition = u_ViewProjection * vec4(WorldPos, 1.0f);
	ProjectedPosition.xyz /= ProjectedPosition.w;
	ProjectedPosition.xyz = ProjectedPosition.xyz * 0.5f + 0.5f;
	return ProjectedPosition.xyz;
}

vec3 ProjectToClipSpace(vec3 WorldPos) 
{
	vec4 ProjectedPosition = u_ViewProjection * vec4(WorldPos, 1.0f);
	ProjectedPosition.xyz /= ProjectedPosition.w;
	return ProjectedPosition.xyz;
}

float LinearizeDepth(float depth)
{
	return (2.0 * u_zNear) / (u_zFar + u_zNear - depth * (u_zFar - u_zNear));
}

// Pseudo rng 
float HASH2SEED = 0.0f;
vec2 hash2() 
{
	return fract(sin(vec2(HASH2SEED += 0.1, HASH2SEED += 0.1)) * vec2(43758.5453123, 22578.1459123));
}

void main()
{
	HASH2SEED = (v_TexCoords.x * v_TexCoords.y) * 64.0;
	HASH2SEED += fract(u_Time) * 64.0f;

	vec4 AlbedoFetch = (u_UsesAlbedoTexture ? texture(u_AlbedoMap, v_TexCoords).xyzw : u_ModelColor);
	vec3 NormalFetch = u_UsesNormalMap ? normalize(v_TBNMatrix * (texture(u_NormalMap, v_TexCoords).xyz * 2.0f - 1.0f)) : v_Normal;

	float Hash = hash2().x;
	float Alpha = AlbedoFetch.w;



}

// Checks if a position lies in screen space (with some bias)
bool SSRayValid(vec2 x) {
	float bias = 0.0001f;
	if (x.x > bias && x.x < 1.0f - bias && x.y > bias && x.y < 1.0f - bias) {
		return true;
	}

	return false;
}

// Fake refraction if trace fails
vec2 FakeRefract(vec2 SS, vec3 Position, vec3 BackPosition, vec3 Normal, vec3 Incident, vec3 Direction) {

	bool TotalInternalRefraction = (dot(Direction, Direction) < 0.01f);
	float RefractionAmount = distance(Position, BackPosition);
    vec3 ApproximateHitPosition = Position + Direction * RefractionAmount;
    vec3 HitCoordinate = ProjectToScreenSpace(ApproximateHitPosition);
    
	if (!SSRayValid(HitCoordinate.xy)) {
		return SS;
	}

	return HitCoordinate.xy;
}

// Screenspace raytracing 
GBufferData ScreenspaceRaytrace(vec3 Origin, vec3 Direction, float Hash)
{
    const float Distance = 196.0f;
	const float ThresholdMultiplier = 0.001f;
	const int Steps = 12;
	const int BinarySteps = 6;
    const float ExpStep = 1.0225f; 

    GBufferData ReturnValue;
    ReturnValue.Position = Origin;
    ReturnValue.ValidMask = false;
	ReturnValue.UV = vec2(-1.0f);

    float StepSize = float(Distance) / float(Steps);
	vec3 RayPosition = Origin + Direction * Hash; 
	vec2 FinalUV = vec2(-1.0f);


	for(int CurrentStep = 0; CurrentStep < Steps; CurrentStep++) 
	{
        float ToleranceStep = mix(2.2f, 6.0f, pow(float(CurrentStep) / float(Steps), 1.5f));

		float Threshold = StepSize * ThresholdMultiplier * ToleranceStep;
		
		vec3 ProjectedRayScreenspace = ProjectToClipSpace(RayPosition); 
		
		if(abs(ProjectedRayScreenspace.x) > 1.0f || abs(ProjectedRayScreenspace.y) > 1.0f || abs(ProjectedRayScreenspace.z) > 1.0f) 
		{
			return ReturnValue;
		}
		
		ProjectedRayScreenspace.xyz = ProjectedRayScreenspace.xyz * 0.5f + 0.5f; 

		if (!SSRayValid(ProjectedRayScreenspace.xy))
		{
			return ReturnValue;
		}
		
		float DepthAt = texture(u_LowResDepth, ProjectedRayScreenspace.xy).x; 
		float CurrentRayDepth = LinearizeDepth(ProjectedRayScreenspace.z); 
		float Error = abs(LinearizeDepth(DepthAt) - CurrentRayDepth);
		
        // Intersected!
		if (Error < Threshold && ProjectedRayScreenspace.z > DepthAt) 
		{
			// Binary search for best intersection point along ray step 

            bool DoBinaryRefinement = true;

            vec3 FinalProjected = vec3(0.0f);
            float FinalDepth = 0.0f;

            if (DoBinaryRefinement) {
			    vec3 BinaryStepVector = (Direction * StepSize) / 2.0f;
                RayPosition -= (Direction * StepSize) / 2.0f;
			    
                for (int BinaryStep = 0 ; BinaryStep < BinarySteps ; BinaryStep++) {
			    		
			    	BinaryStepVector /= 2.0f;
			    	vec3 Projected = ProjectToClipSpace(RayPosition); 
			    	Projected = Projected * 0.5f + 0.5f;
                    FinalProjected = Projected;
                    float Fetch = texture(u_LowResDepth, Projected.xy).x;
                    FinalDepth = Fetch;
			    	float BinaryDepthAt = LinearizeDepth(Fetch); 
			    	float BinaryRayDepth = LinearizeDepth(Projected.z); 

			    	if (BinaryDepthAt < BinaryRayDepth) {
			    		RayPosition -= BinaryStepVector;

			    	}

			    	else {
			    		RayPosition += BinaryStepVector;

			    	}

			    }
            }

            else {
                
                FinalProjected = ProjectToScreenSpace(RayPosition);
                FinalDepth = texture(u_LowResDepth, FinalProjected.xy).x;
            }


            ReturnValue.Position = WorldPosFromDepth(FinalDepth, FinalProjected.xy);
            ReturnValue.ValidMask = true;
			ReturnValue.UV = FinalProjected.xy;

            return ReturnValue;
		}

        // Step 
		RayPosition += StepSize * Direction; 

        StepSize *= ExpStep;
	}

	return ReturnValue;

}
