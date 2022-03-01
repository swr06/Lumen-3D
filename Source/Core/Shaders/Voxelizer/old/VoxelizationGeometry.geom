#version 450 core

layout (triangles) in;

layout (triangle_strip, max_vertices = 3) out;

in vec3 v_WorldPosition[];
in vec3 v_Normal[];

out vec3 g_WorldPosition;
out vec3 g_Normal;
out vec3 g_VolumePosition;

uniform vec3 u_VoxelGridCenter;

uniform int u_VolumeSize;
uniform vec3 u_CoverageSize;

int GetLargestAxis(in vec3 x) {
	int Dominantaxis = x[1] > x[0] ? 1 : 0;
	Dominantaxis = x[2] > x[Dominantaxis] ? 2 : Dominantaxis;
	return Dominantaxis;
}

vec3 ToVoxelSpace(vec3 Position) {
	
	vec3 HalfExtent = u_CoverageSize / 2.0f;
	return (Position - u_VoxelGridCenter) / HalfExtent;
}

void main() {

	vec3 FaceNormal = abs(v_Normal[0] + v_Normal[1] + v_Normal[2]);

	// Figure out index of maximum axis 
	
	int LargestAxis = GetLargestAxis(FaceNormal);


	for (int Vertex = 0 ; Vertex < 3 ; Vertex++) {

		vec3 WorldPosition = v_WorldPosition[Vertex];

		g_WorldPosition = WorldPosition;
		g_Normal = v_Normal[Vertex];

		vec3 VoxelVolumePosition = ToVoxelSpace(WorldPosition);

		// Orthagonal projection to dominant axis 
		
		gl_Position.xyz = VoxelVolumePosition.xyz;

		if (LargestAxis == 0) {
			
			gl_Position.xyz = VoxelVolumePosition.zyx;

		}

		else if (LargestAxis == 1) {
			
			gl_Position.xyz = VoxelVolumePosition.xzy;

		}

		g_VolumePosition = gl_Position.xyz;

		// Transform to clip space 

		gl_Position.xyz /= vec3(u_VolumeSize);

		gl_Position.z = 0.0f;
		gl_Position.w = 1.0f;

		EmitVertex();
	}

	EndPrimitive();
}