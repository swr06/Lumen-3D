#version 330 core

layout (location = 0) in vec2 a_Position;
layout (location = 1) in vec2 a_TexCoords;

// yes
out vec2 v_TexCoords;
out vec2 v_TexCoord;
out vec2 v_Texcoord;
out vec2 v_Texcoords;
out vec2 v_texcoords;
out vec2 texcoord;
out vec2 Texcoord;
out vec2 Texcoords;

void main()
{
	gl_Position = vec4(a_Position, 0.0f, 1.0f);
	v_TexCoords = a_TexCoords;
	v_TexCoord = a_TexCoords;
	v_Texcoord = a_TexCoords;
	v_Texcoords = a_TexCoords;
	v_texcoords = a_TexCoords;
	texcoord = a_TexCoords;
	Texcoord = a_TexCoords;
	Texcoords = a_TexCoords;
}