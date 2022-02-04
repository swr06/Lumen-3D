#version 330 core

layout (location = 0) out vec3 o_Color;

in vec2 v_TexCoords;

void main() {
	o_Color = vec3(1.0f, 0.0f, 0.0f);
}