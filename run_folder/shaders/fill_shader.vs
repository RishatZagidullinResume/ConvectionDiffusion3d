#version 330 core
layout (location = 0) in vec3 aPos;
layout (location = 1) in float aColor;

uniform mat4 projection;
uniform mat4 view;

out float ourColor;

void main()
{
    gl_Position = projection*view*vec4(aPos, 1.0);
    ourColor = aColor;
}
