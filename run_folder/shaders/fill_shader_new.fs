#version 330 core

out vec4 FragColor;

in float ourColor;

vec3 output;

void main()
{
    if (ourColor < 0.0)
    {
        output = vec3(0.0, 1.+ourColor, -ourColor);
    }
    else
    {
        output = vec3(ourColor, 1.0-ourColor, 0.0);
    }

    FragColor = vec4(output, 1.0);
    
}
