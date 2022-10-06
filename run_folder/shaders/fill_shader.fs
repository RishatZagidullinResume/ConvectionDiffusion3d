#version 330 core

out vec4 FragColor;

in float ourColor;

vec3 output;

void main()
{
    if (ourColor < 0.5)
    {
        output = vec3(0.0, 2*ourColor, 2*(0.5-ourColor));
    }
    else
    {
        output = vec3(2*(ourColor-0.5), 2*(1.0-ourColor), 0.0);
    }

    FragColor = vec4(output, 1.0);
    
}
