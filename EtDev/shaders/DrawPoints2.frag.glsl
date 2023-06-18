#version 330

in vec3 vertColor;
in float vertSaliency;

out vec4 fragColor;

void main(void)
{
	//fragColor = vec4(0.6, 0.3, 0.3, 1.0);
	//vec3 temp_c = vec3(0.0f, 0.0f, 1.0f);
	//fragColor = vec4(temp_c, 0.5f);
	fragColor = vec4(vertColor, vertSaliency);

}