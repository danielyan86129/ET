#ifndef GL_UTILS_H
#define GL_UTILS_H

/*------------------------- math utilities -------------------------------*/
#include <oglplus/all.hpp>

using namespace oglplus;

namespace
{
	

	inline Vec3f crossProduct(const Vec3f& _v1, const Vec3f& _v2)
	{
		return Vec3f(
			_v1.y()*_v2.z()-_v1.z()*_v2.y(), 
			_v1.z()*_v2.x()-_v1.x()*_v2.z(), 
			_v1.x()*_v2.y()-_v1.y()*_v2.x()
			);
	}

};

/*------------------------- gl utilities -------------------------------*/
#include <fstream>

namespace
{
	/// load and return entire source code from given shader file
	oglplus::String loadShaderSource(const char * glsl_file)
	{
		std::ifstream file(glsl_file, std::ios::in | std::ios::binary);

		if (!file.is_open())
		{
			std::cout << "Error: can't find file " << glsl_file 
				<< ". Program exiting... "<< std::endl;
			exit(1);
		}

		std::string contents;
		file.seekg(0, std::ios::end);
		contents.resize(file.tellg());
		file.seekg(0, std::ios::beg);
		file.read(&contents[0], contents.size());
		file.close();

		return contents;
	}

};



#endif