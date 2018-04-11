#pragma once
#include "commonTypes.h"
#include "commonUtils.h"

#include <vector>
using std::vector;

class MyMesh : public TriMesh
{
public:
	MyMesh() 
	{
		TriMesh();
	}

	static MyMesh *read( const char *filename );
	static MyMesh *read( const ::std::string &filename );
	static void write( const char *_fname, const MyMesh* _m );
	static void write( const std::string& _fname, const MyMesh* _m );
protected:
	static MyMesh * read_ply_helper( const char * _fname );
	static bool write_ply_helper( const char* _fname, const MyMesh* _m );
public:
	// 1-d line geometry
	vector<TriEdge> lines;
	// measure assoc. on element
	vector<float> msure_face;
};