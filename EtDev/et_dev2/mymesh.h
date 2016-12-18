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
protected:
	static MyMesh * read_ply_helper( const char * _fname );
public:
	// 1-d line geometry
	vector<TriEdge> lines;
};