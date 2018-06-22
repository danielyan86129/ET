// A simple 2-simplicial complex mesh
// 
// Copyright (C) 2018 Yajie Yan <danielyan86129@hotmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.

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
	static MyMesh * readQMATFile( std::string _filename, vector<float>& _radii );
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