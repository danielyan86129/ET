// A simple 2-simplicial complex mesh
// 
// Copyright (C) 2018 Yajie Yan <danielyan86129@hotmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.

#include <string>

using std::string;

#include <plyall.h>
#include "mymesh.h"

struct internal{
	struct Vertex { float x, y, z; };
	struct Edge { int v1, v2; };
	struct Face {
		unsigned char nvts;
		int verts[ 3 ];
		float s;
	};
};
MyMesh * MyMesh::read( const char * _fname )
{
	string fname( _fname );
	if ( fname.substr( fname.size() - string( ".off" ).size() ) == ".off" )
	{
		auto mesh = TriMesh::read( fname.c_str() );
		auto newmesh = new MyMesh;
		newmesh->vertices = mesh->vertices;
		newmesh->faces = mesh->faces;
		delete mesh;
		return newmesh;
	}
	else if ( fname.substr( fname.size() - string( ".ply" ).size() ) == ".ply" )
	{
		return read_ply_helper( fname.c_str() );
	}
	else if ( fname.substr( fname.size() - string( ".ma" ).size() ) == ".ma" )
	{
		vector<float> r;
		return readQMATFile( fname, r );
	}
	return nullptr;
}

MyMesh * MyMesh::read( const::std::string & filename )
{
	return read( filename.c_str() );
}

void MyMesh::write( const char * _fname, const MyMesh * _m )
{
	string fname( _fname );
	if ( fname.find( ".off" ) != string::npos )
		( (TriMesh*)_m )->write( _fname );
	else if ( fname.find( ".ply" ) != string::npos )
		write_ply_helper( fname.c_str(), _m );
}

void MyMesh::write( const std::string & _fname, const MyMesh * _m )
{
	write( _fname.c_str(), _m );
}

MyMesh * MyMesh::read_ply_helper( const char* _fname )
{
	map<string, PlyProperty> v_props_map;
	PlyProperty v_props[] = {
		{ "x", Float32, Float32, offsetof( internal::Vertex, x ), 0, 0, 0, 0 },
		{ "y", Float32, Float32, offsetof( internal::Vertex, y ), 0, 0, 0, 0 },
		{ "z", Float32, Float32, offsetof( internal::Vertex, z ), 0, 0, 0, 0 }
	};
	v_props_map[ "x" ] = v_props[0];
	v_props_map[ "y" ] = v_props[1];
	v_props_map[ "z" ] = v_props[2];
	map<string, PlyProperty> e_props_map;
	PlyProperty e_props[] = {
		{ "vertex1", Int32, Int32, offsetof( internal::Edge, v1 ), PLY_SCALAR, 0, 0, 0 },
		{ "vertex2", Int32, Int32, offsetof( internal::Edge, v2 ), PLY_SCALAR, 0, 0, 0 }
	};
	e_props_map[ "vertex1" ] = e_props[0];
	e_props_map[ "vertex2" ] = e_props[1];
	std::map<std::string, PlyProperty> f_props_map;
	PlyProperty f_props[] = {
		{"vertex_indices", Int32, Int32, offsetof( internal::Face, verts ),
			PLY_LIST, Uint8, Uint8, offsetof( internal::Face,nvts ) },
		{"msure", Float32, Float32, offsetof( internal::Face, s ),
			PLY_SCALAR, 0,0,0}
	};
	f_props_map[ "vertex_indices" ] = f_props[0];
	f_props_map[ "msure" ] = f_props[ 1 ];
	vector<internal::Vertex> vts;
	vector<internal::Edge> edges;
	vector<internal::Face> faces;

	ply::PLYreader ply_reader;
	if ( ply_reader.read( _fname, v_props_map, e_props_map, f_props_map, vts, edges, faces ) != ply::SUCCESS )
	{
		cout << "failed to read mesh file " << _fname << endl;
		return nullptr;
	}

	auto mesh = new MyMesh();
	for ( auto i = 0; i < vts.size(); ++i )
	{
		const auto& p = vts[ i ];
		mesh->vertices.push_back( TriPoint( p.x, p.y, p.z ) );
	}
	vts.clear(); vts.shrink_to_fit();
	for ( auto i = 0; i < edges.size(); ++i )
	{
		const auto& e = edges[ i ];
		mesh->lines.push_back( util::makeEdge( e.v1, e.v2 ) );
	}
	edges.clear(); edges.shrink_to_fit();
	for ( auto i = 0; i < faces.size(); ++i )
	{
		const auto& f = faces[ i ];
		mesh->faces.push_back( util::makeFace(TriFace( f.verts )) );
		mesh->msure_face.push_back( f.s );
	}
	faces.clear(); faces.shrink_to_fit();

	return mesh;
}

bool MyMesh::write_ply_helper( const char * _fname, const MyMesh * _m )
{
	struct Vertex
	{
		float x; float y; float z;
		unsigned char r, g, b;
		float s; // measure
	};
	struct Edge
	{
		int v1; int v2;
		unsigned char r, g, b;
		float s; // measure
	};
	struct Face
	{
		unsigned char nvts;
		int verts[ 3 ];
		unsigned char r, g, b;
		float s; // measure
	};
	std::map<std::string, PlyProperty> vert_props;
	vert_props[ "x" ] = { "x", Float32, Float32, offsetof( Vertex, x ), PLY_SCALAR, 0, 0, 0 };
	vert_props[ "y" ] = { "y", Float32, Float32, offsetof( Vertex, y ), PLY_SCALAR, 0, 0, 0 };
	vert_props[ "z" ] = { "z", Float32, Float32, offsetof( Vertex, z ), PLY_SCALAR, 0, 0, 0 };
	std::map<std::string, PlyProperty> edge_props;
	edge_props[ "vertex1" ] = { "vertex1", Int32, Int32, offsetof( Edge, v1 ), PLY_SCALAR, 0, 0, 0 };
	edge_props[ "vertex2" ] = { "vertex2", Int32, Int32, offsetof( Edge, v2 ), PLY_SCALAR, 0, 0, 0 };
	std::map<std::string, PlyProperty> face_props;
	face_props[ "vertex_indices" ] = {
		"vertex_indices", Int32, Int32, offsetof( Face, verts ),
		PLY_LIST, Uint8, Uint8, offsetof( Face,nvts ) };
	vector<Vertex> output_vts( _m->vertices.size() );
	vector<Edge> output_edges( _m->lines.size() );
	vector<Face> output_faces( _m->faces.size() );
	// fill in the data to output lists (face might contain color)
	for ( auto i = 0; i < output_vts.size(); ++i )
	{
		auto& o_p = output_vts[ i ];
		o_p.x = _m->vertices[ i ][ 0 ];
		o_p.y = _m->vertices[ i ][ 1 ];
		o_p.z = _m->vertices[ i ][ 2 ];
	}
	for ( auto i = 0; i < output_edges.size(); ++i )
	{
		auto& o_e = output_edges[ i ];
		o_e.v1 = _m->lines[ i ][ 0 ];
		o_e.v2 = _m->lines[ i ][ 1 ];
	}
	for ( auto i = 0; i < output_faces.size(); ++i )
	{
		auto& o_f = output_faces[ i ];
		o_f.nvts = 3;
		o_f.verts[ 0 ] = _m->faces[ i ][ 0 ];
		o_f.verts[ 1 ] = _m->faces[ i ][ 1 ];
		o_f.verts[ 2 ] = _m->faces[ i ][ 2 ];
	}

	std::cout << "writing to ply file... " << std::endl;
	std::cout << "# vts/edges/faces: "
		<< output_vts.size() << "/"
		<< output_edges.size() << "/"
		<< output_faces.size() << std::endl;

	ply::PLYwriter ply_writer;
	return ply_writer.write( _fname, true, true, true,
		vert_props, edge_props, face_props,
		output_vts, output_edges, output_faces ) == ply::SUCCESS;
}

MyMesh * MyMesh::readQMATFile( const std::string _filename, vector<float>& _radii )
{
	ifstream infile( _filename );
	if ( !infile.is_open() )
	{
		return nullptr;
	}

	int n_vts, n_edges, n_faces, dummy_int;
	float r, dummy_float;
	std::string data_type_s;
	infile >> n_vts >> n_edges >> n_faces /*>> dummy_int*/;
	vector<TriEdge> edges;
	set<TriEdge> edges_used;
	vector<TriEdge> es_of_f;

	// create mesh
	auto mesh = new MyMesh();

	// read vts
	for ( int i = 0; i < n_vts; ++i )
	{
		point v;
		infile >> data_type_s;
		if ( data_type_s != "v" )
			goto READ_QMAT_FAIL;
		infile >> v[ 0 ] >> v[ 1 ] >> v[ 2 ] >> r/*>>dummy_float*/;
		mesh->vertices.push_back( v );
	}
	// read edges
	for ( int i = 0; i < n_edges; ++i )
	{
		TriEdge e;
		infile >> data_type_s;
		if ( data_type_s != "e" )
			goto READ_QMAT_FAIL;
		infile >> e[ 0 ] >> e[ 1 ];
		edges.push_back( util::makeEdge( e[ 0 ], e[ 1 ] ) );
	}
	// read faces
	for ( int i = 0; i < n_faces; ++i )
	{
		TriFace f;
		infile >> data_type_s;
		if ( data_type_s != "f" )
			goto READ_QMAT_FAIL;
		infile >> f[ 0 ] >> f[ 1 ] >> f[ 2 ];
		mesh->faces.push_back( util::makeFace( f ) );
	}
	infile.close();

	// keep isolated edges, i.e. not used by any faces
	for ( int i = 0; i < n_faces; ++i )
	{
		const auto& f = mesh->faces[ i ];
		es_of_f = util::edgesFromFace( f );
		for ( auto it = es_of_f.begin(); it != es_of_f.end(); ++it )
		{
			edges_used.insert( *it );
		}
	}
	for ( int i = 0; i < n_edges; ++i )
	{
		const auto& e = edges[ i ];
		if ( edges_used.find( e ) == edges_used.end() )
		{
			mesh->lines.push_back( e );
		}
	}
	edges_used.clear();
	edges.clear();

	return mesh;

READ_QMAT_FAIL:
	if ( mesh )
		delete mesh;
	return nullptr;
}