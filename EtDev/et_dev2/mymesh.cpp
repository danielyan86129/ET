#include <string>

using std::string;

#include "mymesh.h"
#include "plyall.h"

struct internal{
	struct Vertex { float x, y, z; };
	struct Edge { int v1, v2; };
	struct Face {
		unsigned char nvts;
		int verts[ 3 ];
	};
};
MyMesh * MyMesh::read( const char * _fname )
{
	string fname( _fname );
	if ( fname.find( ".off" ) != string::npos )
	{
		auto mesh = TriMesh::read( fname.c_str() );
		auto newmesh = new MyMesh;
		newmesh->vertices = mesh->vertices;
		newmesh->faces = mesh->faces;
		delete mesh;
		return newmesh;
	}
	else if ( fname.find( ".ply" ) != string::npos )
	{
		return read_ply_helper( fname.c_str() );
	}
	return nullptr;
}

MyMesh * MyMesh::read( const::std::string & filename )
{
	return read( filename.c_str() );
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
			PLY_LIST, Uint8, Uint8, offsetof( internal::Face,nvts ) }
	};
	f_props_map[ "vertex_indices" ] = f_props[0];
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
		mesh->lines.push_back( TriEdge( e.v1, e.v2 ) );
	}
	edges.clear(); edges.shrink_to_fit();
	for ( auto i = 0; i < faces.size(); ++i )
	{
		const auto& f = faces[ i ];
		mesh->faces.push_back( TriFace( f.verts ) );
	}
	return mesh;
}
