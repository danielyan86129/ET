#ifndef COMMON_TYPES_H
#define COMMON_TYPES_H

// need CGAL's NN search cuz it supports associating point with info like index, etc.
#include <CGAL/Simple_cartesian.h> 
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/K_neighbor_search.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/property_map.h>
#include <boost/iterator/zip_iterator.hpp>
#include <memory>
#include <vector>
#include <cmath>

using std::shared_ptr;
using std::vector;

#include <TriMesh.h>
#include <TriMesh_algo.h>

using namespace trimesh;

class Drawable;
class TrackBall;
class MyGraph;
class SteinerGraph;
struct TopoGraph;
class HybridSkeleton;

typedef trimesh::ivec2 TriEdge;
typedef trimesh::point TriPoint;
typedef TriMesh::Face TriFace;
typedef trimesh::vec3 TriColor;

// Functor to define order of edges.
struct EdgeCompare
{
	bool operator () (const TriEdge& _e1, const TriEdge& _e2) const
	{
		if (_e1[0] == _e2[0])
			return _e1[1] < _e2[1];
		return _e1[0] < _e2[0];
	}
};

struct FaceCompare
{
	bool operator () (const TriFace& _f1, const TriFace& _f2) const
	{
		return _f1[0] == _f2[0] ? (_f1[1] == _f2[1] ? _f1[2] < _f2[2] : _f1[1] < _f2[1]) : _f1[0] < _f2[0];
	}
};

// typedefs for CGAL's NN search
typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3 P3;

typedef boost::tuple<P3,unsigned>	Point_and_uint;
typedef CGAL::Search_traits_3<K>	Traits_base;

typedef CGAL::Search_traits_adapter<Point_and_uint,
	CGAL::Nth_of_tuple_property_map<0, Point_and_uint>,
	Traits_base>    Traits;

// typedefs for CGAL's K-NN search
typedef CGAL::Orthogonal_k_neighbor_search<Traits>	K_neighbor_search;
typedef K_neighbor_search::Tree                     KNNTree;
typedef K_neighbor_search::Distance                 Distance;

// typedefs for CGAL's range search (e.g. using a fuzzy sphere)
typedef CGAL::Kd_tree<Traits> KdTree;
typedef CGAL::Fuzzy_sphere<Traits> FuzySphere3;

#endif