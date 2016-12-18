#ifndef GEOMETRY_REPRESENTATION_H
#define GEOMETRY_REPRESENTATION_H

#include <TriMesh.h>
#include <TriMesh_algo.h>

#include <vector>
#include <set>
#include <map>
#include <cassert>
#include <cstdint>

#include <boost/pool/pool_alloc.hpp>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_walk_along_line_point_location.h>

#include <Eigen/Dense> // use Eigen's linear eq. solver
#include <Eigen/SVD>

#include "commonUtils.h"
#include "allTypes.h"
#include "steinerSubdivision.h"
#include "origGraph.h"
#include "graphDiffusion.h"

using namespace std;
using namespace trimesh;

// give same coordinates of vertices of those faces whose spheres have no intersection
//void findFacesWithNonIsectSpheres();
// preprocess the mesh: remove degenerate triangles
void preprocess(shared_ptr<MyMesh> _m, vector<float>& _radii, float _eps);
// for each face in the given MA mesh, find in the given 3d model the first 2 closest points 
void find2Closest(shared_ptr<MyMesh> _mMA, shared_ptr<MyMesh> _m3d, vector<TriPoint>& _closest);

typedef trimesh::Vec<2, uint8_t> Tri2Byte;

extern TriColor topo_colors[5];
enum TopoType {
	ISOLATE = 0, // isolated 0-d
	MANIFOLD_1D, // in middle of 1d manifold
	JUNCTION_1D, // at 1d junction
	BOUNDARY_2D, // on 2d boundary
	JUNCTION_2D, // at 2d junction
	MANIFOLD_2D, // on 2d manifold
	NONMANIFOLD_2D, // on 2d singular seam
	MIX_2D_1D // mixed of 2d and 1d
};
extern const char* topo_type_str[7];

class SteinerGraph;
struct TopoGraph
{
	// <topo edge rep, list of assoc faces>
	typedef pair<TriEdge, vector<int>> TEdge;
	// list of topo edges
	typedef vector<TEdge> TEdgeList;
	// type for representing topo edge index
	typedef uint16_t EdgeIdx;
	// <char, topo edge idx>
	typedef pair<char, EdgeIdx> ChEidPair; 
	// mapping: face idx -> the sector(topo edge) idx covering the face
	typedef map<unsigned, int> TOrigFace2TopoMap;
public:
	TopoGraph();
	// DO NOT USE DIRECTLY. Use the static versions instead.
	bool is_topo_edge_open(TopoGraph::EdgeIdx _i) const;
	bool has_any_open_topo_edge() const;
	// get the idx of the instance of the given topo edge that has 0 assoc edges
	// return -1:
	// - no such topo edge found
	// - no such instance found
	int getEmptyInstanceIdxOfTEdge(TriEdge _te) const;
	// how many connected components in this graph
	int nConnCmpnt() const;

	friend ostream& operator << (ostream& _os, const TopoGraph& _tg);	
	// return approx. size of this topo graph
	unsigned size();
	// reserve space to facilitate stl operation
	void reserve();
	// free unnecessary space held
	void shrinkToFit();
public:
	static SteinerGraph* stg;
	// return if the topo edge is open
	static bool isOpenSheet( int _vi, int _si );
	// return if the given vertex has any open topo edge
	static bool hasAnyOpenSheet( int _vi );
	// if the vertex has 2d bits in neighborhood
	static bool has2dNbhood( int _vi );
	TopoType type;
	// topo vert list
	vector<int> t_vts;
	// topo edge list
	TEdgeList t_edges;
	// topo vert -> neighbor topo vts
	map<int, vector<int>> t_nbVtsMap;
	// topo edge -> a list of ids pointing to its instance in edge list
	map<TriEdge, vector<EdgeIdx>> t_edgeIdxesMap;
	// graph edge -> the topo vert/edge it's associated with
	TOrigFace2TopoMap t_origFace2TopoMap;
};

// a ternary: <u, v, ti>. 
// <u, v> should be valid graph edge. ti should be valid topo edge of u
typedef trimesh::ivec3 TEdgeTernary;
struct TEdgeTernaryCompare
{
	bool operator () (const TEdgeTernary& _t1, const TEdgeTernary& _t2)
	{
		if (_t1[0] == _t2[0])
			if (_t1[1] == _t2[1])
				return _t1[2] < _t2[2];
			else
				return _t1[1] < _t2[1];
		else
			return _t1[0] < _t2[0];
	}
};

class SteinerGraph
{
public:
	typedef map<TriEdge, vector<int>, EdgeCompare> TriEdgeToStVtsMap; 
	typedef vector<vector<int>> SteinerNbVtsOfVert;

	//
	// the scheme for burning ( & subsequent dualization )
	//
	enum BurnScheme {
		// both original and steiner vts will be burned.
		// later, dual will avoid original triangle vertices.
		ORIGINAL_AND_STEINER,
		// only steiner vts will be burned.
		// dual will pass through original triangle vertices.
		STEINER_ONLY,
		// haven't pick a scheme.
		INVALID
	};

	//
	// how to put dual points on edge and poly region
	//
	enum DualOpt {
		// simply put dual point at center of edge/poly region
		SIMPLE_DUAL, 
		// dualize by solving linear equations
		PSEUDO_INV_DUAL, 
		WEIGHT_CENTER_CLOSENESS_DUAL,
		// dualize by modeling distance function at each point as point source
		WEIGHT_DUAL
	};

	enum FaceField {
		BT2_FINE_TRI,
		BT3_FINE_TRI,
		DIFF_BT2BT3_FINE_TRI
	};

	enum VertField {
		BT2_FINE_TRI_VERT,
		BT3_FINE_TRI_VERT,
		DIFF_BT2BT3_FINE_TRI_VERT
	};

public:
	// construct steiner graph from given graph 
	// by putting steiner points on edges according
	// to the specified scheme
	SteinerGraph(
		const shared_ptr<MyGraph>& _g_ptr, 
		const vector<float>& _dist2surf, 
		SteinerSubdivision::SubdivideScheme _scheme, double _steiner_param,
		int _edgeWeight_idx,
		SteinerGraph::BurnScheme _burn_scheme = BurnScheme::ORIGINAL_AND_STEINER);
	~SteinerGraph();
	//
	// return the burn&dual scheme 
	//
	SteinerGraph::BurnScheme getBurnScheme() const;

	// get the color representing the topo type of vert i's topo graph
	TriColor getTopoColor(int _i);
	// get the topo type of vert i
	TopoType getTopoType(int _i) const;
	// get the assoc.ed faces with the given topo edge of vert v
	void getAssocFaces(int _v, TopoGraph::EdgeIdx _tei, vector<int>& _faces);
	// find the topo edge of _v that _fi mapped to
	// return: 
	// -1 if edge<u, v> / mapping not found; 
	// -2 if topo graphs not constructed yet;
	int mapTopo(int _v, unsigned _fi) const;
	// return whether the vertex has simple topology (i.e. one sector)
	bool simpleTopoAt(int _vi) const;
	// return whether the given topo edge is the leading edge in _v's topo graph
	bool is_leading_topo_edge(unsigned _v, TopoGraph::EdgeIdx _tei) const;
	// return the prev vert for the given vert with the given topo edge
	int getPrevVert(unsigned _v, TopoGraph::EdgeIdx _tei) const;
	int getPrevTopoEdge(unsigned _v, TopoGraph::EdgeIdx _tei) const;
	// analyze and output info. regarding unburnt parts
	// Return a list of bad faces to remove from orig mesh
	vector<int> checkUnburnt(void);
	// is st edge e burnt on assoc face fi? If so, return fire direction _from -> _to
	bool isStEdgeBurnt(const TriEdge& e, int _fi, int& _from, int& _to);
	// return only st. edges that are used as part of some "burn tree"
	void getBurnedStEdgesFromBurnTree(
		vector<std::pair<TriEdge, int> >& _burntE_assocF_list,
		vector< std::pair<size_t,list<TriEdge>> >* _face_with_crossed_paths = nullptr);
	void output_in_mathematica_format(
		const vector< std::pair<size_t,list<TriEdge>> >& _face_with_crossed_paths);
	// compute for each face the 2 intersection points 
	// from the sphere assoc. with each vert.
	// invalid intersections have numeric<float>::max() coordinates
	void computeIntersectionsForFaces(
		const std::shared_ptr<TriMesh> _mesh_orig,
		const vector<float>& _dist2surf);
	// given a kdtree (containing all 3d surface points), an MA face id, 
	// find the 2 "foot points" for that face using the centroid of the face
	// return whether foot points are well defined or not
	util::IsectType find_two_foot_points_for_MA_face(
		const KdTree* _tree, unsigned _fi, 
		TriPoint* _pts
		);
	// for debugging: output info related to foot points of ma faces
	// to a mathematica-readable file
	void output_footPts_info( 
		const vector<TriPoint>& _foot_pts,
		const vector<bool>& _has_valid_isects
		);
	// setup a diffusion system that treats those MA faces without valid foot points as unknowns
	void setupDiffuserForMA( GraphDiffuser& _diffuser );

	///
	/// compute various distance metrics
	///

	// compute the distance to surface for each face using the distance of vertices
	void computeDist2SurfForFaces(const vector<float>& _dist2surf);
	void computeDist2SurfForFaces(const vector<TriPoint>& _closest_2pts);
	// compute the angle-based dist. metric
	void computeAngleMetricForFaces(
		const vector<TriPoint>& _isects_per_face, 
		vector<float>& _angle_per_face);
	// compute the lambda dist. metric
	void computeLambdaMetricForFaces(
		const vector<TriPoint>& _isects_per_face, 
		vector<float>& _lambda_per_face);
	// compute the geodesic metric
	void computeGeodMetricForFaces(
		shared_ptr<TriMesh> _m3d,
		const vector<TriPoint>& _isects_per_face, 
		vector<float>& _geod_per_face, 
		string& _cache_file);
	// compute the difference btw. dist2surf and burn dist for each face
	// store result in vector distFace
	void computeDiffDist();
	// return faces whose diff. dist. is above threshold
	void filterFaceByDiffDist(float _t, vector<unsigned>& _faces);
	// create a dual graph for the planar subdivision obtained from 
	// all burn paths and original faces edges
	void dualizeFacesWithBurnPaths(
		vector<std::pair<TriEdge, int> >& _burntE_assocF_list,
		DualOpt _edge_dual_opt,
		DualOpt _poly_dual_opt,
		float _eps);
	// query if dual is valid to use
	bool isDualValid() const;
	// return if a given dual vert id corresponds to a face-dual.
	bool isFaceDual(int _dual_id) const;
	// return the list of dual vts
	const vector<TriPoint>& getDualVts() const;
	// burn a polyline network
	// start burning from degree-1 vts (boundary vts)
	// whose starting burn time could be passed via _start_time map
	void burnMedialCurveNetwork(
		bool _start_time,
		bool _stop_at_junction,
		bool _protect_bt2,
		bool _under_estimate_dist);
	// relocate face dual to the connected edge dual that's last burned away
	// NOTE: a burning of curve (burnMedialCurveNetwork) must have been performed 
	// before calling this function
	void relocate_face_dual();
	// build neighbor-ship given the vts of edges of a line network
	void build_neighborship(
		const vector<TriPoint>& _vts, const vector<TriEdge>& _edges, 
		vector<vector<unsigned>>& _neighbors,
		vector<int>& _nbs_cnt);
	// perturb the original vertices (keep the relative position of steiner vts) by the given amount
	// store the perturbed vts into the given list
	void perturbVts(vector<TriPoint>& _pert_vts, float _pert_amount);
	// given the network of dual lines, 
	// find a cycle bases and output them & related info to the given ostream
	// used for debugging dual line network.
	void debug_dualline_cycles(
		const unsigned _vts_size, 
		const vector<TriEdge>& _dual_edges,
		vector<int>& _from_which_face_dual_edge,
		bool _output_smallest_cycle = false, // output the smallest of the cycle base?
		std::ostream* _os = &std::cout		// where to output the cycle info?
		);
	// find a cycle base for the given burn line network
	// output them to the given ostream
	void debug_burnline_cycles(
		const vector<TriPoint>& _burn_vts,
		const vector<std::pair<TriEdge, int> >& _burntE_assocF_list,
		bool _output_smallest_cycle = false, // output the smallest of the cycle base?
		std::ostream* _os = &std::cout		// where to output the cycle info?
		);
	// output all connected components of the dual line network to the given ostream
	void debug_dualline_conn_comp(
		const unsigned _vts_size,
		const vector<TriEdge>& _dual_edges,
		vector<int>& _from_which_face_dual_edge,
		bool _output_smallest_cc = false, // output the smallest of the cycle base?
		std::ostream* _os = &std::cout		// where to output the cycle info?
		);
	// find a cycle base given the line network
	void find_cycle_base(
		const unsigned _network_vts_size,
		const vector<TriEdge>& _network_edges, 
		list<vector<unsigned> >& _cycles
		);
	// find all connected components in the line network
	void SteinerGraph::find_conn_cmpnts(
		const unsigned _vts_size, 
		const vector<vector<int>>& _adjacency, 
		vector<vector<TriEdge> >& _cc_eList, 
		vector<vector<unsigned> >& _cc_vList
		);
	void find_conn_cmpnts(
		const unsigned _vts_size,
		const vector<TriEdge>& _network_edges, 
		vector<vector<TriEdge> >& _cc_list, vector<vector<unsigned> >& _cc_vList
		);

	// print a few vertices' burn time and return those chosen vertices' indices
	void printBurnTime(vector<size_t>& _vts_indices);

	/* get finner triangulation related info */
	//
	// given an edge of a tri in dual-induced fine triangulation
	// return true: a dual edge, false: part of orig. tri edge
	//
	inline bool isDualEdgeInFineTri(const TriEdge& _e) const
	{
		int tmp_dual_id;
		return isDualVertInFineTri(_e[0], tmp_dual_id) && isDualVertInFineTri(_e[1], tmp_dual_id);
	}

	//
	// given a vert of a tri in finner triangulation, 
	// return true: a dual point; false: a st vert. and its vert id in the respective vert list
	//
	inline bool isDualVertInFineTri(int _u, int& _v) const
	{
		if (_u >= m_stSubdiv.sizeOfVts())
		{
			_v = _u - m_stSubdiv.sizeOfVts();
			return true;
		}
		else
		{
			_v = _u;
			return false;
		}
	}
	//
	// return the id of the orig face that contains the given fine tri induced by duals
	//
	inline int whichOrigFaceIsDualFineFaceFrom( int _fine_fi ) const
	{
		return m_origTri_for_finerTri_byDual[ _fine_fi ];
	}
	inline TriPoint getVertPosOfFineTri(int _u) const
	{
		int vi;
		if ( isDualVertInFineTri(_u, vi) )
			return this->dual_vts[vi];
		else
			return this->m_stSubdiv.getVert(vi);
	}
	inline float getScalarOfFineTriVert(VertField _scalar_field, int _u, int _fi) const
	{
		int vi;
		int orig_fi = m_origTri_for_finerTri_byDual[_fi];
		float scalar = 0.0f;

		switch (_scalar_field)
		{
		case BT2_FINE_TRI_VERT:
			if ( isDualVertInFineTri(_u, vi) )
			{
				auto find_it = m_bt2_MC_perFace.find(trimesh::ivec2(vi, orig_fi));
				scalar = find_it->second;
			}
			else
			{
				scalar = bt2MA_vert_per_sheet[vi][mapTopo(vi, orig_fi)];
			}
			break;
		default:
			scalar = 0.0f;
			break;
		}

		return scalar;
	}
	inline float getScalarOfFineTri(FaceField _scalar_field, int _fi) const
	{
		int orig_fi = this->m_origTri_for_finerTri_byDual[_fi];
		const auto& f = m_finerTri_byDual[_fi];
		int v;
		bool is_dual_vert; 
		float scalar = 0.0f;

		switch (_scalar_field)
		{
		case BT2_FINE_TRI:
			for (int i = 0; i < 3; ++i)
			{
				is_dual_vert = isDualVertInFineTri(f[i], v);
				if ( is_dual_vert ) // a dual vert
				{
					scalar += m_bt2_MC_perFace.find(trimesh::ivec2(v, orig_fi))->second;
				}
				else // a steiner vert
				{
					int v_te = mapTopo(v, orig_fi);
					scalar += bt2MA_vert_per_sheet[v][v_te];
				}
			}
			scalar /= 3.0f;
			
			break;
		default:
			scalar = 0.0f;
			break;
		}

		return scalar;
	}
	/// set the state of the given face id to "pruned"
	inline void setFacePruned(int _fi, bool _pruned)
	{
		this->m_faceActiveFlag[_fi] = !_pruned;
	}
	/// query "pruned" status for the given face id
	inline bool isFacePruned(int _fi) const
	{
		return !this->m_faceActiveFlag[_fi];
	}

	/* used for comparison against other methods */
	bool outputToWenping(const char* _file_name, const vector<float>& _radii);

	// export bt of original vertices to the specified file
	void exportPerVertexET(std::string _file_name);
	//bool setPerVertexBT(std::string _file_name);
	// exports bt of each sector of each original vertex.
	// the file format is as follows:
	// each line corresponds to a face, followed by 3 of its vertices
	// each vertex is followed by the burntime of the sector of that vertex covering this face.
	// Sure there is redundancy, but this is the quickest and easiest way to structure this kind of info, 
	// which can be read back and rendered easily as well.
	void exportPerSectorET(std::string _file_name);
	//bool setPerSectorBT(std::string _file_name);

private:
	/* steiner graph helpers 
	*/
	// subdivision dispatcher
	void steinerSubdivide(double _steiner_param, SteinerSubdivision::SubdivideScheme _scheme);
	// compute dist-to-surf for steiner vertices and append them to the given list
	// a vector of |vts| 0s will be returned if _dist2Surf is empty
	void computeBT3ForSteinerVts(const vector<float>& _bt3_orig);
	
	/* topo graphs helpers 
	*/	
	// build topo graph at each vertex
	void buildTopoGraphs();
	// get the topo edge idx of the given vert that covers the given face
	// vi has to be either an orig vert of f, or a st vert on one of the edges of f
	/*bool getAssocTEdgeOfVertForFace(const unsigned _vi, const TriFace& _f, TopoGraph::EdgeIdx& _tei);*/
	// identify topo graph type
	void identifyAll();
	TopoType identify(const TopoGraph& _tg);
	// update topo graph stats
	void resetTopoStats();
	void updateTopoStats(TopoType _type);
	void printTopoStats();
	// manage topo graph allocation inside the topo list
	TopoGraph* make_topo(int _vi, TopoGraph* _tg_ptr);
	TopoGraph* make_topo(int _vi);
	void del_topo(int _vi);
	
	/* 
	** burning function helpers 
	*/	
	// set the edge weight function used during burning
	void setEdgeWeightFunc(int _edgeWeight_idx);
	// the edge weight function in effect(set according to specified edge-weight-idx)
	double (SteinerGraph::*edgeWeight)(unsigned _ui, unsigned _vi);
	///
	/// available edge weight function candidates
	///
	inline double euclidean_edge_weight(unsigned _ui, unsigned _vi)
	{
		return trimesh::dist(m_stSubdiv.getAllVts()[_ui], m_stSubdiv.getAllVts()[_vi]);
	}
	// fire goes from u -> v
	inline double bt2bt3diff_edge_weight(unsigned _ui, unsigned _vi)
	{
		return std::max(
			euclidean_edge_weight(_vi, _ui) - (radii_vert[_vi] - radii_vert[_ui]), 
			0.0 
			);
	}
	// burn this graph-represented medial axis
	void burn();
	// triangulate the MA mesh to finner tri faces
	void triangulate();
	//
	// record the mapping: steiner edge -> the dual on it
	// use case: hybrid skeleton
	//
	//void recordDualwrtStEdge(unsigned _fi, TriEdge _st_e, unsigned _e_dual_vi, unsigned _poly_dual_vi);

private:
	typedef CGAL::Quotient<CGAL::MP_Float>							NT;
	typedef CGAL::Cartesian<NT>										Exact_K;
	typedef CGAL::Point_2<Exact_K>									Point2_E; // exact cgal point_2
	typedef CGAL::Arr_segment_traits_2<Exact_K>						Arr_seg_traits;
	typedef Arr_seg_traits::Point_2									Arr_point2; // point in 2d arr.
	typedef Arr_seg_traits::X_monotone_curve_2						Arr_seg2; // segment in 2d arr.
	typedef CGAL::Arrangement_2<Arr_seg_traits>						Arrangement2; 
	typedef CGAL::Arr_walk_along_line_point_location<Arrangement2>	Walk_pl; 
	typedef Arrangement2::Vertex_handle								Arr_vert_hdl;
	typedef Arrangement2::Halfedge_handle							Arr_halfedge_hdl;
	typedef Arrangement2::Face_handle								Arr_face_hdl;

private:
	// do initialize chores
	void init();
	// reset statistics related to burn function
	void resetVtsBurnStats();
	// for each topo edge, merge its assoc faces to _assoc
	void mergeAssocFaces(
		// merge associated faces in here
		vector<int>& _assoc, 
		// indices of topo edges whose assoc.ed faces with be merged
		const vector<int>& _idxes, 
		// topo edges that _idxes index into
		const TopoGraph::TEdgeList& _tedges
		);
	// finalize the given topo edge and other topo edges that are "openned" by it.
	void openTopo(
		const TopoGraph& _tg, // within this topo graph:
		const TopoGraph::EdgeIdx _te, // the topo edge to final
		vector<bool>& _final, // the "final" flag list for tg
		vector< std::pair<TopoGraph::EdgeIdx,int> >& _fnIdxes // (return) a list of edgeIdxes finalized
		);
	// initialize boundary vts start burn dist to given dist value (e.g. balls radii)
	// if _dist2surf is emtpy, 0.0f is assigned to each boundary vert.
	void assignBoundaryStartDist(const vector<float>& _dist2surf/*, vector<int>& _bndry*/);
	//
	// critical logic wrapped here:
	// given a vertex id and its current burntime value, 
	// return the updated burntime that's at least radius at that vertex (the protection)
	//
	float update_with_protection( int _vi, float _bt ) const;
	//
	// steiner_only burnscheme helper:
	// recover original vts' burn times from those of st vts
	//
	void recover_bt_for_orig_vts();
	// <sum of dual vts assoc. w. an edge, count, dual index for that edge in dual_vts>
	struct movingSum_edgeDual_struct
	{
		movingSum_edgeDual_struct() 
			: count(0) {};
		TriPoint sum; // sum of dual coordinates
		unsigned count; // count of dual added to sum
		float bt2; // the bt2 distance for corresponding edge dual
		unsigned dual_v_id; // the id representing the only instance in dual_vts
	};
	// assign a burn direction for boundary vertices
	// use this function to set up these directions 
	// if there are no dist2surf info available
	void findBoundaryBurnDir(
		vector<trimesh::vec>& _bndry_dirs,
		map<unsigned, unsigned>& _vIdx_burnDirIdx_map);
	// find the dual point of the given poly subdivision 
	// and the edges of each poly
	// using a very simple strategy
	void dualize_poly_subdivision_simple(
		unsigned _fi, 
		const std::vector<vector<unsigned>>& _poly_subdivision,
		const map<unsigned, std::pair<trimesh::vec, float>>& _vert_burnDirDistPair_map,
		const map<TriEdge, vector<int> >& _burntEdge_assocFaces_map, 
		map<TriEdge, movingSum_edgeDual_struct>& edge_dual_map,
		DualOpt _dual_method,
		bool _verbose);
	// find the dual point of the given polygon face
	// modeled and solved as a quadratic opt. problem
	void dualize_poly_subdivision(
		const std::vector<unsigned>& _vts, unsigned _fi, 
		const std::vector<vector<unsigned>>& _poly_subdivision,
		const map<unsigned, std::pair<trimesh::vec, float>>& _vert_burnDirDistPair_map,
		const map<TriEdge, vector<int> >& _burntEdge_assocFaces_map,
		map<TriEdge, movingSum_edgeDual_struct>& _edge_dual_map,
		DualOpt _edge_dual_opt = WEIGHT_CENTER_CLOSENESS_DUAL,
		DualOpt _poly_dual_opt = WEIGHT_CENTER_CLOSENESS_DUAL,
		float _w = 0.08f, float _eps = 1e-7f,
		bool _verbose = false, 
		ostream* _os = NULL );
	//
	// create dual vts at orig vts' location
	void process_orig_vts_as_duals();
	// find the incoming burn dir. and distance 
	// for the given st vts on a tri face
	void find_burn_dir_dist_of_face(
		const vector<unsigned>& _vts, unsigned _fi, 
		// burn dir can be retrieved directly from this map
		const vector<trimesh::vec>& _bndry_dirs,
		const map<unsigned, unsigned>& _vIdx_burnDirIdx_map,
		map<unsigned, std::pair<trimesh::vec, float>>& _vert_burnDirDistPair_map,
		bool _verbose );
	// finalize dual vts
	void finalize_dual_vts(
		map<TriEdge, movingSum_edgeDual_struct>& _edge_dual_map);
	// return if undirected edge <u, v> is burnt or not w.r.t. the give face _fi
	bool is_edge_dualizable(
		unsigned _u, unsigned _v, unsigned _fi,
		const map<TriEdge, vector<int> >& _burntEdge_assocFaces_map,
		bool _verbose);
	// solve the linear equation _ax=_b robustly using SVD
	// min ||_c-_x||^2
	// if _w not NULL, use it to weight closeness of _x to the center _c
	template<std::size_t D>
	void linearSolveSVD(
		const Eigen::Matrix<double,D,D>& _a, const Eigen::Matrix<double,D,1>& _b, 
		const Eigen::Matrix<double,D,1>& _c, 
		Eigen::Matrix<double,D,1>& _x)
	{
		typedef Eigen::Matrix<double,D,D> MatrixDD;
		double eps = 1e-2;
		Eigen::JacobiSVD<MatrixDD> svd(_a, Eigen::ComputeFullU | Eigen::ComputeFullV);
		auto singular_vals = svd.singularValues();
		MatrixDD singular_diag = MatrixDD::Zero();
		for (unsigned i = 0; i < D; ++i )
		{
			if (singular_vals(i) / singular_vals(0) > eps)
				singular_diag(i, i) = 1.0 / singular_vals(i);
			else
				singular_diag(i, i) = 0.0;
		}

		/*_x = _c + svd.matrixV()*singular_diag*(svd.matrixU().transpose())*(_b-_a*_c);*/
		/*_x = svd.matrixV()*singular_diag*(svd.matrixU().transpose())*_b;*/
		MatrixDD a2 = svd.matrixV()*singular_diag*(svd.matrixU().transpose());
		_x = (MatrixDD::Identity()-a2*_a)*_c + a2*_a*a2*_b;
	}
	// snap in place the given list of weights (sum up to 1) to the range [0, 1]
	inline void snap_to_boundary(vector<float>& _w)
	{
		float sum = 0.0f;
		for (auto w_it = _w.begin(); w_it != _w.end(); ++w_it)
		{
			*w_it = std::max(0.0f, *w_it);
			sum += *w_it;
		}
		for (auto w_it = _w.begin(); w_it != _w.end(); ++w_it)
			*w_it /= sum;
	}

public:
	// the orig graph corresponding the input MA mesh
	shared_ptr<MyGraph> m_origG;
	// the steiner subdivision structure 
	SteinerSubdivision m_stSubdiv;
	/*vector<TriPoint> st_vts;
	vector<TriEdge> st_edges;
	MyGraph::WeightMap st_weights;*/
	SteinerNbVtsOfVert st_nbVtsOfVert;
	/*TriEdgeToStVtsMap triE_stVts_map;*/
	// topo graphs. 
	// one topo graph for each orig vertex
	// one topo graph for all st vts on the same edge
	vector<TopoGraph*> t_graphs;
	// since steiner vertices on the same edge share only one topo graph
	// we need to know whether actual topo graph is held at a orig/st vertex vi
	// if true, t_graphs[vi] points to a part of memory holding a topo graph
	vector<bool> t_holds_actual_storage;

	/// Important states that affect general behavior of the steiner graph
	// how will burning and dualization perform
	BurnScheme m_curBurnScheme;

	/// burning function related
	enum EdgeWeight {EW_EUCLIDEAN, EW_BT2BT3DIFF};
	// current edge weight of choice
	EdgeWeight m_curEdgeWt;
	//  for each vert final status of topo edge (id)
	vector<vector<bool>> final;
	//  for each vert the topo edge (id) that captures its cur burn path
	vector<TopoGraph::EdgeIdx> min_tedge;
	//  for each vert the dist from every sheet (topo edge id)
	vector<vector<float>> bt2MA_vert_per_sheet;
	//  for each vert its burn dist 
	vector<float> bt2MA_vert;
	//  for each vert prev vert on each sheet (topo edge id)
	vector<vector<int>> prev_vert;
	//  for each vert the prev sheet that fire comes from
	vector<vector<TopoGraph::EdgeIdx>> prev_tedge;
	//  for each vert its burnt-or-not status
	vector<bool> burnt;
	//  for each vert |bt1 - bt2|
	vector<float> bt2bt3diffMA_vert;
	// for each face/vert its dist. to surface
	vector<float> bt3MA_face;
	vector<float> bt3MA_vert;
	vector<float> radii_vert; // a radius for each vert of MA
	// for each face its burn dist.
	vector<float> bt2MA_face;
	// for each face |burn dist - dist2surf|
	vector<float> bt2bt3diffMA_face;
	// burn distance stats
	float min_burnDistOrigVts, max_burnDistOrigVts;
	float min_burnDistAllVts, max_burnDistAllVts;
	float min_diffDistOrigFaces, max_diffDistOrigFaces;

	/// medial curve related
	// whether dual structure has been computed
	bool m_dual_valid;
	// vts of the dual structure
	vector<TriPoint> dual_vts;
	// edges of the dual structure
	vector<TriEdge> dual_edges;
	// the face from which a dual edge is generated (-1 for N/A)
	vector<int> m_dualEdge_from_which_face;
	// dual vert -> the tri edge it's on (-1 for face dual)
	vector<int> m_dualV_triEdgeIdx;
	// -1: not face dual (edge dual)
	// >= 0: the original tri face the face dual is from
	vector<int> m_is_face_dual;
	// the prev vert of each vert of medial curve determined after burning
	vector<int> burnPrev_medialCurve;
	// the next vert of each vert determined after burning. by default -1.
	vector<int> burnNext_medialCurve;
	// the dist at each vert of medial curve after burning
	vector<float> bt1_medialCurve;
	// the dist at each vert of medial curve to the boundary of the medial curve
	vector<float> bt2_MC;
	// <edge dual id, assoc tri face id> -> the bt2 w.r.t. that face
	map<trimesh::ivec2, float> m_bt2_MC_perFace;
	// the dist at each MC vert to the orig 3D surface
	vector<float> bt3_MC;

	/// possible use cases for the following: Hybrid skeleton
	// a st edge -> the dual point on it
	vector<int> m_stEdgeIdx_dual_map;
	// TODO: remove it! <edge dual, face id> -> poly dual that's connected to the edge dual
	//map<pair<unsigned, unsigned>, unsigned> m_eDualFaceIdPr_polyDual_map;

	// face flag indicating whether it's "active" or not (possibly meaning it's not pruned, etc.)
	vector<bool> m_faceActiveFlag;

	// finer triangulation induced from st vts and face dual vts (generated in the dualization process)
	// the vts of this triangulation will be: Steiner + face duals. 
	// their indices within respective list will stay untouched.
	// only the new tri faces need to be tracked here. 
	// each face simply a steiner edge + a poly face dual
	// Use case: 
	// - TODO: construct 2-complex in skeleton
	// - in *video section*, for smooth display of iso-contour on MA 
	//   (TODO: get rid of this use case and just use the m_MA_finer_tris instead!!!
	//   see the TODO note in IsoContourOnMA.cpp)
	vector<TriFace> m_finerTri_byDual;
	vector<int> m_origTri_for_finerTri_byDual;
	// mapping dual edge id -> the fine tri it's from
	// Used by:
	// - pruneMedialCurveConstrainedByIsoContour()
	// - uploadFromRemainedMC()
	vector<int> m_fromFineTri_for_dualE;

	// Store any triangulation of the original MA faces here
	// along with the face id of the original face for each fine tri
	// Use cases: during fine-pruning, to extract iso-contour on MA
	vector<std::pair<TriFace, unsigned>> m_MA_finer_tris;

	///
	/// competing measures related (e.g. angle, lambda, etc)
	///
	// store foot points of each MA face
	vector<TriPoint> m_footPtsForMAFaces;
	// indicating for each face whether foot points are valid or not
	vector<bool> has_isects;

	/* private data member */
private:
	// topo type for each vertex
	vector<TopoType> t_type;
	// num of sectors (topo edges) at each vertex.
	// If a vertex has only one sector, then after its topo type is determined,
	// we can delete its topo graph. 
	// Supposedly this will save a big chunk of memory.
	vector<int> t_n_sector;
	// num of each type of vertices
	vector<int> t_cnt;
};

#endif