#ifndef HYBRID_SKELETON_H
#define HYBRID_SKELETON_H

#include <map>
#include <set>
#include <queue>
#include "allTypes.h"
#include "origGraph.h"
#include "steinerSubdivision.h"
#include "geometryRepresentation.h"

using namespace std;

//
// utitility class for computing a subdivision of a tri face
// by its dual vts and orig vts and the edges among them
//
class HybridSkeleton
{
public:
	struct simple_pair
	{
		simple_pair(int _t, unsigned _idx0, unsigned _idx1)
		{
			type = _t;
			idx0 = _idx0;
			idx1 = _idx1;
		}
		simple_pair(const simple_pair& _spr)
		{
			type = _spr.type;
			idx0 = _spr.idx0;
			idx1 = _spr.idx1;
		}
		simple_pair& operator = (const simple_pair& _spr)
		{
			type = _spr.type;
			idx0 = _spr.idx0;
			idx1 = _spr.idx1;

			return *this;
		}
		// type = 1: face-edge pair, idx0 is face id, idx1 is edge id
		// type = 0: edge-vert pair, idx0 is edge id, idx1 is vert id
		int type; 
		unsigned idx0;
		unsigned idx1;
	};

public:
	HybridSkeleton();
	~HybridSkeleton();

	/* initialize h.s. states */
	void setup(
		const std::shared_ptr<SteinerGraph>& _stg
		);
	/* clear current per-face-subdiv */
	void reset();
	
	/**/
	void preprocess();
	///*prepare for the face extraction function*/
	//void prepareFaceExtraction(
	//	const vector<TriPoint>& dual_vts, 
	//	const vector<TriEdge>& dual_edges);
	///* derive poly faces from subdivision */
	//void extractPolyFaces();

	/* get vts coordinates, edges, and tri faces of H.S.
	   Note: edge and face are list of vertex indices */
	const vector<TriPoint>& getVts() const;
	const vector<TriEdge>& getEdges() const;
	const vector<TriFace>& getTriFaces() const;
	const vector<bool>& getEdgeRemovedTags() const;
	const vector<bool>& getFaceRemovedTags() const;
	//const vector<unsigned>& getTriFromPolyID() const;

	void getRemainedVts(vector<int>& _vts_ids) const;
	void getRemainedEdges(vector<TriEdge>& _edges) const;
	void getRemainedFaces(vector<TriFace>& _tris) const;
	void getDualEdges( vector<TriEdge>& _edges ) const;
	// get per edge color. orig edges and dual edges are assigned color
	// from different space
	void getEdgeDiffValColor(float* edge_color);
	void getEdgeRelDiffValColor(float* edge_color);
	// get per face color
	void getTriColor(float* tri_color);

	/* assign values to elements: faces, edges
	   note: no need to maintain values for vertices 
	   input values contain information required to derive values for 
	   edges and faces elements (1, 2-cell)
	 */
	void assignElementValues(
		const vector<vector<float>>& _orig_persheet_bt2bt3, 
		const vector<vector<float>>& _orig_persheet_bt2bt3rel, 
		const vector<float>& _orig_vts_bt2bt3,
		const vector<float>& _orig_vts_bt2bt3rel,
		const vector<float>& _dual_vts_bt2bt3,
		const vector<float>& _dual_vts_bt2bt3rel,
		const vector<float>& _dual_vts_bt1bt2,
		const vector<float>& _dual_vts_bt1bt2rel
		);
	/* 
	set the threshold for detecting degenerate poly faces
	*/
	void setPolyFaceDegenerateThresh(float _t);
	void setComponentFaceNumberThresh(int _t);

	/* prune HS in a cell-complex retraction manner
	   remaining edges, faces are returned in the lists
	*/
	void prune(
		float _f_diff_r, float _f_reldiff_r, 
		float _l_diff_r, float _l_reldiff_r,
		bool _use_inputs_directly,
		bool _remove_small_components
		);

	/* export the currently remaining parts of the skeleton to file
	*/
	void exportSkeleton(std::string _skel_name, const trimesh::xform& _transform);

	/* transforming threshold: given ratio -> absolute value */
	void absToRatio(
		float _abs_f_diff_t = 0.0f,
		float* _f_diff_r = nullptr,
		float   _abs_f_reldiff_t = 0.0f,
		float* _f_reldiff_r = nullptr,
		float   _abs_l_diff_t = 0.0f,
		float* _l_diff_r = nullptr,
		float   _abs_l_reldiff_t = 0.0f,
		float* _l_reldiff_r = nullptr
		);
	/* transforming threshold: absolute value -> given ratio */
	void ratioToAbs(
		float _f_diff_r = 0.0f,
		float* _abs_f_diff_t = nullptr,
		float _f_reldiff_r = 0.0f,
		float* _abs_f_reldiff_t = nullptr,
		float _l_diff_r = 0.0f,
		float* _abs_l_diff_t = nullptr,
		float _l_reldiff_r = 0.0f,
		float* _abs_l_reldiff_t = nullptr);
private:
	//
	// convert a edge-rep tri to a vert-rep
	//
	void to_vts_list( const TriFace& _tri_e_rep, int* _vts_l ) const;
	//
	// convert a list of edges to a list of verts
	//
	void to_vts_list( const vector<unsigned>& _edge_l, vector<unsigned>& _vts_l) const;
	//
	// if the given edge id corresponds to a st edge on orig tri edge
	// 
	bool part_of_orig_edge(unsigned _ei) const;
	//
	// given a vert-id in fine-tri, map it to the correct id within vts list of the skeleton
	//
	int transform_fine_vert_id( int _vi ) const;
	//
	// given an isolated edge in the orig graph, map the vts indices so that they refer to cur v-list
	TriEdge transform_isolate_edge( const TriEdge& _e ) const;
	//
	// given a vert-id of skeleton, map it to within the vts range of the fine-triangulation in stg
	//
	int transform_skel_vert_id( int _vi ) const;
	//
	// given a skeleton edge, map it so that vts indices are in the vts range of the fine-triangulation in stg
	//
	inline TriEdge transform_skel_edge( const TriEdge& _e ) const
	{
		return TriEdge( transform_skel_vert_id( _e[ 0 ] ), transform_skel_vert_id( _e[ 1 ] ) );
	}
	//
	// if edge has a value below threshold
	//
	bool edge_below_threshold(unsigned _ei, 
		float _bt2bt3_t, float _bt2bt3rel_t, float _bt1bt2_t, float _bt1bt2rel_t);
	//
	//
	//
	bool face_edge_below_threshold(unsigned _fi, unsigned _ei, 
		float _bt2bt3_t, float _bt2bt3rel_t, float _bt1bt2_t, float _bt1bt2rel_t);
	//
	//
	//
	void prune_while_iteration(
		const set<unsigned>& vts_to_debug,
		float bt2bt3_t, float bt2bt3rel_t,
		float bt1bt2_t, float bt1bt2rel_t,
		std::queue<simple_pair>& q
		);
	void mark_components(
		const vector<bool>* const _removed, 
		const vector<int>& _ref_edge, const vector<int>& _ref_vert, 
		const vector<vector<unsigned>>& _nb_faces, const vector<vector<unsigned>>& _nb_edges,
		float _cmpnt_geomSize_t,
		int _num_face_t
		);

	//
	// debug routines
	//
	void print_edge( int _ei ) const;

private:
	// basic information H.S. needs from the outside
	std::shared_ptr<SteinerGraph> m_stg;
	SteinerSubdivision* m_stSubdiv;
	MyGraph* m_orig_g;

	/// orig edge id -> the dual vts on the edge
	//vector<vector<unsigned>> m_triEdgeIdx_dualVts_map; // TODO: get rid of this.

	vector<TriPoint> m_vts;
	vector<TriEdge> m_edges;
	//map<TriEdge, unsigned> m_edge_id_map; // TODO: get rid of this.
	//vector<int> m_edge_fromOrigEdgeID; // TODO: get rid of this.
	/// the polygon faces from subdivision (deprecated)
	//vector<vector<unsigned>> m_poly_faces; // TODO: get rid of this.
	//vector<unsigned> m_poly_fromFaceID; // TODO: get rid of this.
	//vector<unsigned> m_tri_fromPolyID; // TODO: get rid of this.

	/// the collection of faces (triangles, from triangulation induced by duals)
	vector<TriFace> m_tri_faces;
	vector<float> m_face_area;
	
	/* information for pruning */
	// threshold under which faces are considered as degenerate
	float m_face_degen_thresh;
	// threshold above which a component will be considered as *too big*
	int m_cmpnt_num_faces_thresh;
	// values used in pruning for HS elements
	vector<float> face_bt2bt3;
	vector<float> face_bt2bt3rel;
	vector<float> edges_measure;
	vector<float> edges_rel_measure;
	float min_bt2bt3, max_bt2bt3;
	float min_bt2bt3rel, max_bt2bt3rel;
	float min_bt1bt2, max_bt1bt2;
	float min_bt1bt2rel, max_bt1bt2rel;
	// remove tag for face [1], edge [0], and vertex [2]
	vector<bool> m_removed[3];
	// ref count and neighboring information for verts and edges
	// model specific ref count
	// doesn't change from prune to prune
	vector<int> m_ref_vert_const; 
	vector<int> m_ref_edge_const;
	vector<vector<unsigned>> m_nb_faces_for_edge;
	vector<vector<unsigned>> nb_edges_for_vert;
	// prune-specific ref count
	// valid after each prune
	vector<int> m_ref_vert_per_prune; 
	vector<int> m_ref_edge_per_prune;

	// remove flag for faces
	vector<bool> m_to_remove_face;
};
#endif