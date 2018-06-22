// Generate & display iso-surface of the distance function 
// defined by a 2-manfifold triangle mesh
// 
// Copyright (C) 2018 Yajie Yan <danielyan86129@hotmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef ISOSURF_FROM_2_MANIFOLD_H
#define ISOSURF_FROM_2_MANIFOLD_H

#include "allTypes.h"
#include "commonUtils.h"
#include "origGraph.h"
#include "drawable.h"
#include "geometryRepresentation.h"
#include <boost/heap/fibonacci_heap.hpp> // implementaion for fibonacci queue
#include <map>

using std::map;
class IsoSurfFrom2Manifold
{
public:
	IsoSurfFrom2Manifold ();
	~IsoSurfFrom2Manifold();

	void reset();
	/* init current iso-surf by passing in the steiner graph struct, 
	the 2 manifold which is to be deformed to obtain the iso-surf, 
	and the mesh drawer object */
	void setup(
		const std::shared_ptr<SteinerGraph>& _stg_ptr, 
		const std::shared_ptr<TriMesh>& _3d_surf, 
		std::shared_ptr<Drawable>& _surf_drawer);
	/* precompute: 
	- correspondence between orig vert and MC/MA vert */
	void precompute();
	/* given iso value, update the current iso-surf.
	Update vertices' position, and upload data to gpu for render*/
	void updateIsoSurf(float _iso_ratio, bool _hide_snapped_faces);
	void updateIsoSurf(float _iso, bool _hide_snapped_faces, int);
	void refineTriangulation(bool _propagate);

private:
	/// face refinement related
	// edge_err_tuple
	typedef tuple<TriEdge, float> edge_err_pair;
	// give face with longer longest-edge higher priority
	struct edge_err_pair_comp
	{
		bool operator () (
			const edge_err_pair& _a, 
			const edge_err_pair& _b) const
		{
			return std::get<1>(_a) < std::get<1>(_b);
		}
	};
	typedef boost::heap::fibonacci_heap<edge_err_pair, boost::heap::compare<edge_err_pair_comp> > fib_heap;
	// used to prioritize faces to refine
	fib_heap m_refine_q;
	//vector<faceId_edge_edgeId> m_refine_q;
	// vertex types (also defines their priority)
	// for original tri vert, steiner vert, edge dual, face dual
	static const int vType_priority[4];
	// handler entry for each pair of possible vert types
	enum vTypePair {
		ORIG_ORIG = 0x33, ORIG_STEINER = 0x32, ORIG_EDUAL = 0x31, ORIG_FDUAL = 0x30,
		STEINER_STEINER = 0x22, STEINER_EDUAL = 0x21, STEINER_FDUAL = 0x20,
		EDUAL_EDUAL = 0x11, EDUAL_FDUAL = 0x10,
		FDUAL_FDUAL = 0x00
	};
	inline int get_vert_type(int& _u);
	inline int encode_type_for_vts(int& _u, int& _v);
	typedef bool (IsoSurfFrom2Manifold::*vTypePairHandler)(const TriEdge& _e);
	//static vTypePairHandler vType_pair_entry[0x33+1];

	// handler for each pair of possible vert types
	// return whether the pair of vts are on the same edge/face of MA
	inline bool same_face_orig_orig(int _u, int _v);
	inline bool same_face_orig_steiner(int _u, int _v);
	inline bool same_face_orig_edual(int _u, int _v);
	inline bool same_face_orig_fdual(int _u, int _v);
	inline bool same_face_steiner_steiner(int _u, int _v);
	inline bool same_face_steiner_edual(int _u, int _v);
	inline bool same_face_steiner_fdual(int _u, int _v);
	inline bool same_face_edual_edual(int _u, int _v);
	inline bool same_face_edual_fdual(int _u, int _v);
	inline bool same_face_fdual_fdual(int _u, int _v);

	// update iso surf geometry and upload to gpu
	void update_iso_surf(bool _hide_snapped_faces);
	// estimate iso-surface edge error
	float estimate_error(const TriEdge& _e);
	// whether edge needs to split
	bool need_split_1(const TriEdge& _e);
	bool need_split_2(float _err);
	bool need_split_3(const TriEdge& _e);
	// whether face needs to split
	bool need_refine(const TriFace& _f, TriEdge& _e_to_split);
	// return the new vertex used to split the edge.
	TriPoint split_isoSurf_edge(const TriEdge& _e);
	// return the new vertex used to split the orig surf edge
	TriPoint split_origSurf_edge(const TriEdge& _e);
	// refine the face using the _vidx on edge _e. return the two new faces.
	void refine_face(const TriFace& _f, const TriEdge& _e, int _vidx, TriFace& _f1, TriFace& _f2);
	// update the corresp. for iso surf vert v_idx
	void addTo_closest(int _v_idx, const  TriEdge& _e);
	// update the nb face for edge e
	void update_nb_face(TriEdge _e, int _old_f, int _new_f);

	/// low-level helpers
private:
	// is _u on edge <_eu, _ev> or do they form a valid tri face in MA?
	inline bool on_same_edge_or_face(int _u, int _eu, int _ev);
	inline bool same_edge_or_on_same_face(const TriEdge& _e1, const TriEdge& _e2);
	inline bool is_face_vertex(int _u, const TriFace& _f);
	inline bool is_face_edge(const TriEdge& _e, const TriFace& _f);

private:
	std::shared_ptr<SteinerGraph> m_stg_ptr;
	std::shared_ptr<TriMesh> m_origSurf;
	std::shared_ptr<Drawable> m_isoSurf_drawer;
	bool m_ready_to_draw;
	// eps used to measure closeness of two vts on MA/MC
	float m_eps_1;
	// eps used to measure closeness of two vts on iso-surface
	float m_eps_2;
	
	// the original surf vts extended with new vts on split edges
	vector<TriPoint> m_vts_origSurf;
	// the current iso-surf vts
	vector<TriPoint> m_vts_isoSurf;
	// edges of cur iso-surf (same for the extended original surf)
	vector<TriEdge> m_edges;
	// faces of cur iso-surf (same for the extended original surf)
	vector<TriFace> m_faces;
	map<TriEdge, vector<int>> m_nbFaces_for_edge;
	KNNTree* m_tree;

	// the corresp. orig vert -> MC/MA vert
	vector<int> m_closest_vert_for_orig;
	// the range of iso values (specific to the field type, e.g. bt3)
	float m_iso_max;
	// cur iso value
	float m_iso_value;
	// pre iso ratio param
	float m_pre_iso_ratio;
};

inline int
	IsoSurfFrom2Manifold::get_vert_type(int& _u)
{
	const auto& MA_vts = m_stg_ptr->m_stSubdiv.getAllVts();
// 	if (_u >= MA_vts.size())
// 	{
// 		_u -= MA_vts.size();
// 		if (m_stg_ptr->is_face_dual[_u])
// 			return 0;
// 		else
// 			return 1;
// 	}
// 	else if (m_stg_ptr->m_stSubdiv.isSteinerVert(_u))
// 		return 2;
// 	else
// 		return 3;

	return _u >= MA_vts.size() ? ( // a dual vert
		( _u -= MA_vts.size(), m_stg_ptr->m_is_face_dual[_u] >= 0 ) ? 
		0 : 1
		)
		: ( // original or steiner vert
		m_stg_ptr->m_stSubdiv.isSteinerVert(_u) ? 
		2 : 3
		);
}

inline int
	IsoSurfFrom2Manifold::encode_type_for_vts(int& _u, int& _v)
{
	int t_u = get_vert_type(_u);
	int t_v = get_vert_type(_v);
	if (t_u < t_v)
	{
		std::swap(t_u, t_v);
		std::swap(_u, _v);
	}
	return ( (t_u << 4) | t_v );
}

inline bool 
	IsoSurfFrom2Manifold::same_face_orig_orig(int _u, int _v)
{
	auto ret_code = m_stg_ptr->m_origG->getEdgeIdx(util::makeEdge(_u, _v));
	return ret_code != -1;
}
inline bool IsoSurfFrom2Manifold::same_face_orig_steiner(int _u, int _v)
{
	auto reside_e = m_stg_ptr->m_stSubdiv.getResidingEdge(_v);
	return on_same_edge_or_face(_u, reside_e[0], reside_e[1]);
}
inline bool IsoSurfFrom2Manifold::same_face_orig_edual(int _u, int _v)
{
	int ei = m_stg_ptr->m_dualV_triEdgeIdx[_v];
	auto reside_e = 
		m_stg_ptr->m_stSubdiv.getOrigTriEdge(
		ei );
	return on_same_edge_or_face(_u, reside_e[0], reside_e[1]);
}
inline bool IsoSurfFrom2Manifold::same_face_orig_fdual(int _u, int _v)
{
	auto reside_f = 
		m_stg_ptr->m_origG->faces[ m_stg_ptr->m_is_face_dual[_v] ];
	return is_face_vertex(_u, reside_f);
}
inline bool IsoSurfFrom2Manifold::same_face_steiner_steiner(int _u, int _v)
{
	return same_edge_or_on_same_face(
		m_stg_ptr->m_stSubdiv.getResidingEdge(_u),
		m_stg_ptr->m_stSubdiv.getResidingEdge(_v)
		);
}
inline bool IsoSurfFrom2Manifold::same_face_steiner_edual(int _u, int _v)
{
	return same_edge_or_on_same_face(
		m_stg_ptr->m_stSubdiv.getResidingEdge(_u),
		m_stg_ptr->m_stSubdiv.getOrigTriEdge( m_stg_ptr->m_dualV_triEdgeIdx[_v] )
		);
}
inline bool IsoSurfFrom2Manifold::same_face_steiner_fdual(int _u, int _v)
{
	auto reside_e = 
		m_stg_ptr->m_stSubdiv.getResidingEdge(_u);
	auto reside_f = 
		m_stg_ptr->m_origG->faces[ m_stg_ptr->m_is_face_dual[_v] ];
	return is_face_edge(reside_e, reside_f);
}
inline bool IsoSurfFrom2Manifold::same_face_edual_edual(int _u, int _v)
{
	return same_edge_or_on_same_face(
		m_stg_ptr->m_stSubdiv.getOrigTriEdge( m_stg_ptr->m_dualV_triEdgeIdx[_u] ),
		m_stg_ptr->m_stSubdiv.getOrigTriEdge( m_stg_ptr->m_dualV_triEdgeIdx[_v] )
		);
}
inline bool IsoSurfFrom2Manifold::same_face_edual_fdual(int _u, int _v)
{
	auto reside_e = 
		m_stg_ptr->m_stSubdiv.getOrigTriEdge( m_stg_ptr->m_dualV_triEdgeIdx[_u] );
	auto reside_f = 
		m_stg_ptr->m_origG->faces[ m_stg_ptr->m_is_face_dual[_v] ];
	return is_face_edge(reside_e, reside_f);
}
inline bool IsoSurfFrom2Manifold::same_face_fdual_fdual(int _u, int _v)
{
	return m_stg_ptr->m_is_face_dual[_u] == 
		m_stg_ptr->m_is_face_dual[_v] ;
}

/// low-level helpers
inline bool IsoSurfFrom2Manifold::on_same_edge_or_face(int _u, int _eu, int _ev)
{
	return (_u == _ev || _u == _ev) || 
		(m_stg_ptr->m_origG->getFaceIdx(util::makeFace(_u, _eu, _ev)) != -1
		);
}

inline bool IsoSurfFrom2Manifold::same_edge_or_on_same_face(const TriEdge& _e1, const TriEdge& _e2)
{
	return _e1 == _e2 || _e1[0] == _e2[0] || _e1[0] == _e2[1];
}

inline bool IsoSurfFrom2Manifold::is_face_vertex(int _u, const TriFace& _f)
{
	return _f[0] == _u || _f[1] == _u || _f[2] == _u;
}

inline bool IsoSurfFrom2Manifold::is_face_edge(const TriEdge& _e, const TriFace& _f)
{
	auto edges = util::edgesFromFace(_f);
	return _e == edges[0] || _e == edges[1] || _e == edges[2];
}

#endif