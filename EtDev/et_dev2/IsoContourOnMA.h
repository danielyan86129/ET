// Generate & display iso-contour on the finer triangulation of MA
// 
// Copyright (C) 2018 Yajie Yan <danielyan86129@hotmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef ISO_CONT_ON_MA
#define ISO_CONT_ON_MA

#include "all_gl.h"

#include "allTypes.h"
#include "commonUtils.h"
#include "geometryRepresentation.h"

//////////////////////////////////////////////////////////////////////////
/// generate & display iso-contour on the finer triangulation of MA
//////////////////////////////////////////////////////////////////////////
class IsoContourOnMA
{
public:
	IsoContourOnMA() {}
	~IsoContourOnMA() {}

	/* setup the iso-cont object by passing in the steiner graph representing the MA */
	void setup(
		const shared_ptr<SteinerGraph>& _stg, 
		bool _is_transparent, float _min_alpha, float _exp_alpha, 
		const float* const _new_min, const float* const _new_max);
	void setupDrawers(
		int _w, int _h, 
		shared_ptr<oglplus::Program> _simpleProg, 
		shared_ptr<oglplus::Program> _edgeProg,
		shared_ptr<oglplus::Program> _linesProg,
		shared_ptr<TrackBall> _trackBall
		);
	void setupGeometryScale(float _s);
	/* reset this object (release held shared_ptr, clean states, data, etc.) */
	void reset();
	/* does as much as possible here 
	before interactive contour generation and rendering */
	void precompute();
	/* return the iso value corresponding to the _ratio */
	inline float getIsoValue(float _iso_ratio)
	{
		return (m_iso_max-m_iso_min) * _iso_ratio + m_iso_min;
	}
	/* generate & upload contour geometry to the drawer */
	void genContour(float _iso_ratio);
	void genContour(float _iso, int);
	/* render the contour line on MA's finer tessellation */
	void render(bool _show_static_MA);
	void submitDrawCalls( 
		vector<std::pair<std::shared_ptr<Drawable>,bool>>& _transparent_drawCalls, 
		vector<std::pair<std::shared_ptr<Drawable>,bool>>& _opaque_drawCalls);
	/* reshape func (usually called by render loop) */
	void reshape(int _w, int _h);
	/* set the wire-frame state */ 
	void enableWireFrame(bool _enable);

	/* rendering related */
	// control transparency for the faces burned already
	void setTransparencyEnabled(bool _transparent);
	// set offset for lines to fix z-fighting
	void setZFightOffset(float _r);
	// control lighting for the faces burned already
	void setLightingEnabled(bool _enable_lighting);
	void setTransparencyParam(float _min, float _exp);

	/// public data members
public:
	// BEHIND: behind iso-curve; CROSSED: crosses iso-curve; AHEAD: ahead of iso-curve;
	enum FaceStatus {BEHIND, CROSSED, AHEAD};
	vector<int> m_face_status_static;

	/// private helpers
private:
	/* fill in data member with updated iso-contour geometry */
	void update_contour(float _iso);
	/* upload static data to drawer that doesn't change from frame to frame */
	void upload_static_data();
	/* upload contour data to the drawer that may change frequently */
	void upload_cont_data();
	/* Does iso-curve crosses the given edge? If so, return the iso-point on the edge. */
	inline bool iso_curve_crosses(TriEdge _fine_e, int _MA_fi, float _iso, TriPoint& _p);
	inline bool iso_curve_crosses(int _p1, int _p2, float _v1, float _v2, float _iso, TriPoint& _p);
	/* position relative to iso-curve. -1: behind, 0: intersected, 1: ahead */
	inline int iso_curve_crosses(TriFace _f, int _MA_fi, float _iso, float* _value_at_vert);

	/// private data members
private:
	shared_ptr<SteinerGraph> m_stg;
	/// drawer for contour line geometry
	shared_ptr<Drawable> m_cont_line_drawer;
	/// drawer for static & dynamic part of MA's finer tessellation
	shared_ptr<Drawable> m_fineMA_staticBehindIso_drawer;
	shared_ptr<Drawable> m_fineMA_staticAheadIso_drawer;
	shared_ptr<Drawable> m_fineMA_dynamicBehindIso_drawer;
	shared_ptr<Drawable> m_fineMA_dynamicAheadIso_drawer;
	/// contour vts & edges
	vector<TriPoint> m_vts_cont;
	vector<TriEdge> m_edges_cont;
	/// the dynamic part: small faces left after proper subdivision of the faces containing the iso-contour
	vector<TriPoint> m_vts_dynamic;
	vector<TriFace> m_faces_dynamic;
	vector<float> m_scalarField_dynamic;
	vector<FaceStatus> m_face_status_dynamic;
	vector<int> m_dynamicFace_fromFineFace;
	vector<TriFace> dynamic_behindIsoFaces;
	vector<TriFace> dynamic_aheadIsoFaces;
	/// the faces that are ahead of iso-contour
	vector<unsigned> m_faces_ahead_iso;

	/// possibly min/max iso value (max bt2 here)
	float m_iso_min, m_iso_max;

	/// rendering related
	// transparency
	float m_min_alpha;
	float m_alpha_exp;
	// scalar range for vis.
	float m_vis_min, m_vis_max;
};

/// inline private helpers

inline bool IsoContourOnMA::iso_curve_crosses(TriEdge _fine_e, int _MA_fi, float _iso, TriPoint& _p)
{
// 	if (_fine_e == TriEdge(5, 128))
// 		int stop = 1; // debug

	TriPoint e0, e1;
	float bt2_e0, bt2_e1;
	int tei;
	const auto& st_vts = m_stg->m_stSubdiv.getAllVts();
	int u;
	int v;

	if ( m_stg->isDualVertInFineTri(_fine_e[0], u) ) // a face dual point
	{
		bt2_e0 = m_stg->bt2_MC[u];
		e0 = m_stg->dual_vts[u];
	}
	else // a steiner vert
	{
		tei = m_stg->mapTopo(u, _MA_fi);
		assert(tei >= 0);
		bt2_e0 = m_stg->bt2MA_vert_per_sheet[u][tei];
		e0 = st_vts[u];
	}

	if (m_stg->isDualVertInFineTri(_fine_e[1], v)) // a face dual point
	{
		bt2_e1 = m_stg->bt2_MC[v];
		e1 = m_stg->dual_vts[v];
	}
	else // a steiner vert
	{
		tei = m_stg->mapTopo(v, _MA_fi);
		assert(tei >= 0);
		bt2_e1 = m_stg->bt2MA_vert_per_sheet[v][tei];
		e1 = st_vts[v];
	}

	if ( (bt2_e0 - _iso) * (bt2_e1 - _iso) < 0.0f )
	{
		float l = trimesh::dist(e0, e1);
		if (::almost_zero(l, 0.000000001f))
			_p = e0;
		else
			_p = trimesh::mix(e0, e1, std::abs(bt2_e0 - _iso) / std::abs(bt2_e0-bt2_e1));
		return true;
	}
	else
		return false;
}

inline bool IsoContourOnMA::iso_curve_crosses(int _u, int _v, float _val1, float _val2, float _iso, TriPoint& _p)
{
	TriPoint e0, e1;
	const auto& st_vts = m_stg->m_stSubdiv.getAllVts();

	e0 = m_stg->getVertPosOfFineTri(_u);
	e1 = m_stg->getVertPosOfFineTri(_v);

	if ( (_val1 - _iso) * (_val2 - _iso) < 0.0f )
	{
		float l = std::abs(_val1 - _val2);
		if (::almost_zero(l, 0.000000001f))
			_p = e0;
		else
			_p = trimesh::mix( e0, e1, std::abs(_val1 - _iso) / l );
		return true;
	}
	else
		return false;
}

inline int IsoContourOnMA::iso_curve_crosses(TriFace _f, int _MA_fi, float _iso, float* _value_at_vert)
{
	float bt2_v;
	int tei;
	for (int i = 0; i < 3; ++i)
	{
		int vi = _f[i];
		int new_vi;
		if (m_stg->isDualVertInFineTri(vi, new_vi)/*vi < 0*/) // a face dual point
		{
			bt2_v = m_stg->bt2_MC[new_vi]; //m_stg->bt2_MC_perFace.find(trimesh::ivec2(-vi, _MA_fi))->second
		}
		else // a steiner vert
		{
			tei = m_stg->mapTopo(new_vi, _MA_fi);
			assert(tei >= 0);
			bt2_v = m_stg->bt2MA_vert_per_sheet[new_vi][tei];
		}
		_value_at_vert[i] = bt2_v/* - _iso*/;
	}

	if (_value_at_vert[0] >= _iso && _value_at_vert[1] >= _iso && _value_at_vert[2] >= _iso )
		return AHEAD;
	else if (_value_at_vert[0] < _iso && _value_at_vert[1] < _iso && _value_at_vert[2] < _iso )
		return BEHIND;
	else
		return CROSSED;
}

#endif