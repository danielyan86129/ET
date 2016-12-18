#include "IsoContourOnMA.h"
#include "drawable.h"
#include "meshDrawer.h"
#include "lineDrawer.h"
#include "geometryRepresentation.h"
#include <cmath>
/*
**TODO:
**The implementation of a few methods here relies on a data structure inside steiner graph
**  called m_finerTri_byDual, and some of its auxiliary info. See .h file of SteinerGraph for details.
**  That's not a right decision! 
**  All we want is a meaningful finer triangulation of the original MA for smoother display.
**  Presumably, we can just do well with the pre-existing triangulation, SteinerGraph::m_MA_finer_tris
**  to achieve the same goal. That structure is much simpler to use than m_finerTri_byDual
**  (only composed of original + steiner vts, not mixed with dual vts).
**  m_finerTri_byDual is designed for something different (i.e. used as 2-complex hybrid skeleton).
**Therefore, get rid of the dependence on m_finerTri_byDual, and use m_MA_finer_tris instead.
*/
void IsoContourOnMA::setup(const shared_ptr<SteinerGraph>& _stg, 
	bool _is_transparent, float _min_alpha, float _exp_alpha, 
	const float* const _new_min, const float* const _new_max)
{
	m_stg = _stg;
	
	// initially all faces ahead of iso-curve
	m_face_status_static.assign(m_stg->m_finerTri_byDual.size(), AHEAD);
	m_face_status_dynamic.clear();

	m_iso_max = 0.0f; m_iso_min = 0.0f;
	//m_iso_max = *( std::max_element(m_stg->bt2MA_face.begin(), m_stg->bt2MA_face.end()) );
	m_iso_min = 0.0f;
	for (int fi = 0; fi < m_stg->m_finerTri_byDual.size(); ++fi)
	{
		m_iso_max = std::max(m_iso_max, 
			m_stg->getScalarOfFineTri(SteinerGraph::BT2_FINE_TRI, fi) 
			);
	}

	if ( _new_min != nullptr && _new_max != nullptr ) // use new scale?
	{
		m_vis_min = *_new_min;
		m_vis_max = *_new_max;
	}
	else // range of iso (bt2 over ma face)
	{
		m_vis_min = m_iso_min;
		m_vis_max = m_iso_max;
	}
	
	m_fineMA_staticBehindIso_drawer->setTransparencyEnabled(true);
	m_fineMA_staticAheadIso_drawer->setTransparencyEnabled(false);
	m_fineMA_dynamicBehindIso_drawer->setTransparencyEnabled(true);
	m_fineMA_dynamicAheadIso_drawer->setTransparencyEnabled(false);
	this->setTransparencyEnabled(_is_transparent);
	m_min_alpha = _min_alpha; m_alpha_exp = _exp_alpha;
	this->upload_static_data();
	cout << "IsoContourOnMA::precompute() done." << endl;
}

void IsoContourOnMA::setupDrawers( 
	int _w, int _h, 
	shared_ptr<oglplus::Program> _simpleProg, 
	shared_ptr<oglplus::Program> _edgeProg, 
	shared_ptr<oglplus::Program> _linesProg, 
	shared_ptr<TrackBall> _trackBall )
{
	m_fineMA_staticBehindIso_drawer = std::shared_ptr<MeshDrawer>(
		new MeshDrawer(_w, _h, _simpleProg, _edgeProg, _trackBall)
		);
	m_fineMA_staticAheadIso_drawer = std::shared_ptr<MeshDrawer>(
		new MeshDrawer(_w, _h, _simpleProg, _edgeProg, _trackBall)
		);
	m_fineMA_dynamicBehindIso_drawer = std::shared_ptr<MeshDrawer>(
		new MeshDrawer(_w, _h, _simpleProg, _edgeProg, _trackBall)
		);
	m_fineMA_dynamicAheadIso_drawer = std::shared_ptr<MeshDrawer>(
		new MeshDrawer(_w, _h, _simpleProg, _edgeProg, _trackBall)
		);
	m_cont_line_drawer = std::shared_ptr<LineDrawer>(
		new LineDrawer(_linesProg, _trackBall)
		);

}

void IsoContourOnMA::setupGeometryScale(float _s)
{
	std::dynamic_pointer_cast<MeshDrawer>(m_fineMA_staticBehindIso_drawer)->setScale(_s);
	std::dynamic_pointer_cast<MeshDrawer>(m_fineMA_staticAheadIso_drawer)->setScale(_s);
	std::dynamic_pointer_cast<MeshDrawer>(m_fineMA_dynamicBehindIso_drawer)->setScale(_s);
	std::dynamic_pointer_cast<MeshDrawer>(m_fineMA_dynamicAheadIso_drawer)->setScale(_s);
	std::dynamic_pointer_cast<LineDrawer>(m_cont_line_drawer)->setScale(_s);
}

void IsoContourOnMA::reset()
{
	m_stg.reset();
	/*m_cont_line_drawer.reset();
	m_fineMA_static_drawer.reset();
	m_fineMA_dynamic_drawer.reset();*/
}

void IsoContourOnMA::precompute()
{
	// initially all faces ahead of iso-curve
	m_face_status_static.assign(m_stg->m_finerTri_byDual.size(), AHEAD);
	m_face_status_dynamic.clear();

	// range of iso (bt2 over ma face)
	m_iso_max = 0.0f; m_iso_min = 0.0f;
	//m_iso_max = *( std::max_element(m_stg->bt2MA_face.begin(), m_stg->bt2MA_face.end()) );
 	m_iso_min = 0.0f;
 	for (int fi = 0; fi < m_stg->m_finerTri_byDual.size(); ++fi)
 	{
 		m_iso_max = std::max(m_iso_max, 
			m_stg->getScalarOfFineTri(SteinerGraph::BT2_FINE_TRI, fi) 
			);
 	}

	m_fineMA_staticBehindIso_drawer->setTransparencyEnabled(true);
	m_fineMA_staticAheadIso_drawer->setTransparencyEnabled(false);
	m_fineMA_dynamicBehindIso_drawer->setTransparencyEnabled(true);
	m_fineMA_dynamicAheadIso_drawer->setTransparencyEnabled(false);
	this->upload_static_data();
	cout << "IsoContourOnMA::precompute() done." << endl;
}

void IsoContourOnMA::genContour(float _iso_ratio)
{
	float iso = getIsoValue(_iso_ratio);
	genContour(iso, 1);
}

void IsoContourOnMA::genContour(float _iso, int)
{
	this->update_contour(_iso);
	this->upload_cont_data();
}

void IsoContourOnMA::render(bool _show_static_MA)
{
	m_fineMA_staticBehindIso_drawer->render(0.0);	
	m_fineMA_staticAheadIso_drawer->render(0.0);	
	m_cont_line_drawer->render(0.0);
	m_fineMA_dynamicAheadIso_drawer->render(0.0);
	m_fineMA_dynamicBehindIso_drawer->render(0.0);
	//cout << "IsoContourOnMA::render() !" << endl;
}

void IsoContourOnMA::submitDrawCalls(
	vector<std::pair<std::shared_ptr<Drawable>,bool>>& _transparent_drawCalls, 
	vector<std::pair<std::shared_ptr<Drawable>,bool>>& _opaque_drawCalls
	)
{
	pair<std::shared_ptr<Drawable>, bool> drawers[5] = {
		std::make_pair(m_fineMA_staticBehindIso_drawer, true),
		std::make_pair(m_fineMA_staticAheadIso_drawer, true),
		std::make_pair(m_fineMA_dynamicBehindIso_drawer, true),
		std::make_pair(m_fineMA_dynamicAheadIso_drawer, true),
		std::make_pair(m_cont_line_drawer, false)
	};
	for (int i = 0; i < 5; ++i)
	{
		const auto& draw_pr = drawers[i];
		draw_pr.first->getTransparencyEnabled() ? 
			_transparent_drawCalls.push_back(draw_pr) 
			: _opaque_drawCalls.push_back(draw_pr);
	}
}

void IsoContourOnMA::reshape(int _w, int _h)
{
	m_fineMA_staticBehindIso_drawer->reshape(_w, _h);	
	m_fineMA_staticAheadIso_drawer->reshape(_w, _h);	
	m_cont_line_drawer->reshape(_w, _h);
	m_fineMA_dynamicBehindIso_drawer->reshape(_w, _h);
	m_fineMA_dynamicAheadIso_drawer->reshape(_w, _h);
}

void IsoContourOnMA::enableWireFrame(bool _enable)
{
	std::dynamic_pointer_cast<MeshDrawer>(m_fineMA_staticBehindIso_drawer)->setDrawEdge(_enable);
	std::dynamic_pointer_cast<MeshDrawer>(m_fineMA_staticAheadIso_drawer)->setDrawEdge(_enable);
	std::dynamic_pointer_cast<MeshDrawer>(m_fineMA_dynamicBehindIso_drawer)->setDrawEdge(_enable);
	std::dynamic_pointer_cast<MeshDrawer>(m_fineMA_dynamicAheadIso_drawer)->setDrawEdge(_enable);
}

void IsoContourOnMA::setTransparencyEnabled(bool _transparent)
{
	std::dynamic_pointer_cast<MeshDrawer>(
		m_fineMA_dynamicBehindIso_drawer
		)->setTransparencyEnabled(_transparent);

	std::dynamic_pointer_cast<MeshDrawer>(
		m_fineMA_staticBehindIso_drawer
		)->setTransparencyEnabled(_transparent);
}

void IsoContourOnMA::setZFightOffset(float _r)
{
	std::dynamic_pointer_cast<LineDrawer>(
		m_cont_line_drawer
		)->setZFightOffset(_r);
}

void IsoContourOnMA::setLightingEnabled(bool _enable_lighting)
{
	std::dynamic_pointer_cast<MeshDrawer>(
		m_fineMA_staticBehindIso_drawer
		)->setLightingEnabled(_enable_lighting);
	std::dynamic_pointer_cast<MeshDrawer>(
		m_fineMA_dynamicBehindIso_drawer
		)->setLightingEnabled(_enable_lighting);
}

void IsoContourOnMA::setTransparencyParam(float _min, float _exp)
{
	this->m_min_alpha = _min; 
	this->m_alpha_exp = _exp;
}

///private helpers

void IsoContourOnMA::update_contour(float _iso)
{
	m_vts_cont.clear();
	m_edges_cont.clear();
	m_faces_ahead_iso.clear();
	std::fill(m_face_status_static.begin(), m_face_status_static.end(), AHEAD);
	m_vts_dynamic.clear();
	m_scalarField_dynamic.clear();
	m_faces_dynamic.clear();
	m_dynamicFace_fromFineFace.clear();
	m_face_status_dynamic.clear();

	// for each face, test if it contains part of iso-curve
	// if so, extract that part via triangle marching 
	// and tessellate the face properly using the iso-curve line
	const auto& faces = m_stg->m_finerTri_byDual;
	TriFace f; // cur face
	int from_tri_fi; // which original tri face of MA is f from?
	int f_status;
	TriEdge e; // edge of cur face
	vector<bool> e_has_iso(3, false); // does any edge contains iso-curve bit?
	vector<TriEdge> edges_that_cross; // the edges that cross the iso-curve. At most 2.
	vector<TriPoint> iso_points(3); // the iso point on each edge
	float value_at_vert[3]; // the distance value at three vertices
	TriPoint u, v;

	int cross_cnt_all = 0; // debug
	int past_cnt_all = 0; // debug
	for (unsigned fi = 0; fi < faces.size(); ++fi)
	{
		edges_that_cross.clear();

		from_tri_fi = m_stg->m_origTri_for_finerTri_byDual[fi];
		// skip this face if it's pruned
		if (m_stg->isFacePruned(from_tri_fi))
		{
			m_face_status_static[fi] = BEHIND;
			continue;
		}

		f = faces[fi];
		f_status = iso_curve_crosses(f, from_tri_fi, _iso, value_at_vert);
		if (f_status != CROSSED) // not intersected by iso-curve
		{
			if (f_status == AHEAD) // face ahead of iso-curve
			{
				// add face to to-display list
				m_faces_ahead_iso.push_back(fi);
				m_face_status_static[fi] = AHEAD;
			}
			else
			{
				past_cnt_all ++;
				m_face_status_static[fi] = BEHIND;
			}
			continue;
		}
		m_face_status_static[fi] = CROSSED;
		cross_cnt_all ++;

		for (unsigned i = 0; i < 3; ++i)
		{
			//e_has_iso[ei] = iso_curve_crosses(e, from_tri_fi, _iso, iso_points[ei]);
			e_has_iso[i] = iso_curve_crosses(
				f[i], f[(i+1)%3], 
				value_at_vert[i], value_at_vert[(i+1)%3], 
				_iso, iso_points[i]);
			if (e_has_iso[i])
				m_vts_cont.push_back(iso_points[i]);
		}

		// it is not possible that all three edges cross the iso-curve
		// or only one edge crosses the iso-curve
		int cross_cnt = e_has_iso[0] + e_has_iso[1] + e_has_iso[2] ;
		assert (cross_cnt != 3/* && cross_cnt != 1*/);

		if (cross_cnt == 0)
			continue;
		else if (cross_cnt == 1)
			m_edges_cont.push_back(util::makeEdge(m_vts_cont.size()-1, m_vts_cont.size()-1));
		else 
		{
			//2 edges intersected iso-curve. 
			// extract the iso-curve line
			m_edges_cont.push_back(util::makeEdge(m_vts_cont.size()-1, m_vts_cont.size()-2));

			// tessellate the face properly to retain the part that's ahead of iso-curve
			TriEdge e1, e2;
			int i, j;
			int shared_v;
			if (e_has_iso[0] && e_has_iso[1]) // 0 & 1th edges of face 
			{
				e1 = TriEdge(f[0], f[1]);
				e2 = TriEdge(f[1], f[2]);
				i = 0; j = 1;
				shared_v = 1;
			}
			else if (e_has_iso[1] && e_has_iso[2]) // 1 & 2th edges of face 
			{
				e1 = TriEdge(f[1], f[2]);
				e2 = TriEdge(f[2], f[0]);
				i = 1; j = 2;
				shared_v = 2;
			}
			else // 0 & 2th edges of face 
			{
				e1 = TriEdge(f[0], f[1]);
				e2 = TriEdge(f[0], f[2]);
				i = 0; j = 2;
				shared_v = 0;
			}

			// put all dynamic faces in the list
			int e_other_end = util::otherEnd(e1, f[shared_v]);
			int f_other_end = util::otherEnd(e2, f[shared_v]);
			float scalar_e_other_end = m_stg->getScalarOfFineTriVert(
				SteinerGraph::BT2_FINE_TRI_VERT, e_other_end, fi);
			float scalar_f_other_end = m_stg->getScalarOfFineTriVert(
				SteinerGraph::BT2_FINE_TRI_VERT, f_other_end, fi);

			// a dynamic tri. formed by shared v and two iso vts.
			m_vts_dynamic.push_back( m_stg->getVertPosOfFineTri(f[shared_v]) ); 
			m_vts_dynamic.push_back(iso_points[i]);
			m_vts_dynamic.push_back(iso_points[j]);
			m_scalarField_dynamic.push_back(m_stg->getScalarOfFineTriVert(
				SteinerGraph::BT2_FINE_TRI_VERT, f[shared_v], fi)
				);
			m_scalarField_dynamic.push_back(_iso);
			m_scalarField_dynamic.push_back(_iso);
			m_faces_dynamic.push_back(util::makeFace(
				m_vts_dynamic.size()-3, m_vts_dynamic.size()-2, m_vts_dynamic.size()-1));
			m_dynamicFace_fromFineFace.push_back(fi);

			// a dynamic quad formed by two ends of the edge and two iso vts.
			m_vts_dynamic.push_back(iso_points[i]);
			m_vts_dynamic.push_back( m_stg->getVertPosOfFineTri(e_other_end) );
			m_vts_dynamic.push_back(iso_points[j]);
			m_scalarField_dynamic.push_back(_iso);
			m_scalarField_dynamic.push_back( scalar_e_other_end );
			m_scalarField_dynamic.push_back(_iso);
			m_faces_dynamic.push_back(util::makeFace(
				m_vts_dynamic.size()-3, m_vts_dynamic.size()-2, m_vts_dynamic.size()-1));
			m_vts_dynamic.push_back( m_stg->getVertPosOfFineTri(e_other_end) );
			m_vts_dynamic.push_back( m_stg->getVertPosOfFineTri(f_other_end) );
			m_vts_dynamic.push_back(iso_points[j]);
			m_scalarField_dynamic.push_back( scalar_e_other_end );
			m_scalarField_dynamic.push_back( scalar_f_other_end );
			m_scalarField_dynamic.push_back(_iso);
			m_faces_dynamic.push_back(util::makeFace(
				m_vts_dynamic.size()-3, m_vts_dynamic.size()-2, m_vts_dynamic.size()-1));
			m_dynamicFace_fromFineFace.push_back(fi);
			m_dynamicFace_fromFineFace.push_back(fi);

			if (value_at_vert[shared_v] > _iso) // shared vert's value ahead of iso. 
			{
				// the tri face ahead, the quad face behind.
				m_face_status_dynamic.push_back(AHEAD);
				m_face_status_dynamic.push_back(BEHIND);
				m_face_status_dynamic.push_back(BEHIND);
			}
			else // it's a quad to retain
			{
				// the quad face ahead, the tri face behind.
				m_face_status_dynamic.push_back(BEHIND);
				m_face_status_dynamic.push_back(AHEAD);
				m_face_status_dynamic.push_back(AHEAD);
			}

			/*
			// only put the dynamic face ahead of iso curve in the list
			if (value_at_vert[shared_v] > _iso) // shared vert's value ahead of iso. it's a triangle to retain
			{
				m_vts_dynamic.push_back( m_stg->getVertPosOfFineTri(f[shared_v]) ); 
				m_vts_dynamic.push_back(iso_points[i]);
				m_vts_dynamic.push_back(iso_points[j]);
				m_faces_dynamic.push_back(util::makeFace(
					m_vts_dynamic.size()-3, m_vts_dynamic.size()-2, m_vts_dynamic.size()-1));
				m_dynamicFace_fromFace.push_back(fi);
			}
			else // it's a quad to retain
			{
				int e_other_end = util::otherEnd(e1, f[shared_v]);
				int f_other_end = util::otherEnd(e2, f[shared_v]);
				m_vts_dynamic.push_back(iso_points[i]);
				m_vts_dynamic.push_back( m_stg->getVertPosOfFineTri(e_other_end) );
				m_vts_dynamic.push_back(iso_points[j]);
				m_faces_dynamic.push_back(util::makeFace(
					m_vts_dynamic.size()-3, m_vts_dynamic.size()-2, m_vts_dynamic.size()-1));
				m_vts_dynamic.push_back( m_stg->getVertPosOfFineTri(e_other_end) );
				m_vts_dynamic.push_back( m_stg->getVertPosOfFineTri(f_other_end) );
				m_vts_dynamic.push_back(iso_points[j]);
				m_faces_dynamic.push_back(util::makeFace(
					m_vts_dynamic.size()-3, m_vts_dynamic.size()-2, m_vts_dynamic.size()-1));
				m_dynamicFace_fromFace.push_back(fi);
				m_dynamicFace_fromFace.push_back(fi);
			} */
		}
	}

	//cout << "iso-cont vts/edges: " << m_vts_cont.size() << "/"<<m_edges_cont.size()<<endl; //debug
	/*cout << "# faces all/ahead/crossed/past: " 
		<< faces.size() <<"/"
		<< m_faces_ahead_iso.size() <<"/"
		<< cross_cnt_all << "/"
		<< past_cnt_all
		<< endl;*/
	assert( m_vts_dynamic.size() == 3 * m_faces_dynamic.size() );
	//cout << "ma_dynamic vts/faces: " << m_vts_dynamic.size() << "/"<<m_faces_dynamic.size()<<endl; //debug
}

void IsoContourOnMA::upload_cont_data()
{
	// upload contour data. use per vert color.
	auto& line_drawer = std::dynamic_pointer_cast<LineDrawer>(m_cont_line_drawer);
	float* color_data = new float[3 * m_vts_cont.size()];
	TriColor contour_color(1.0f, 0.0f, 0.0f);
	for (unsigned i = 0; i < m_vts_cont.size(); ++i)
	{
		color_data[3*i+0] = contour_color[0];
		color_data[3*i+1] = contour_color[1];
		color_data[3*i+2] = contour_color[2];
	}
	line_drawer->setPoints(m_vts_cont);
	line_drawer->setLines(m_edges_cont);
	line_drawer->setPerVertColor(color_data, m_vts_cont.size());
	delete [] color_data;

	//const auto& ma_face_bt2 = m_stg->bt2MA_face;
	/// BEGIN: upload the DYNAMIC faces data
	dynamic_behindIsoFaces.clear();
	dynamic_aheadIsoFaces.clear();
	auto& dynamicBehindIso_drawer = std::dynamic_pointer_cast<MeshDrawer>(m_fineMA_dynamicBehindIso_drawer);
	auto& dynamicAheadIso_drawer = std::dynamic_pointer_cast<MeshDrawer>(m_fineMA_dynamicAheadIso_drawer);
	dynamicBehindIso_drawer->setRenderMode(MeshDrawer::PER_VERT);
	dynamicBehindIso_drawer->setPoints(m_vts_dynamic);
	dynamicAheadIso_drawer->setRenderMode(MeshDrawer::PER_VERT);
	dynamicAheadIso_drawer->setPoints(m_vts_dynamic);

	// since shared vertices are duplicated, it's OK to 
	// collect colors for BEHIND and AHEAD faces' vertices into one array.
	// Correct IBO for each drawer (BEHIND and AHEAD, resp) will fetch colors from correctly indexed positions.
	color_data = new float[m_vts_dynamic.size() * 3];
	float* saliency_data = new float[m_vts_dynamic.size()];
	TriColor constant_color(0.4f, 0.4f, 0.4f);
	TriColor vert_color;
	float saliency;
	float* nml_data = new float[m_vts_dynamic.size() * 3];
	trimesh::vec3 nml;
	for (unsigned fi = 0; fi < m_faces_dynamic.size(); ++fi)
	{
		const auto& f = m_faces_dynamic[fi];
		float scalar;
		// color data
		if (m_face_status_dynamic[fi] == BEHIND)
		{
			dynamic_behindIsoFaces.push_back(f);
			for (unsigned j = 0; j < 3; ++j)
			{
				scalar = m_scalarField_dynamic[f[j]];
				vert_color = util::GetColour(
					scalar, m_vis_min, m_vis_max
					);

				color_data[3*f[j]+0] = vert_color[0];
				color_data[3*f[j]+1] = vert_color[1];
				color_data[3*f[j]+2] = vert_color[2];

				saliency_data[f[j]] = util::rescale(
					scalar, m_vis_min, m_vis_max, m_alpha_exp, m_min_alpha, 1.0f
					);
			}
		}
		else
		{
			dynamic_aheadIsoFaces.push_back(f);
			saliency = 1.0;
			for (unsigned j = 0; j < 3; ++j)
			{
				color_data[3*f[j]+0] = constant_color[0];
				color_data[3*f[j]+1] = constant_color[1];
				color_data[3*f[j]+2] = constant_color[2];
			}
		}
		
		// normal data
		nml = trimesh::normalize(
			trimesh::trinorm(
			m_vts_dynamic[f[0]], m_vts_dynamic[f[1]], m_vts_dynamic[f[2]])
			);
		for (unsigned j = 0; j < 3; ++j)
		{
			nml_data[3*f[j]+0] = nml[0];
			nml_data[3*f[j]+1] = nml[1];
			nml_data[3*f[j]+2] = nml[2];
		}
	}
	dynamicBehindIso_drawer->setFaces(dynamic_behindIsoFaces, false);
	dynamicBehindIso_drawer->setPerVertColor(color_data, m_vts_dynamic.size());
	dynamicBehindIso_drawer->setPerVertSaliency(saliency_data, m_vts_dynamic.size());
	dynamicBehindIso_drawer->setPerVertNormal(nml_data, m_vts_dynamic.size());
	
	dynamicAheadIso_drawer->setFaces(dynamic_aheadIsoFaces, false);
	dynamicAheadIso_drawer->setPerVertColor(color_data, m_vts_dynamic.size());
	dynamicAheadIso_drawer->setPerVertSaliency(saliency_data, m_vts_dynamic.size());
	dynamicAheadIso_drawer->setPerVertNormal(nml_data, m_vts_dynamic.size());

	delete [] nml_data;
	delete [] color_data;
	delete [] saliency_data;
	/// END: upload the DYNAMIC faces data

	/// BEGIN: specify STATIC faces that are BEHIND and AHEAD into the corresponding drawer, resp.
	const auto& fine_faces = m_stg->m_finerTri_byDual;
	const auto& fine_tri_fromFace = m_stg->m_origTri_for_finerTri_byDual;
	
	vector<unsigned> static_behindIso_faces_indices_to_draw;
	vector<unsigned> static_aheadIso_faces_indices_to_draw;
	for (int fi = 0; fi < fine_faces.size(); ++fi)
	{
		if (m_stg->isFacePruned(fine_tri_fromFace[fi]))
			continue;

		if (m_face_status_static[fi] == BEHIND)
		{
			static_behindIso_faces_indices_to_draw.push_back(fi);
		}
		else if (m_face_status_static[fi] == AHEAD)
		{
			static_aheadIso_faces_indices_to_draw.push_back(fi);
		}
	}

	auto& staticBehindIsoFaces_drawer = std::dynamic_pointer_cast<MeshDrawer>(m_fineMA_staticBehindIso_drawer);
	staticBehindIsoFaces_drawer->setFaces(static_behindIso_faces_indices_to_draw, int(0));
	auto& staticAheadIsoFaces_drawer = std::dynamic_pointer_cast<MeshDrawer>(m_fineMA_staticAheadIso_drawer);
	staticAheadIsoFaces_drawer->setFaces(static_aheadIso_faces_indices_to_draw, int(0));
	/// BEGIN: collect STATIC faces
}

void IsoContourOnMA::upload_static_data()
{
	auto& static_drawer_behindIso = std::dynamic_pointer_cast<MeshDrawer>(m_fineMA_staticBehindIso_drawer);
	static_drawer_behindIso->setRenderMode(MeshDrawer::PER_VERT);
	auto& static_drawer_aheadIso = std::dynamic_pointer_cast<MeshDrawer>(m_fineMA_staticAheadIso_drawer);
	static_drawer_aheadIso->setRenderMode(MeshDrawer::PER_VERT);

	// upload vts & faces of the finner tessellation of MA
	vector<TriPoint> vts;
	vts.reserve(m_stg->m_finerTri_byDual.size() * 3);
	vector<float> scalar_field;
	scalar_field.reserve(m_stg->m_finerTri_byDual.size() * 3);
	auto faces = m_stg->m_finerTri_byDual;
	for (int fi = 0; fi < faces.size(); ++fi)
	{
		const auto& f = faces[fi];
		// duplicate shared vertices
		for (int j = 0; j < 3; ++j)
		{
			int vi = f[j];
			TriPoint p = m_stg->getVertPosOfFineTri(vi);
			vts.push_back(p);
			scalar_field.push_back(
				m_stg->getScalarOfFineTriVert(SteinerGraph::BT2_FINE_TRI_VERT, vi, fi)
				);
		}
		faces[fi] = util::makeFace(vts.size()-1, vts.size()-2, vts.size()-3);
	}
	static_drawer_behindIso->setPoints(vts);
	static_drawer_behindIso->setFaces(faces, false);
	static_drawer_aheadIso->setPoints(vts);
	static_drawer_aheadIso->setFaces(faces, false);
	
	// now set colors and normals info into drawer for static faces AHEAD of iso contour
	float scalar;
	TriColor vert_color;
	TriColor constant_color(0.4f, 0.4f, 0.4f);
	trimesh::vec3 face_normal;
	float* color_data = new float[ vts.size()*3 ];
	float* normal_data = new float[ vts.size()*3 ];
	for (int fi = 0; fi < faces.size(); ++fi)
	{
		const auto& f = faces[fi];
		face_normal = trimesh::normalize(trimesh::trinorm(vts[f[0]], vts[f[1]], vts[f[2]]));
		// duplicate shared vertices
		for (int j = 0; j < 3; ++j)
		{
			int vi = f[j];
			// per vert color constant for ahead faces 
			color_data[vi * 3 + 0] = constant_color[0];
			color_data[vi * 3 + 1] = constant_color[1];
			color_data[vi * 3 + 2] = constant_color[2];

			normal_data[3*vi + 0] = face_normal[0];
			normal_data[3*vi + 1] = face_normal[1];
			normal_data[3*vi + 2] = face_normal[2];
		}
	}
	static_drawer_aheadIso->setPerVertColor(color_data, vts.size());
	static_drawer_aheadIso->setPerVertNormal(normal_data, vts.size());

	// now set colors and normals info into drawer for static faces BEHIND the iso contour
	float* saliency_data = new float[ vts.size() ];
	/*float* color_data = new float[ faces.size()*3 ];
	float* normal_data = new float[ faces.size()*3 ];*/
	for (int fi = 0; fi < faces.size(); ++fi)
	{
		const auto& f = faces[fi];
		face_normal = trimesh::normalize(trimesh::trinorm(vts[f[0]], vts[f[1]], vts[f[2]]));
		// duplicate shared vertices
		for (int j = 0; j < 3; ++j)
		{
			int vi = f[j];
			// per vert color according to properly selected scalar value 
			scalar = scalar_field[vi];
			vert_color = util::GetColour(scalar, m_vis_min, m_vis_max);
			color_data[vi * 3 + 0] = vert_color[0];
			color_data[vi * 3 + 1] = vert_color[1];
			color_data[vi * 3 + 2] = vert_color[2];

			// per vert normal (use face normal)
			normal_data[3*vi + 0] = face_normal[0];
			normal_data[3*vi + 1] = face_normal[1];
			normal_data[3*vi + 2] = face_normal[2];

			// importance of this face (passed to gpu as opaqueness)
			saliency_data[vi] = util::rescale(scalar, m_vis_min, m_vis_max, m_alpha_exp, m_min_alpha, 1.0f);
		}
	}
	static_drawer_behindIso->setPerVertColor(color_data, vts.size());
	static_drawer_behindIso->setPerVertNormal(normal_data, vts.size());
	static_drawer_behindIso->setPerVertSaliency(saliency_data, vts.size());

	delete [] color_data;
	delete [] normal_data;
	delete [] saliency_data;

	cout << "IsoContourOnMA::upload_static_data() - finner MA uploaded." << endl;
}