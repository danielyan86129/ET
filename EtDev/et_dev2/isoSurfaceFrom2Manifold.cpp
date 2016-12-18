#include "isoSurfaceFrom2Manifold.h"
#include "geometryRepresentation.h"
#include "meshDrawer.h"

#include <queue>

const int IsoSurfFrom2Manifold::vType_priority[] = {3, 2, 1, 0};
//IsoSurfFrom2Manifold::vTypePairHandler 
//	IsoSurfFrom2Manifold::vType_pair_entry[0x33] = 

IsoSurfFrom2Manifold::IsoSurfFrom2Manifold()
{
	m_ready_to_draw = false;
	m_tree = nullptr;
}

IsoSurfFrom2Manifold::~IsoSurfFrom2Manifold()
{
	if (m_tree)
		delete m_tree;
}

void IsoSurfFrom2Manifold::reset()
{
	m_stg_ptr.reset();
	m_origSurf.reset();
	//m_isoSurf_drawer.reset();
	m_ready_to_draw = false;
}

void IsoSurfFrom2Manifold::setup(
	const std::shared_ptr<SteinerGraph>& _stg_ptr, 
	const std::shared_ptr<TriMesh>& _3d_surf, 
	std::shared_ptr<Drawable>& _surf_drawer)
{
	m_stg_ptr = _stg_ptr;
	m_origSurf = _3d_surf;
	m_isoSurf_drawer = _surf_drawer;

	m_vts_origSurf = m_origSurf->vertices;
	m_vts_isoSurf = m_origSurf->vertices;
	m_faces = m_origSurf->faces;
	m_nbFaces_for_edge.clear();

	m_ready_to_draw = false;
	m_pre_iso_ratio = 0.0f;
}

void IsoSurfFrom2Manifold::precompute()
{
	try{
		// collect edges, and nb faces for each edge
		for (unsigned fi = 0; fi < m_faces.size(); ++fi)
		{
			const auto f = m_faces[fi];
			const auto edges_f = util::edgesFromFace(f);
			for (auto eit = edges_f.begin(); eit != edges_f.end(); ++eit)
			{
				m_nbFaces_for_edge[*eit].push_back(fi);
			}
		}

		// build kd-tree on MA & MC points
		if (m_tree)
		{
			delete m_tree;
		}
		m_tree = new KNNTree();

		const auto& MA_vts = m_stg_ptr->m_stSubdiv.getAllVts();
		for (unsigned i = 0; i < MA_vts.size(); ++i)
		{
			const auto& p = MA_vts[i];
			m_tree->insert(Point_and_uint(P3(p[0], p[1], p[2]), i));
		}
		const auto& MC_vts = this->m_stg_ptr->dual_vts;
		for (unsigned i = 0; i < MC_vts.size(); ++i)
		{
			const auto& p = MC_vts[i];
			// NOTE: the vert index of MC vert should be offset by size of MA_vts
			m_tree->insert(Point_and_uint(P3(p[0], p[1], p[2]), i + MA_vts.size()));
		}

		// query the kd-tree to obtain corresp. orig vert -> MA/MC vert
		
		float temp_diff;
		m_closest_vert_for_orig.assign(m_origSurf->vertices.size(), -1);
		for (unsigned i = 0; i < m_origSurf->vertices.size(); ++i)
		{
			//cout << "# of queries done for orig -> MA/MC vert correspondence: "<<i<<endl; // debug

			const auto& query_p = m_origSurf->vertices[i];
			unsigned vi_min = 0;
			float diff_min = std::numeric_limits<float>::max(); // min absolute difference
			K_neighbor_search search1(*m_tree, P3(query_p[0], query_p[1], query_p[2]), 2);
			//cout << "search done."<<endl;

			for (auto q_it = search1.begin(); q_it != search1.end(); q_it++)
			{
				//cout << "# NN points examined: "<<q_it-search1.begin()<<endl;

				auto vi = boost::get<1>(q_it->first);
				if (vi >= MA_vts.size())
				{
					temp_diff = std::abs(
						m_stg_ptr->bt3_MC[vi - MA_vts.size()] - 
						trimesh::dist(query_p, m_stg_ptr->dual_vts[vi - MA_vts.size()])	);
				}
				else
				{
					temp_diff = std::abs(
						m_stg_ptr->bt3MA_vert[vi] - 
						trimesh::dist(query_p, MA_vts[vi])	);
				}

				//cout << "temp_diff = "<<temp_diff<<endl;

				if (temp_diff < diff_min)
				{
					diff_min = std::min(temp_diff, diff_min);
					vi_min = vi;
				}

				//cout << "next point."<<endl;
			}
			m_closest_vert_for_orig[i] = vi_min;

			//cout << "next query."<<endl;
		}

		// get the range of iso values
		m_iso_max = std::max(
			*( std::max_element(m_stg_ptr->bt3_MC.begin(), m_stg_ptr->bt3_MC.end()) ),
			*( std::max_element(m_stg_ptr->bt3MA_vert.begin(), m_stg_ptr->bt3MA_vert.end()) )
			);

		// average length of a subset of edges
		m_eps_1 = m_eps_2 = 0.0f;

		for (unsigned fi = 0; fi < m_faces.size(); ++fi)
		{
			const auto& f = m_faces[fi];
			m_eps_2 += trimesh::dist(m_vts_isoSurf[f[0]], m_vts_isoSurf[f[1]]);
		}
		m_eps_2 /= m_faces.size();
		m_eps_2 *= 1.2f;

		for (unsigned ei = 0; ei < m_stg_ptr->m_origG->edges.size(); ++ei)
		{
			const auto& e = m_stg_ptr->m_origG->edges[ei];
			m_eps_1 += trimesh::dist(m_stg_ptr->m_origG->vts[e[0]], m_stg_ptr->m_origG->vts[e[1]]);
		}
		m_eps_1 /= m_stg_ptr->m_origG->edges.size();
		m_eps_1 *= 5;

		cout << "m_eps_1 = "<<m_eps_1<< endl;
		cout << "m_eps_2 = "<<m_eps_2<< endl;

		m_ready_to_draw = false;

		cout << "iso-surface precompute done." << endl;
	}
	catch(const std::exception& se)
	{
		std::cerr <<
			"General error: " <<
			se.what() << std::endl;
	}
}

void IsoSurfFrom2Manifold::updateIsoSurf(float _iso_ratio, bool _hide_snapped_faces)
{
	//// refine the triangulation if necessary
	//if (_iso_ratio > m_pre_iso_ratio)
	//	refineTriangulation();
	m_pre_iso_ratio = _iso_ratio;
	m_iso_value = (m_iso_max - 0.0f) * _iso_ratio + 0.0f;

	// update iso-surf vertices to isovalue = _iso surface
	updateIsoSurf(m_iso_value, _hide_snapped_faces, 1);

	// upload geometry to drawer
	/*upload_data_to_drawer();*/
}

void IsoSurfFrom2Manifold::updateIsoSurf(float _iso, bool _hide_snapped_faces, int)
{
	m_iso_value = _iso;
	update_iso_surf(_hide_snapped_faces);
}

/// private members

void IsoSurfFrom2Manifold::refineTriangulation(bool _propagate)
{
	// prioritize edges based on their estimated error
	float err;
	for (auto e_iter = m_nbFaces_for_edge.begin(); e_iter != m_nbFaces_for_edge.end(); ++e_iter)
	{
		err = estimate_error(e_iter->first);
		m_refine_q.push(std::make_tuple(e_iter->first, err));
	}

	vector<int> face_tag; // 1 = unvisited; 2 = keep; 3 = discard;
	face_tag.assign(m_faces.size(), 1);
	int nb_f1_id, nb_f2_id;
	TriFace nb_f1, nb_f2;
	vector<int> nb_faces;
	TriEdge e1, e2, e3, e4;
	TriFace f1, f2, f3, f4;
	int f1_i, f2_i, f3_i, f4_i;
	TriEdge nbF1_e1, nbF1_e2, nbF2_e1, nbF2_e2;

	// refine each face if necessary
	while (!m_refine_q.empty())
	{
		if (m_refine_q.size() % 200000 == 0 /*1*/ ) // debug
			cout << "refine_q size: " << m_refine_q.size() << endl;

		// 1. check if this edge needs refine. if so, refine it & its nb face sharing the edge to split
		auto top_elem = m_refine_q.top();
		const auto e_split = std::get<0>(top_elem);
		float err = std::get<1>(top_elem);
		m_refine_q.pop();
		if ( /*!need_split_1(e_split)*/ /*!need_split_2(err)*/ !need_split_3(e_split) )
			continue;		

		nb_faces = m_nbFaces_for_edge.find(e_split)->second;
		nb_f1_id = nb_faces[0];
		nb_f2_id = nb_faces[1];
		nb_f1 = m_faces[nb_f1_id];
		nb_f2 = m_faces[nb_f2_id];
		bool bad_face = false;

		/// splitting this edge

		// refine f and its nb by splitting e in the middle
		// add the new faces to refine_q
		m_vts_isoSurf.push_back(split_isoSurf_edge(e_split));
		m_vts_origSurf.push_back(split_origSurf_edge(e_split));
		int v_idx = m_vts_isoSurf.size() - 1;

		refine_face(nb_f1, e_split, v_idx, f1, f2);
		m_faces.push_back(f1);
		f1_i = m_faces.size() - 1;
		m_faces.push_back(f2);
		f2_i = m_faces.size() - 1;

		refine_face(nb_f2, e_split, v_idx, f3, f4);
		m_faces.push_back(f3);
		f3_i = m_faces.size() - 1;
		m_faces.push_back(f4);
		f4_i = m_faces.size() - 1;

		// add new edge-face adjacency 
		e1 = util::makeEdge(v_idx, e_split[0]);
		e2 = util::makeEdge(v_idx, e_split[1]);
		int nbF1_oppoV = util::oppositeVert(nb_f1, e_split);
		int nbF2_oppoV = util::oppositeVert(nb_f2, e_split);
		e3 = util::makeEdge(v_idx, nbF1_oppoV);
		e4 = util::makeEdge(v_idx, nbF2_oppoV);

		auto& new_nbs_1 = m_nbFaces_for_edge[e1];
		new_nbs_1.push_back(f1_i);
		new_nbs_1.push_back(f3_i);
		auto& new_nbs_2 = m_nbFaces_for_edge[e2];
		new_nbs_2.push_back(f2_i);
		new_nbs_2.push_back(f4_i);
		auto& new_nbs_3 = m_nbFaces_for_edge[e3];
		new_nbs_3.push_back(f1_i);
		new_nbs_3.push_back(f2_i);
		auto& new_nbs_4 = m_nbFaces_for_edge[e4];
		new_nbs_4.push_back(f3_i);
		new_nbs_4.push_back(f4_i);
		nbF1_e1 = util::makeEdge(nbF1_oppoV, e_split[0]);
		nbF1_e2 = util::makeEdge(nbF1_oppoV, e_split[1]);
		nbF2_e1 = util::makeEdge(nbF2_oppoV, e_split[0]);
		nbF2_e2 = util::makeEdge(nbF2_oppoV, e_split[1]);
		update_nb_face(nbF1_e1, nb_f1_id, f1_i);
		update_nb_face(nbF1_e2, nb_f1_id, f2_i);
		update_nb_face(nbF2_e1, nb_f2_id, f3_i);
		update_nb_face(nbF2_e2, nb_f2_id, f4_i);

		// remove old edge-face adjacency
		m_nbFaces_for_edge.erase(m_nbFaces_for_edge.find(e_split));
		// find corresp. for new vertex
		addTo_closest(v_idx, e_split);

		bad_face = true;

		if (_propagate) // recursively split new edges?
		{
			// add four new edges properly to refine_q
			// b.c. they need to be checked in later iterations
			m_refine_q.push(std::make_tuple(e1, estimate_error(e1)));
			m_refine_q.push(std::make_tuple(e2, estimate_error(e2)));
			m_refine_q.push(std::make_tuple(e3, estimate_error(e3)));
			m_refine_q.push(std::make_tuple(e4, estimate_error(e4)));
		}

		//// 2. after splitting edge, mark orig nb faces as "discard"
		////    b.c. they no longer exists
		//face_tag[nb_f1_id] = 3;
		//face_tag[nb_f2_id] = 3;
	} // end while refine_q !empty

	// finalize all faces. (all nb faces in adjacency map)
	vector<TriFace> good_faces;
	good_faces.reserve(m_faces.size());
	map<int, int> oldFi_newFi; // used for renaming face index

	for (auto iter = m_nbFaces_for_edge.begin(); iter != m_nbFaces_for_edge.end(); ++ iter)
	{
		nb_faces = iter->second;
		auto find_iter = oldFi_newFi.find(nb_faces[0]);
		if (find_iter == oldFi_newFi.end()) // add this face & record renaming rule
		{
			good_faces.push_back(m_faces[nb_faces[0]]);
			oldFi_newFi[nb_faces[0]] = good_faces.size() - 1;
		}
		find_iter = oldFi_newFi.find(nb_faces[1]);
		if (find_iter == oldFi_newFi.end()) // add this face & record renaming rule
		{
			good_faces.push_back(m_faces[nb_faces[1]]);
			oldFi_newFi[nb_faces[1]] = good_faces.size() - 1;
		}
	}
	// renaming nb face index
	for (auto iter = m_nbFaces_for_edge.begin(); iter != m_nbFaces_for_edge.end(); ++iter)
	{
		auto& nbfaces = iter->second;
		auto find_it = oldFi_newFi.find(nbfaces[0]);
		if (find_it != oldFi_newFi.end())
		{
			nbfaces[0] = find_it->second; // rename
		}
		find_it = oldFi_newFi.find(nbfaces[1]);
		if (find_it != oldFi_newFi.end())
		{
			nbfaces[1] = find_it->second; // rename
		}
	}

	m_faces = good_faces;

	cout << "iso-surface refine_triangulation() done." << endl;
}

void IsoSurfFrom2Manifold::update_iso_surf(bool _hide_snapped_faces)
{
	// update iso surface geometry
	const auto& MC_vts = m_stg_ptr->dual_vts;
	const auto& MA_vts = m_stg_ptr->m_stSubdiv.getAllVts();
	int closest_vi = -1;
	float bt3;
	TriPoint target_v;
	trimesh::vec3 dir;
	// has the vert snapped onto the MA?
	vector<bool> vert_on_MA(m_vts_isoSurf.size(), false);
	for (unsigned vi = 0; vi < m_vts_isoSurf.size(); ++vi)
	{
		closest_vi = m_closest_vert_for_orig[vi];
		//if (_ratio < m_pre_iso_ratio) // move toward orig surf
		//{
		//	target_v = m_vts_origSurf[vi];
		//}
		//else // move toward MC/MA vert
		//{
		if (closest_vi >= MA_vts.size())
		{
			closest_vi -= MA_vts.size();
			//bt3 = m_stg_ptr->bt3_MC[closest_vi];
			bt3 = trimesh::dist(MC_vts[closest_vi], m_vts_origSurf[vi]);
			target_v = MC_vts[closest_vi];
		}
		else
		{
			//bt3 = m_stg_ptr->bt3MA_vert[closest_vi] ;
			bt3 = trimesh::dist(MA_vts[closest_vi], m_vts_origSurf[vi]);
			target_v = MA_vts[closest_vi];
		}
		//}

		dir = trimesh::normalize(trimesh::vec3(m_vts_origSurf[vi] - target_v));
		//dir = trimesh::normalize(trimesh::vec3(m_vts_isoSurf[vi] - target_v));
		float t = bt3 - m_iso_value ;
		vert_on_MA[vi] = (t < 0.0f ? true : false);
		m_vts_isoSurf[vi] = target_v + std::max(t, 0.0f) * dir;
	}

	// upload iso surface geometry to gpu
	auto drawer = std::dynamic_pointer_cast<MeshDrawer>(m_isoSurf_drawer);
	drawer->setRenderMode(MeshDrawer::PER_FACE);
	vector<TriFace> faces_to_draw;

	if (_hide_snapped_faces)
	{
		faces_to_draw.reserve(m_faces.size());
		// only show faces snapped onto MA
		for (unsigned fi = 0; fi < m_faces.size(); ++fi)
		{
			const auto& f = m_faces[fi];
			if ( (int)vert_on_MA[f[0]] + (int)vert_on_MA[f[1]] + (int)vert_on_MA[f[2]] < 2 )
				faces_to_draw.push_back(f);
		}
	}
	else
		faces_to_draw = m_faces;
	drawer->setPointsPerFace(m_vts_isoSurf, faces_to_draw);
	/*drawer->setPoints(m_vts_isoSurf);
	drawer->setFaces(faces_to_draw, true);*/

	// set isosurf color
	TriColor iso_color;
	TriColor cur_level_color = util::GetColour(m_iso_value - 0.0f, 0.0f, 1.0f);
	//float* vts_color_data = new float[m_vts_isoSurf.size() * 3];
	float* vts_scalar = new float[m_vts_isoSurf.size()];
	float* face_color_data = new float[faces_to_draw.size() * 3];
	for (unsigned i = 0; i < m_vts_isoSurf.size(); ++i)
	{
		closest_vi = m_closest_vert_for_orig[i];
		if (closest_vi >= MA_vts.size())
		{
			closest_vi -= MA_vts.size();
			target_v = MC_vts[closest_vi];
		}
		else
		{
			target_v = MA_vts[closest_vi];
		}
		bt3 = trimesh::dist(target_v, m_vts_origSurf[i]);
		if (m_iso_value > bt3)
			vts_scalar[i] = bt3;
		else
			vts_scalar[i] = m_iso_value;
		/*if (m_iso_value > bt3)
			iso_color = util::GetColour(bt3, 0.0f, m_iso_max);
		else
			iso_color = cur_level_color;

		vts_color_data[3*i+0] = iso_color[0];
		vts_color_data[3*i+1] = iso_color[1];
		vts_color_data[3*i+2] = iso_color[2];*/
	}
	for (unsigned fi = 0; fi < faces_to_draw.size(); ++fi)
	{
		const auto& f = faces_to_draw[fi];
		TriColor c = util::GetColour(
			(vts_scalar[f[0]]+vts_scalar[f[1]]+vts_scalar[f[2]]) / 3.0f, 
			0.0f, 
			m_iso_max );
		face_color_data[3*fi+0] = c[0];
		face_color_data[3*fi+1] = c[1];
		face_color_data[3*fi+2] = c[2];
	}
	//drawer->setPerVertColor(color_data, m_vts_isoSurf.size());
	//delete [] color_data;
	//color_data = nullptr;
	drawer->setPerFaceColor(face_color_data, faces_to_draw.size());
	//delete [] vts_color_data;
	delete [] vts_scalar;
	delete [] face_color_data;

	// set isosurf normal
	/* float* normal_data = new float[m_vts_isoSurf.size() * 3];
	memset((void*)normal_data, 0, m_vts_isoSurf.size()*3*sizeof(float));
	int* nbFace_cnt = new int[m_vts_isoSurf.size()];
	memset((void*)nbFace_cnt, 0, m_vts_isoSurf.size()*sizeof(int));
	trimesh::vec3 nml;
	for (unsigned fi = 0; fi < m_faces.size(); ++fi)
	{
		const auto& f = m_faces[fi];
		nml = trimesh::normalize(trimesh::trinorm(
			m_vts_isoSurf[f[0]], m_vts_isoSurf[f[1]], m_vts_isoSurf[f[2]] ));
		for (unsigned j = 0; j < 3; ++j)
		{
			normal_data[3*f[j] + 0] += nml[0];
			normal_data[3*f[j] + 1] += nml[1];
			normal_data[3*f[j] + 2] += nml[2];

			nbFace_cnt[f[j]] ++;
		}
	}
	for (unsigned vi = 0; vi < m_vts_isoSurf.size(); ++vi)
	{
		normal_data[3*vi+0] /= nbFace_cnt[vi];
		normal_data[3*vi+1] /= nbFace_cnt[vi];
		normal_data[3*vi+2] /= nbFace_cnt[vi];
	} 

	drawer->setPerVertNormal(normal_data, m_vts_isoSurf.size());
	delete [] normal_data;
	delete [] nbFace_cnt;
	*/

	float* normal_data = new float[faces_to_draw.size() * 3];
	trimesh::vec3 nml;
	
	for (unsigned fi = 0; fi < faces_to_draw.size(); ++fi)
	{
		const auto& f = faces_to_draw[fi];
		/*nml = trimesh::normalize(trimesh::trinorm(
			m_vts_isoSurf[f[0]], m_vts_isoSurf[f[1]], m_vts_isoSurf[f[2]] ));*/
		nml = trimesh::normalize(
			(m_vts_isoSurf[f[0]] - m_vts_isoSurf[f[1]]).cross(m_vts_isoSurf[f[2]] - m_vts_isoSurf[f[1]])  
			);
		normal_data[3*fi + 0] = nml[0];
		normal_data[3*fi + 1] = nml[1];
		normal_data[3*fi + 2] = nml[2];
	}

	drawer->setPerFaceNormal(normal_data, faces_to_draw.size());
	delete [] normal_data;

	this->m_ready_to_draw = true;

	//cout << "iso-surface update_iso_surf() done." << endl;
}

float IsoSurfFrom2Manifold::estimate_error(const TriEdge& _e)
{
	TriPoint target_1, target_2;
	int v1, v2;
	v1 = m_closest_vert_for_orig[_e[0]];
	v2 = m_closest_vert_for_orig[_e[1]];
	const auto& ma_vts = m_stg_ptr->m_stSubdiv.getAllVts();
	const auto& mc_vts = m_stg_ptr->dual_vts;
	if (v1 >= ma_vts.size())
		target_1 = mc_vts[v1 - ma_vts.size()];
	else
		target_1 = ma_vts[v1];
	if (v2 >= ma_vts.size())
		target_2 = mc_vts[v2 - ma_vts.size()];
	else
		target_2 = ma_vts[v2];
	float d1 = trimesh::dist(target_1, target_2);

	return d1;
}

bool IsoSurfFrom2Manifold::need_split_1(const TriEdge& _e)
{
	// function on_same_edge_or_face(u, e0):
	//   u is one of e0's ends, or there existing face <u, p, q>
	// function is_face_vertex(u, fi):
	//   u is one of face fi's vertices
	// function stVts_on_same_edge_or_face(u, v):
	//   return if st. vts u and v are on the same edge or 
	//   different edges of the same face

	// <u, v> <- vts on MC/MA e's ends mapped to
	// if u & v are all orig vts
	//   check if there exists orig edge <u, v>;
	// if u is orig, v is steiner
	//   e0<p, q> <- v's residing edge;
	//   return ! on_same_edge_or_face(u, e0);
	// if u is orig, v is edge dual
	//   e0<p, q> <- v's residing edge;
	//   return ! on_same_edge_or_face(u, e0);
	// if u is orig, v is face dual
	//   fi <- v's from_face;
	//   return ! is_face_vertex(u, fi);
	// if u & v are all steiner vts
	//   e0, e1 <- u and v's residing edge resp.
	//   return ! same_edges_or_on_same_face(e0, e1);
	// if u is steiner and v is edge dual
	//   e0, e1 <- u and v's residing edge resp.
	//   return ! same_edges_or_on_same_face(e0, e1);
	// if u is steiner and v is face dual
	//   e0 <- u's residing edge
	//   fi <- v's from face
	//   return ! is_face_edge(e0, fi);
	// if u & v are all edge duals
	//   e0, e1 <- u and v's residing edge resp.
	//   return ! same_edges_or_on_same_face(e0, e1);
	// if u is edge dual and v is face dual
	//   e0 <- u's residing edge
	//   fi <- v's from face
	//   return ! is_face_edge(e0, fi);
	// if u and v are all face duals
	//   return ! from_the_same_face(u's residing f, v's residing f)

	// 
	// OR
	//
	// if (mappings on MA of e's two ends are farther apart than some eps value)
	// then return true
	// else return false

	TriPoint target_1, target_2;
	int v1, v2;
	v1 = m_closest_vert_for_orig[_e[0]];
	v2 = m_closest_vert_for_orig[_e[1]];
	const auto& ma_vts = m_stg_ptr->m_stSubdiv.getAllVts();
	const auto& mc_vts = m_stg_ptr->dual_vts;
	if (v1 >= ma_vts.size())
		target_1 = mc_vts[v1 - ma_vts.size()];
	else
		target_1 = ma_vts[v1];
	if (v2 >= ma_vts.size())
		target_2 = mc_vts[v2 - ma_vts.size()];
	else
		target_2 = ma_vts[v2];
	float d1 = trimesh::dist(target_1, target_2);

	float d2 = trimesh::dist(m_vts_isoSurf[_e[0]], m_vts_isoSurf[_e[1]]);
	return (d1 > m_eps_1 && d2 > m_eps_2);
}

bool IsoSurfFrom2Manifold::need_split_2(float _err)
{
	return _err > m_eps_1;
}

bool IsoSurfFrom2Manifold::need_split_3(const TriEdge& _e)
{

	int tgt_1, tgt_2;
	tgt_1 = m_closest_vert_for_orig[_e[0]];
	tgt_2 = m_closest_vert_for_orig[_e[1]];

// 	if (_e[0] == 25 && _e[1] == 6138)
// 		int stop = 1;
	auto vert_pair = encode_type_for_vts(tgt_1, tgt_2);
	bool ret = false;
	switch ( (vTypePair)vert_pair )
	{
	case ORIG_ORIG:
		ret = !same_face_orig_orig(tgt_1, tgt_2);
		break;
	case ORIG_STEINER:
		ret = !same_face_orig_steiner(tgt_1, tgt_2);
		break;
	case ORIG_EDUAL:
		ret = !same_face_orig_edual(tgt_1, tgt_2);
		break;
	case ORIG_FDUAL:
		ret = !same_face_orig_fdual(tgt_1, tgt_2);
		break;
	case STEINER_STEINER:
		ret = !same_face_steiner_steiner(tgt_1, tgt_2);
		break;
	case STEINER_EDUAL:
		ret = !same_face_steiner_edual(tgt_1, tgt_2);
		break;
	case STEINER_FDUAL:
		ret = !same_face_steiner_fdual(tgt_1, tgt_2);
		break;
	case EDUAL_EDUAL:
		ret = !same_face_edual_edual(tgt_1, tgt_2);
		break;
	case EDUAL_FDUAL:
		ret = !same_face_edual_fdual(tgt_1, tgt_2);
		break;
	case FDUAL_FDUAL:
		ret = !same_face_fdual_fdual(tgt_1, tgt_2);
		break;
	default:
		// doesn't recognize this encoding. rather not split.
		ret = false;
		assert(0);
		break;
	}

	float d2 = trimesh::dist(m_vts_isoSurf[_e[0]], m_vts_isoSurf[_e[1]]);
	return d2 > m_eps_2 && ret;
}

bool IsoSurfFrom2Manifold::need_refine(const TriFace& _f, TriEdge& _e_to_split)
{
	float dummy;
	util::longestEdge(_f[0], _f[1], _f[2], 
		m_vts_isoSurf[_f[0]], m_vts_isoSurf[_f[1]], m_vts_isoSurf[_f[2]], 
		_e_to_split, dummy);

	return need_split_1(_e_to_split);
}

TriPoint IsoSurfFrom2Manifold::split_isoSurf_edge(const TriEdge& _e) 
{
	return trimesh::mix(m_vts_isoSurf[_e[0]], m_vts_isoSurf[_e[1]], 0.5f);
}

TriPoint IsoSurfFrom2Manifold::split_origSurf_edge(const TriEdge& _e) 
{
	return trimesh::mix(m_vts_origSurf[_e[0]], m_vts_origSurf[_e[1]], 0.5f);
}

void IsoSurfFrom2Manifold::refine_face(const TriFace& _f, const TriEdge& _e, int _vidx, TriFace& _f1, TriFace& _f2)
{
	int oppo_v = util::oppositeVert(_f, _e);
	_f1 = util::makeFace(oppo_v, _e[0], _vidx);
	_f2 = util::makeFace(oppo_v, _e[1], _vidx);
}

void IsoSurfFrom2Manifold::addTo_closest(int _v_idx, const  TriEdge& _e)
{
	// _v_idx will not be mapped to the verts that the ends of _e are mapped to
	const auto& query_p = m_vts_isoSurf[_v_idx];
	K_neighbor_search search1(*m_tree, P3(query_p[0], query_p[1], query_p[2]), 3);
	auto iter = search1.begin();
	int vi = boost::get<1>(iter->first);
	int target_1 = m_closest_vert_for_orig[_e[0]];
	int target_2 = m_closest_vert_for_orig[_e[1]];
	while (vi == target_1 || vi == target_2)
	{
		iter++;
		vi = boost::get<1>(iter->first);
	}
	
	m_closest_vert_for_orig.push_back(vi);
}

void IsoSurfFrom2Manifold::update_nb_face(TriEdge _e, int _old_f, int _new_f)
{
	auto& nbs = m_nbFaces_for_edge.find(_e)->second;
	assert(nbs[0] == _old_f || nbs[1] == _old_f);
	nbs[0] == _old_f ? nbs[0] = _new_f : nbs[1] = _new_f;
}