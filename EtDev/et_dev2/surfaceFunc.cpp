#include "surfaceFunc.h"
#include "geometryRepresentation.h"

const float SurfaceFunc::INVALID_VAL = -999.0f;

SurfaceFunc::SurfaceFunc()
{

}

SurfaceFunc::~SurfaceFunc()
{

}

void SurfaceFunc::reset()
{
	m_stg_ptr.reset();
	m_origSurf_graph.reset();
}

void SurfaceFunc::setup(
	const std::shared_ptr<SteinerGraph>& _stg_ptr, 
	const std::shared_ptr<TriMesh>& _3d_surf)
{
	this->m_stg_ptr = _stg_ptr;
	this->m_origSurf_mesh = _3d_surf;
	vector<TriEdge> dummy;
	this->m_origSurf_graph = std::shared_ptr<MyGraph>(
		new MyGraph(_3d_surf->vertices, dummy, _3d_surf->faces)
		);

	this->m_closest_vert_for_MC.clear();
	this->m_closest_vert_for_origSurf.clear();
	this->m_close_vts_for_origSurf.clear();
}

void SurfaceFunc::setCurrentFunction()
{

}

void SurfaceFunc::getSurfaceFunction(
	SurfFuncType _func, 
	vector<float>& _scalar_vals, 
	float _smooth_ratio, 
	bool _do_diffuse/* = false*/
	)
{
	if (/*m_closest_vert_for_MC.empty() || 
		*/m_closest_vert_for_origSurf.empty() && 
		m_close_vts_for_origSurf.empty() )
	{
		//compute_dualVert_origVert_corresp_with_range_search();
		computeCorrespondenceWithKNNSearch(6);

		/*cout << "computing next-vertex along bt2 gradient..."<<endl;
		compute_next_vert_along_bt2_gradient();
		cout << "next-vertex along bt2 gradient computed."<<endl;*/
	}

	float eps = 0.0f;
	switch (_func)
	{
	case SHAPE_DIAM:
		convert_to_shape_diam_eps_range(_smooth_ratio, eps);
		get_shape_diam(_scalar_vals, eps);
		break;
	case SHAPE_WIDTH:
		convert_to_shape_width_eps_range(_smooth_ratio, eps);
		get_shape_width(_scalar_vals, eps);
		break;
	case SHAPE_WIDTH_EXTREMITY:
		convert_to_shape_width_extremity_eps_range(_smooth_ratio, eps);
		get_shape_width_extremity(_scalar_vals, eps);
		break;
	case SHAPE_LENGTH_EXTREMITY:
		eps = numeric_limits<float>::max();
		get_shape_width_extremity(_scalar_vals, eps);
		break;
	default:
		convert_to_shape_width_eps_range(_smooth_ratio, eps);
		get_shape_width(_scalar_vals, eps);
		break;
	}

	float scalar_min = numeric_limits<float>::max();
	float scalar_max = numeric_limits<float>::min();
	for (auto it = _scalar_vals.begin(); it != _scalar_vals.end(); ++it)
	{
		scalar_min = std::min(scalar_min, *it);
		scalar_max = std::max(scalar_max, *it);
	}
	cout << "scalar field obtained. min/max "<<scalar_min<<"/"<<scalar_max<<endl;

	if ( _do_diffuse )
	{
		//diffuse_scalar_field( _scalar_vals );
		m_diffuser.solve(_scalar_vals);
		//this->outputDiffuseSystem(_scalar_vals);
	}

	//cout << "surf func obtained." << endl;
}

int SurfaceFunc::getOrigToMC(int _i)
{
	if (m_closest_vert_for_origSurf.size() <= _i)
		return -1;
	return m_closest_vert_for_origSurf[_i];
}

void SurfaceFunc::outputDiffuseSystem( const vector<float>& _scalars )
{
	ofstream o_file("diffuse_on_surf"); 
	std::stringstream out_string;
	out_string << "{"<< endl; // Big {

	// Begin: output mesh vts
	out_string << "{ (*Begin: mesh vts*)"<< endl; 
	for (unsigned i = 0; i < m_origSurf_graph->vts.size(); ++i)
	{
		const auto & p = m_origSurf_graph->vts[i];
		out_string << p[0]<<","<<p[1]<<","<<p[2];
		if (i != m_origSurf_graph->vts.size()-1)
		{
			out_string <<",";
		}
	}
	out_string << endl;
	out_string << "}, (*End: mesh vts*)"<< endl; 
	// End: output mesh vts

	// Begin: output mesh faces
	out_string << "{ (*Begin: mesh faces*)"<< endl; 
	for (unsigned i = 0; i < m_origSurf_graph->faces.size(); ++i)
	{
		const auto & f = m_origSurf_graph->faces[i];
		out_string << f[0]+1<<","<<f[1]+1<<","<<f[2]+1;
		if (i != m_origSurf_graph->faces.size()-1)
		{
			out_string <<",";
		}
	}
	out_string << endl;
	out_string << "}, (*End: mesh faces*)"<< endl; 
	// End: output mesh faces

	// Begin: output A
	out_string << "{ (*Begin: A*)"<< endl; 
	for ( int i=0; i<A_sm.outerSize(); ++i )
	{
		for ( SpMat::InnerIterator it(A_sm,i); it; ++it )
		{
			out_string << "{"<<it.row()+1<<","<<it.col()+1<<","<<it.value()<<"},";
		}
	}
	out_string.seekp(-1, out_string.cur);
	out_string << "}, (*End: A*)"<< endl; 
	// End: output A

	// Begin: output C
	out_string << "{ (*Begin: C*)"<< endl; 
	for ( int i=0; i<C_sm.outerSize(); ++i )
	{
		for ( SpMat::InnerIterator it(C_sm,i); it; ++it )
		{
			out_string << "{"<<it.row()+1<<","<<it.col()+1<<","<<it.value()<<"},";
		}
	}
	out_string.seekp(-1, out_string.cur);
	out_string << "}, (*End: C*)"<< endl; 
	// End: output C

	// Begin: output f, i.e. the known values
	out_string << "{ (*Begin: f (known values) *)"<< endl; 
	for ( int  i = 0; i < f.size(); ++i )
	{
		out_string << f[i] << ",";
	}
	out_string.seekp(-1, out_string.cur);
	out_string << "}, (*End: f (known values) *)" <<endl;
	// End: output f

	// Begin: output unknown flags
	out_string << "{ (*Begin: flags (0: unknown 1: known) *)"<< endl; 
	for ( int  i = 0; i < m_closest_vert_for_origSurf.size(); ++i )
	{
		out_string << (m_closest_vert_for_origSurf[i] == -1 ? 0 : 1) << ",";
	}
	out_string.seekp(-1, out_string.cur);
	out_string << "}, (*End: flags*)" <<endl;
	// End: output unknown flags

	// Begin: output scalar values
	out_string << "{ (*Begin: scalar values *)"<< endl; 	
	for ( int  i = 0; i < _scalars.size(); ++i )
	{
		out_string << _scalars[i] << ",";
	}
	out_string.seekp(-1, out_string.cur);
	out_string << "} (*End: scalar values*)" <<endl;
	// End: output f

	out_string << "}" << endl;// Big }

	o_file << out_string.rdbuf();
	o_file.close();
}

/*
void SurfaceFunc::compute_next_vert_along_bt2_gradient()
{
	this->m_next_vert_by_bt2.assign(m_stg_ptr->bt2MA_vert.size(), -1);
	const auto& allvts = m_stg_ptr->m_stSubdiv.getAllVts();
	vector<unsigned> nbs;
	for (unsigned vi = 0; vi < allvts.size(); ++vi)
	{
		nbs.clear();
		m_stg_ptr->m_stSubdiv.getNbVerts(vi, nbs);
		float max_bt2_in_nbs = -1.0f;
		unsigned next_nb = nbs[0];
		for (unsigned j = 0; j < nbs.size(); ++j)
		{
			auto bt2 = m_stg_ptr->bt2MA_vert[nbs[j]];
			if ( max_bt2_in_nbs < bt2 )
			{
				max_bt2_in_nbs = bt2;
				next_nb = nbs[j];
			}
		}

		if (m_stg_ptr->bt2MA_vert[vi] <= max_bt2_in_nbs)
		{
			m_next_vert_by_bt2[vi] = next_nb;
		}
		else
		{
			m_next_vert_by_bt2[vi] = vi;
			cout << "local maximum/bt2: "<<vi<<"/"<<m_stg_ptr->bt2MA_vert[vi]<<endl;
		}
	}
}

int SurfaceFunc::find_local_maximum(unsigned _vi, float _eps)
{
	set<unsigned> visited;
	const auto& adjacency = m_stg_ptr->m_origG->nbVtsOfVert;
	int cur_vi = _vi;
	float cur_bt2 = m_stg_ptr->bt2MA_vert[cur_vi];

	if ( cur_bt2 >= _eps )
		return cur_vi;

	visited.insert(cur_vi);
	while (true)
	{
		unsigned next_vi = this->m_next_vert_by_bt2[cur_vi];
		float next_bt2 = m_stg_ptr->bt2MA_vert[next_vi];
		if ( next_vi == cur_vi || visited.count(next_vi) > 0 || next_bt2 >= _eps )
		{
			break;
		}
		visited.insert(next_vi);
		cur_bt2 = next_bt2;
		cur_vi = next_vi;
	}
	visited.clear();

	return cur_vi;
}
*/

// This implementation starts from 3d surface to find corresponding knn on MC/MA
//void SurfaceFunc::compute_dualVert_origVert_corresp()
//{
//	cout << "computing correspondences... "<<endl;
//
//	// build kd-tree on MC points
//	KNNTree* tree = new KNNTree();
//	const auto& ma = *m_stg_ptr->m_origG;
//	for (unsigned i = 0; i < ma.vts.size(); ++i)
//	{
//		const auto& p = ma.vts[i];
//		tree->insert(Point_and_uint(P3(p[0], p[1], p[2]), i));
//	}
//
//	// for each orig vert, find its closest MA point
//	m_closest_vert_for_origSurf.assign(m_origSurf->vertices.size(), -1);
//	for (unsigned i = 0; i < m_origSurf->vertices.size(); ++i)
//	{
//		const auto& query_p = m_origSurf->vertices[i];
//		K_neighbor_search search1(*tree, P3(query_p[0], query_p[1], query_p[2]), 20);
//		unsigned vi_min = 0;
//		float diff_min = std::numeric_limits<float>::max(); // min absolute difference
//		for (auto q_it = search1.begin(); q_it != search1.end(); ++q_it)
//		{
//			auto vi_on_ma = boost::get<1>(q_it->first);
//			auto temp_diff = std::abs(
//				m_stg_ptr->bt3MA_vert[vi_on_ma] - 
//				trimesh::dist(query_p, ma.vts[vi_on_ma])
//				);
//			if (temp_diff < diff_min)
//			{
//				diff_min = std::min(temp_diff, diff_min);
//				vi_min = vi_on_ma;
//			}
//		}
//		m_closest_vert_for_origSurf[i] = vi_min;
//		/*K_neighbor_search::iterator nb_it = search1.begin();
//		unsigned vi = boost::get<1>(nb_it->first);
//		m_closest_vert_for_origSurf[i] = vi;*/
//	}
//
//	delete tree;
//
//	cout << "correspondence computed."<<endl;
//}

float SurfaceFunc::compute_mvc_weight(unsigned _vi, unsigned _vj, std::shared_ptr<MyGraph> _g)
{
	auto e = util::makeEdge(_vi, _vj);
	const auto& nb_faces = _g->getNbFaces( e );
	float theta[2] = {0.0f, 0.0f};
	TriPoint p0 = _g->vts[_vi];
	TriPoint p1 = _g->vts[_vj];
	for ( int i = 0; i < 2; ++i )
	{
		const auto& f = _g->faces[nb_faces[i]];
		auto vk = util::oppositeVert( f, e );
		
		TriPoint p2 = _g->vts[vk];
		theta[i] = acos( 
			trimesh::normalize(vec3(p1 - p0)) DOT trimesh::normalize(vec3(p2 - p0)) 
			);
	}

	return ( std::tan(theta[0] / 2) + tan(theta[1] / 2) ) / trimesh::dist(p0, p1);
}

// the system: A.x = b, where b = C.f, and f 
// is only determined in run time (the known scalars)
void SurfaceFunc::setupDiffusionSystem()
{
	vector<SpTriplet> A,C;
	A.clear();
	map<unsigned, unsigned> vert_newOrder_map;

	// separate vts with and w/o correspondence
	int unknown_vts_cnt = 0, known_vts_cnt = 0;
	for ( unsigned vi = 0; vi < m_origSurf_graph->vts.size(); ++vi )
	{
		if ( m_close_vts_for_origSurf[vi].empty() )
		{
			vert_newOrder_map[vi] = unknown_vts_cnt;
			unknown_vts_cnt ++;
		}
		else
		{
			vert_newOrder_map[vi] = known_vts_cnt;
			known_vts_cnt ++;
		}
	}

	if ( unknown_vts_cnt == 0 ) // there is no unknown vars.
		return;

	// construct A & C
	float wij_max = 0.0f;
	float wij_min = numeric_limits<float>::max();
	for ( unsigned vi = 0; vi < m_origSurf_graph->vts.size(); ++vi )
	{
		if ( m_close_vts_for_origSurf[vi].empty() )
		{
			auto i = vert_newOrder_map[vi]; // new order for vi
			float sum_w = 0.0f; // weight for wii
			const auto& nbs = m_origSurf_graph->nbVtsOfVert[vi];
			for ( auto it = nbs.begin(); it != nbs.end(); ++it )
			{
				auto j = vert_newOrder_map[*it]; // new order for nb *it

				// compute the weight for this nb
				float wij = compute_mvc_weight(vi, *it, m_origSurf_graph);
				if( ::is_nan(wij) == true )
					cout << "w("<<i<<","<<j<<")="<<wij<<", ";
				sum_w += wij;
				wij_min = std::min(wij_min, wij);
				wij_max = std::max(wij_max, wij);

				// put the weight for nb in the proper sparse matrix (A or C)
				if ( m_close_vts_for_origSurf[*it].empty() )
				{
					A.push_back( SpTriplet(i, j, -wij) );
				}
				else
				{
					C.push_back( SpTriplet(i, j, wij) );
				}
			}

			// put wii into _A
			A.push_back( SpTriplet(i, i, sum_w) );
		}
	}
	cout << endl;
	cout <<"wij range: "<<wij_min<<"/"<<wij_max<<endl;

	f.resize(known_vts_cnt);
	A_sm.resize(unknown_vts_cnt, unknown_vts_cnt);
	A_sm.setFromTriplets( A.begin(), A.end() );
	C_sm.resize(unknown_vts_cnt, known_vts_cnt);
	C_sm.setFromTriplets( C.begin(), C.end() );
	this->decomp_of_A.compute(A_sm);
	auto ret = decomp_of_A.info();
	if (ret == Eigen::Success)
	{
		cout << "cholesky computation success."<<endl;
	}
	else if (ret == Eigen::NumericalIssue)
	{
		cout << "cholesky numerical issue!"<<endl;
	}
}

void SurfaceFunc::solve_diffusion_system(
	const vector<float>& _known,
	vector<float>& _unknown
	)
{
	//cout << "f: ";
	for (auto it = _known.begin(); it != _known.end(); ++it)
	{
		f[it - _known.begin()] = *it;
		//cout << *it <<", ";
	}
	//cout << endl;

	auto b = C_sm * f;
	Eigen::VectorXd x = decomp_of_A.solve( b );
	/*vector<float> x(10);*/
	_unknown.resize( x.size() );
	cout << "x: ";
	for ( unsigned i = 0; i < x.size(); ++i )
	{
		_unknown[i] = x[i];
		//cout << x[i]<<", ";
	}
	cout << endl;
}

void SurfaceFunc::diffuse_scalar_field( vector<float>& _scalar_vals )
{
	// known scalars
	vector<float> known_scalars;
	for (unsigned i = 0; i < _scalar_vals.size(); i ++)
	{
		if ( _scalar_vals[i] != INVALID_VAL )
		{
			known_scalars.push_back( _scalar_vals[i] );
		}
	}

	if (known_scalars.empty())
	{
		cout << "Cannot perform diffusion: there is no known scalar."<<endl;
		return;
	}

	// solve the diffusion system
	vector<float> unknown_scalars;
	solve_diffusion_system(known_scalars, unknown_scalars);

	int unknown_cnt = 0;
	float val_max = numeric_limits<float>::min();
	float val_min = numeric_limits<float>::max();
	for (unsigned i = 0; i < _scalar_vals.size(); i ++)
	{
		if ( _scalar_vals[i] == INVALID_VAL )
		{
			auto val = unknown_scalars[unknown_cnt];
			_scalar_vals[i] = val;
			unknown_cnt ++;
			val_max = std::max(val_max, val);
			val_min = std::min(val_min, val);
		}
	}

	cout << "unknown scalars (solved) range: "<< val_min<<"/"<<val_max<<endl;
}

//
// This implementation starts from MC/MA to find corresponding knn on 3d surface
//
//void SurfaceFunc::compute_dualVert_origVert_corresp()
//{
//	cout << "computing correspondences... "<<endl;
//	/* compute intersection for MA faces */
//	const auto& bt3_vts = m_stg_ptr->bt3MA_vert;
//	const auto& bt2_vts = m_stg_ptr->bt2MA_vert;
//	vector<TriPoint> isects_per_face;
//	m_stg_ptr->computeIntersectionsForFaces(bt3_vts, &isects_per_face, NULL);
//	const auto& from_face = m_stg_ptr->m_dualEdge_from_which_face;
//
//	// kd tree built upon m3d's vertices
//	KNNTree* tree = new KNNTree();
//	TriPoint p;
//	for (unsigned i = 0; i < this->m_origSurf_mesh->vertices.size(); ++i)
//	{
//		p = m_origSurf_mesh->vertices[i];
//		tree->insert( Point_and_uint(P3(p[0], p[1], p[2]), i) );
//	}
//
//	this->m_closest_vert_for_MC.assign(2*m_stg_ptr->dual_edges.size(), -1);
//	//this->m_close_vts_for_origSurf.resize(m_stg_ptr->dual_edges.size());
//	const auto& ma = m_stg_ptr->m_origG;
//	const auto& ma_vts = ma->vts;
//	for (unsigned i  =0; i < m_stg_ptr->dual_edges.size(); ++i)
//	{
//		auto e = m_stg_ptr->dual_edges[i];
//		auto query_p = trimesh::mix(m_stg_ptr->dual_vts[e[0]], m_stg_ptr->dual_vts[e[0]], 0.5f);
//		K_neighbor_search search1(*tree, P3(query_p[0], query_p[1], query_p[2]), 1);
//		K_neighbor_search::iterator nb_it = search1.begin();
//		unsigned vi = boost::get<1>(nb_it->first);
//		m_closest_vert_for_MC[2*i + 0] = vi;
//		auto first_closest = m_origSurf_mesh->vertices[vi];
//
//		// is the face containing the dual edge good enough to be used as supporting face?
//		if (isects_per_face[from_face[i]*2][0] == numeric_limits<float>::max())
//		{
//			// bad face. can't derive another closest point.
//			m_closest_vert_for_MC[2*i + 1] = vi;
//		}
//		else
//		{
//			// good face. derive another closest point using face normal.
//			auto MA_f = ma->faces[from_face[i]];
//			auto nml = trimesh::normalize(
//				(ma_vts[MA_f[0]] - ma_vts[MA_f[1]]).cross(ma_vts[MA_f[0]] - ma_vts[MA_f[2]]));
//			auto mirror_len = 2*(first_closest - ma_vts[MA_f[0]]).dot(nml);
//			// locate the second closest
//			query_p = first_closest - mirror_len * nml; 
//			K_neighbor_search search2(*tree, P3(query_p[0], query_p[1], query_p[2]), 1);
//			nb_it = search2.begin();
//			vi = boost::get<1>(nb_it->first);
//			m_closest_vert_for_MC[2*i + 1] = vi;
//		}
//	}
//	delete tree;
//
//	// then assign each vert v on orig surf a vert u of the MC edge mapped to v 
//	// s.t. v is an approx. closest vert for u
//	m_closest_vert_for_origSurf.assign(m_origSurf_mesh->vertices.size(), -1);
//	const auto& mc_next = m_stg_ptr->burnNext_medialCurve;
//	for (unsigned i = 0; i < m_stg_ptr->dual_edges.size(); ++i)
//	{
//		TriEdge e = m_stg_ptr->dual_edges[i];
//		int orig_vi = m_closest_vert_for_MC[2*i + 0];
//		int mc_vi_mapped = mc_next[e[0]] == e[1] ? e[0] : e[1];
//		m_closest_vert_for_origSurf[orig_vi] = mc_vi_mapped;
//
//		orig_vi = m_closest_vert_for_MC[2*i + 1];
//		m_closest_vert_for_origSurf[orig_vi] = mc_vi_mapped;
//	}
//
//	cout << "correspondence computed."<<endl;
//}

void SurfaceFunc::setup_diffuser()
{
	//setupDiffusionSystem();
	//UniformWeight w_func;
	MVCWeight w_func;
	w_func.setup(m_origSurf_graph);
	vector<bool> is_known(m_close_vts_for_origSurf.size(), true);
	for (auto i = 0; i < m_close_vts_for_origSurf.size(); ++i)
	{
		if ( m_close_vts_for_origSurf[i].empty() )
			is_known[i] = false;
	}
	m_diffuser.setupDiffuseSystem(m_origSurf_graph.get(), w_func, is_known);
}

void SurfaceFunc::computeCorrespondenceWithKNNSearch(
	int _n
	)
{
	cout << "computing correspondences... "<<endl;
	/* compute intersection for MA faces */
	const auto& bt3_vts = m_stg_ptr->bt3MA_vert;
	const auto& bt2_vts = m_stg_ptr->bt2MA_vert;

	// kd tree built upon m3d's vertices
	KdTree* tree = new KdTree();
	TriPoint p;
	for (unsigned i = 0; i < m_origSurf_graph->vts.size(); ++i)
	{
		p = m_origSurf_graph->vts[i];
		tree->insert( Point_and_uint(P3(p[0], p[1], p[2]), i) );
	}

	this->m_close_vts_for_origSurf.clear();
	this->m_close_vts_for_origSurf.resize(m_origSurf_graph->vts.size());
	const auto& bt3_mc = m_stg_ptr->bt3_MC;
	const auto& ma = m_stg_ptr->m_origG;
	const auto& ma_vts = ma->vts;
	const auto& mc_vts = m_stg_ptr->dual_vts;
	m_origSurf_mesh->need_bbox();
	for (unsigned i  =0; i < mc_vts.size(); ++i)
	{
		float r = bt3_mc[i];
		P3 query_p(mc_vts[i][0], mc_vts[i][1], mc_vts[i][2]);
		K_neighbor_search search1( *tree, query_p, std::max(_n, 1) );

		float r_max = 0.0f;
		for (auto nb_it = search1.begin(); nb_it != search1.end(); ++nb_it)
		{
			r_max = std::max((float)nb_it->second, r_max);
		}
		for ( auto nb_it = search1.begin(); nb_it != search1.end(); ++nb_it )
		{
			const auto& pr = nb_it->first;
			auto v_on_surf = boost::get<0>(pr);
			auto vid_surf = boost::get<1>(pr);
			auto d = CGAL::sqrt(CGAL::squared_distance( v_on_surf, query_p ));
			float w = gaussian( d, r, std::abs(r_max-r) * 2.0f );
			m_close_vts_for_origSurf[vid_surf].push_back(
				std::make_pair(i, w) );
		}
	}
	delete tree;

	int cnt_zero_close = 0;
	for (auto close_vts_it = m_close_vts_for_origSurf.begin();
		close_vts_it != m_close_vts_for_origSurf.end(); ++close_vts_it)
	{
		float w_sum = 0.0f;
		for (auto it = close_vts_it->begin(); it != close_vts_it->end(); ++it)
		{
			w_sum += it->second;
		}
		if ( ::almost_zero( w_sum, 0.000000001f ) )
			close_vts_it->clear();

		if (close_vts_it->empty())
			cnt_zero_close ++;
	}
	cout << "# vts that have zero close vts on MC: "<<cnt_zero_close<<endl;

	cout << "correspondence computed." <<endl;

	cout << "setting up diffusion system..."<<endl;
	setup_diffuser();
	cout << "diffusion system setup."<<endl;
}

// Use sphere-search to build correspondence from MC to original 3d surface
void SurfaceFunc::computeCorrespondenceWithRangeSearch(
	float _eps_ratio, int _k /*= -1*/
	)
{
	cout << "computing correspondences... "<<endl;
	/* compute intersection for MA faces */
	const auto& bt3_vts = m_stg_ptr->bt3MA_vert;
	const auto& bt2_vts = m_stg_ptr->bt2MA_vert;

	// kd tree built upon m3d's vertices
	KdTree* tree = new KdTree();
	TriPoint p;
	for (unsigned i = 0; i < m_origSurf_graph->vts.size(); ++i)
	{
		p = m_origSurf_graph->vts[i];
		tree->insert( Point_and_uint(P3(p[0], p[1], p[2]), i) );
	}

	this->m_close_vts_for_origSurf.clear();
	this->m_close_vts_for_origSurf.resize(m_origSurf_graph->vts.size());
	const auto& bt3_mc = m_stg_ptr->bt3_MC;
	const auto& ma = m_stg_ptr->m_origG;
	const auto& ma_vts = ma->vts;
	const auto& mc_vts = m_stg_ptr->dual_vts;
	m_origSurf_mesh->need_bbox();
	float eps_r = m_origSurf_mesh->bbox.radius() * _eps_ratio;
	vector<Point_and_uint> knn_pts;
	unsigned knn_cnt_min = 999;
	unsigned knn_cnt_max = 0;
	for (unsigned i  =0; i < mc_vts.size(); ++i)
	{
		float r = bt3_mc[i];
		float r_puffed = std::max(r + eps_r, 0.0f); // puffed up to include more points
		P3 center(mc_vts[i][0], mc_vts[i][1], mc_vts[i][2]);
		FuzySphere3 query_sphere( center, r_puffed );
		knn_pts.clear();
		tree->search(std::back_inserter(knn_pts), query_sphere);
		auto most_knn = 
			_k <= 0 ? knn_pts.size() : std::min((size_t)_k, knn_pts.size());

		for ( unsigned j = 0; j < most_knn; ++j )
		{
			const auto& pr = knn_pts[j];
			auto v_on_surf = boost::get<0>(pr);
			auto vid_surf = boost::get<1>(pr);
			auto d = CGAL::sqrt(CGAL::squared_distance( v_on_surf, center ));
			float w = gaussian( d, r, std::abs(r_puffed-r) * 2.0f );
			m_close_vts_for_origSurf[vid_surf].push_back(
				std::make_pair(i, w) );
		}

		knn_cnt_min = std::min(knn_cnt_min, (unsigned)knn_pts.size());
		knn_cnt_max = std::max(knn_cnt_max, (unsigned)knn_pts.size());
	}
	delete tree;

	int cnt_zero_close = 0;
	for (auto close_vts_it = m_close_vts_for_origSurf.begin();
		close_vts_it != m_close_vts_for_origSurf.end(); ++close_vts_it)
	{
		float w_sum = 0.0f;
		for (auto it = close_vts_it->begin(); it != close_vts_it->end(); ++it)
		{
			w_sum += it->second;
		}
		if ( ::almost_zero( w_sum, 0.000000001f ) )
			close_vts_it->clear();

		if (close_vts_it->empty())
			cnt_zero_close ++;
	}
	cout << "# vts that have zero close vts on MC: "<<cnt_zero_close<<endl;

	cout << "correspondence computed. knn range: "<<knn_cnt_min<<"/"<<knn_cnt_max<<endl;

	cout << "setting up diffusion system..."<<endl;
	setup_diffuser();
	cout << "diffusion system setup."<<endl;
}

void SurfaceFunc::convert_to_shape_diam_eps_range(const float _smooth_ratio, float& _eps)
{
	// compute the range for shape diameter's eps paramter
	float max_val = numeric_limits<float>::min();
	float min_val = numeric_limits<float>::max();
	const auto& bt1 = m_stg_ptr->bt1_medialCurve;
	const auto& bt3 = m_stg_ptr->bt3_MC;
	for (unsigned i = 0; i < m_stg_ptr->dual_vts.size(); ++i)
	{
		float val = bt1[i] - bt3[i];
		max_val = std::max(max_val, val);
		min_val = std::min(min_val, val);
	}

	_eps = _smooth_ratio * (max_val - min_val) + min_val;
}

void SurfaceFunc::convert_to_shape_width_eps_range(const float _smooth_ratio, float& _eps)
{
	// compute the range for shape width's eps paramter
	float max_val = numeric_limits<float>::min();
	float min_val = numeric_limits<float>::max();
	const auto& bt2 = m_stg_ptr->bt2_MC;
	const auto& bt1 = m_stg_ptr->bt1_medialCurve;
	for (unsigned i = 0; i < m_stg_ptr->dual_vts.size(); ++i)
	{
		float val = bt1[i] - bt2[i];
		max_val = std::max(max_val, val);
		min_val = std::min(min_val, val);
	}

	_eps = _smooth_ratio * (max_val - min_val) + min_val;

	/*float max_bt2 = numeric_limits<float>::min();
	float min_bt2 = numeric_limits<float>::max();
	const auto& bt2 = m_stg_ptr->bt2MA_vert;
	for (unsigned i = 0; i < bt2.size(); ++i)
	{
		max_bt2 = std::max(max_bt2, bt2[i]);
		min_bt2 = std::min(min_bt2, bt2[i]);
	}

	_eps = _smooth_ratio * (max_bt2 - min_bt2);*/
}

void SurfaceFunc::convert_to_shape_width_extremity_eps_range(const float _smooth_ratio, float& _eps)
{
	float max_bt1 = numeric_limits<float>::min();
	float min_bt1 = numeric_limits<float>::max();

	const auto& bt1 = m_stg_ptr->bt1_medialCurve;
	for (unsigned i = 0; i < m_stg_ptr->dual_vts.size(); ++i)
	{
		max_bt1 = std::max(max_bt1, bt1[i]);
		min_bt1 = std::min(min_bt1, bt1[i]);
	}

	_eps = _smooth_ratio * (max_bt1 - min_bt1);
}

void SurfaceFunc::get_shape_diam(vector<float>& _scalar_vals, float _eps)
{
	const auto& bt3_mc = m_stg_ptr->bt3_MC;
	obtain_smooth_field_using_mc(SHAPE_DIAM, bt3_mc, _scalar_vals, _eps);

	//cout << "shape diameter obtained." << endl;
}

void SurfaceFunc::get_shape_width(vector<float>& _scalar_vals, float _eps)
{
	/*
	const auto& bt2_ma = m_stg_ptr->bt2MA_vert;

	// smooth the bt2 field on the orig surface using bt2 gradient
	// estimated by the one-ring neighborhood
	_scalar_vals.assign(m_closest_vert_for_origSurf.size(), 0.0f);
	for (unsigned i = 0; i < m_closest_vert_for_origSurf.size(); i ++)
	{
		auto vi_on_ma = m_closest_vert_for_origSurf[i];
		if (m_closest_vert_for_origSurf[i] == -1)
			continue;

		int local_minimum_vi = this->find_local_maximum(vi_on_ma, _eps);
		if ( local_minimum_vi != -1 )
		{
			_scalar_vals[i] = bt2_ma[local_minimum_vi];
		}
	}
	*/

	const auto& bt2_mc = m_stg_ptr->bt2_MC;
	obtain_smooth_field_using_mc(SHAPE_WIDTH, bt2_mc, _scalar_vals, _eps);

	//cout << "shape width obtained." << endl;
}

void SurfaceFunc::get_shape_width_extremity(vector<float>& _scalar_vals, float _eps)
{
	const auto& bt2_mc = m_stg_ptr->bt2_MC;
	obtain_smooth_field_using_mc(SHAPE_WIDTH_EXTREMITY, bt2_mc, _scalar_vals, _eps);
}

void SurfaceFunc::obtain_smooth_field_using_mc(
	SurfFuncType _type, 
	const vector<float>& _mc_field, 
	vector<float>& _scalar_vals, 
	float _eps
	)
{
	const auto& mc_next = m_stg_ptr->burnNext_medialCurve;
	vector<float> cur_final_quantity_map; // on MC: vi -> the final vi/quantity vi ends up with
	cur_final_quantity_map.assign(m_stg_ptr->dual_vts.size(), -1.0f);

	// smooth the bt2 field on the orig surface using difference btw tree-dist and bt2
	// as we travel on the MC, from the MC vert that the orig vert mapped to
	_scalar_vals.assign(m_origSurf_graph->vts.size(), INVALID_VAL);
	int final_vi;
	float final_treeDist;
	for ( unsigned i = 0; i < this->m_close_vts_for_origSurf.size(); ++i )
	{
		if (i%5000 == 0)
			cout << "smoothing done for "<<i<<" vertices."<<endl;

		const auto& close_vts = m_close_vts_for_origSurf[i];
		if ( close_vts.empty() )
			continue;

		float weighted_scalar = 0.0f;
		float w_sum = 0.0f;
		for ( unsigned j = 0; j < close_vts.size(); ++j )
		{
			const auto& v_w_pr = close_vts[j];
			float w = v_w_pr.second;
			w_sum += w;
			unsigned cur_vi = v_w_pr.first;
			unsigned pre_vi = cur_vi;
			float tree_dist;

			if ( _type == SHAPE_WIDTH || _type == SHAPE_DIAM )
			{
				final_vi = cur_final_quantity_map[cur_vi];
				if (final_vi < 0)
				{
					tree_dist = _mc_field[cur_vi];
					while ( tree_dist - _mc_field[cur_vi] < _eps && mc_next[cur_vi] != -1 )
					{
						// goto next vert & increment tree_dist by edge len <cur_vi, next_vi>
						tree_dist += trimesh::dist(
							m_stg_ptr->dual_vts[cur_vi], 
							m_stg_ptr->dual_vts[mc_next[cur_vi]]
						);
						tree_dist = std::max(tree_dist, _mc_field[mc_next[cur_vi]]);
						pre_vi = cur_vi;
						cur_vi = mc_next[cur_vi];
					}
					final_vi = pre_vi;
					cur_final_quantity_map[v_w_pr.first] = final_vi;
				}

				if ( _type == SHAPE_WIDTH )
				{
					weighted_scalar += w * m_stg_ptr->bt2_MC[final_vi];
				}
				else if ( _type == SHAPE_DIAM )
				{
					weighted_scalar += w * m_stg_ptr->bt3_MC[final_vi];
				}
			}
			else if ( _type == SHAPE_WIDTH_EXTREMITY || _type == SHAPE_LENGTH_EXTREMITY )
			{
				final_treeDist = cur_final_quantity_map[cur_vi];
				if (final_treeDist < 0.0f)
				{
					tree_dist = 0.0f;
					float tree_dist_2 = _mc_field[cur_vi];
					while (tree_dist_2 - _mc_field[cur_vi] < _eps && mc_next[cur_vi] != -1 )
					{
						//goto next vert & increment tree_dist by edge len <cur_vi, next_vi>
						float e_len = trimesh::dist(
							m_stg_ptr->dual_vts[cur_vi], m_stg_ptr->dual_vts[mc_next[cur_vi]]);
						tree_dist += e_len;
						tree_dist_2 += e_len;
						tree_dist_2 = std::max(tree_dist_2, _mc_field[mc_next[cur_vi]]);
						pre_vi = cur_vi;
						cur_vi = mc_next[cur_vi];
					}
					final_treeDist = tree_dist;

					cur_final_quantity_map[v_w_pr.first] = final_treeDist;
				}

				weighted_scalar += w * final_treeDist;
			}
			else
			{
				cout << "Error: couldn't recognize SurfFuncType." <<endl;
				exit(-1);
			}
		}
		weighted_scalar /= w_sum;
		_scalar_vals[i] = weighted_scalar;
	}

	cur_final_quantity_map.clear();
}