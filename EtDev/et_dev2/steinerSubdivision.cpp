#include "steinerSubdivision.h"
#include "steinerSubdivision_inline.h"
#include "origGraph.h"
#include "commonUtils.h"

SteinerSubdivision::SteinerSubdivision()
{
	this->clear();
}

SteinerSubdivision::~SteinerSubdivision()
{

}

/***************************/	
/***** public interfaces *****/	
/***************************/	
bool SteinerSubdivision::fixedSubdivide(const shared_ptr<MyGraph>& _g_ptr, int _nSamples)
{
	this->clear();
	m_origG_ptr = _g_ptr;
	this->m_subd_scheme = FIXED;

	if (_nSamples <= 0)
		_nSamples = 0;
	m_nFixed = _nSamples;

	// increase space to accommodate new steiner vts
	m_stVts = m_origG_ptr->vts;
	m_stVts.reserve(m_stVts.size() + _nSamples*m_origG_ptr->edges.size());
	m_triEdge_stVtsRange_list.resize(m_origG_ptr->edges.size());

	// generate new steiner vts on each edge
	// and connect each steiner vert to its neighbor on this edge
	float interp_weight;
	vector<unsigned> stVtsOnE;
	for (unsigned ei = 0; ei < m_origG_ptr->edges.size(); ++ei)
	{
		TriEdge orig_e = m_origG_ptr->edges[ei];
		auto& pr = m_triEdge_stVts_map[orig_e]; // TODO: OPTIMIZE_AWAY_THIS
		m_triEdge_stVtsRange_list[ei] = std::make_pair(m_stVts.size(), m_nFixed);
		pr = m_triEdge_stVtsRange_list[ei];
		stVtsOnE.clear();

		for (int i = 0; i < _nSamples; ++i)
		{
			interp_weight = float(i+1) / (_nSamples+1);
			m_stVts.push_back(
				trimesh::mix(this->m_stVts[orig_e[0]], m_stVts[orig_e[1]], interp_weight)
				);
			stVtsOnE.push_back(m_stVts.size() - 1);
			m_stV_triEdge_map.push_back(make_pair(ei, i));
		}

		if (stVtsOnE.empty())
			m_stEdges.push_back(orig_e);
		//else
		//{
		//	// v0 to steiner vert 1
		//	m_stEdges.push_back(util::makeEdge(orig_e[0], stVtsOnE[0]));

		//	// last steiner vert to v1
		//	m_stEdges.push_back(util::makeEdge(stVtsOnE[stVtsOnE.size() - 1], orig_e[1]));

		//	// steiners in between v0 & v1
		//	if ( stVtsOnE.size() > 1 )
		//		for (unsigned st_i = 0; st_i < stVtsOnE.size()-1; ++st_i)
		//		{
		//			m_stEdges.push_back(util::makeEdge(stVtsOnE[st_i], stVtsOnE[st_i+1]));
		//		}
		//}
	}

	set<TriEdge> visited_origedges;
	for (unsigned fi = 0; fi < m_origG_ptr->faces.size(); ++fi)
	{
		const auto& f = m_origG_ptr->faces[fi];
		connect_st_edges_of_face(f, visited_origedges, this->m_stEdges);
	}
	visited_origedges.clear();

	build_st_edge_mapping();

	return true;
}

bool SteinerSubdivision::adaptiveSubdivide(const shared_ptr<MyGraph>& _g_ptr, const double _param)
{
	this->clear();
	m_origG_ptr = _g_ptr;
	this->m_subd_scheme = ADAPTIVE;
	
	// compute average edge length of graph and put _d samples on that avg edge
	/*float avg_len = 0.0f;
	for (unsigned ei = 0; ei < m_origG_ptr->edges.size(); ++ei)
	{
		const auto& e = m_origG_ptr->edges[ei];
		avg_len += trimesh::dist(m_origG_ptr->vts[e[0]], m_origG_ptr->vts[e[1]]);
	}
	this->m_d = avg_len * _param;*/
	//this->m_d = (avg_len / m_origG_ptr->edges.size()) / (std::max(_d, 0) + 1);
	this->m_d = _param;

	m_stVts = m_origG_ptr->vts;
	m_triEdge_stVtsRange_list.resize( 
		m_origG_ptr->edges.size(), 
		pair<unsigned, unsigned>( 0, 0 ) ); // no steiner vts by default

	// generate new steiner vts on each edge
	// and connect each steiner vert to its neighbor on this edge
	float interp_weight;
	vector<unsigned> stVtsOnE;
	int n_samples;
	for (unsigned ei = 0; ei < m_origG_ptr->edges.size(); ++ei)
	{
		// skip isolated line geometry (used by no face)
		if ( m_origG_ptr->getNbFaces( ei ).size() == 0 )
			continue;
		// determine how many samples on this edge. At least there will be 1.
		TriEdge orig_e = m_origG_ptr->edges[ei];
		n_samples = (int)(
			std::max(
				trimesh::dist(
				this->getVert(orig_e[0]), this->getVert(orig_e[1])
				) / m_d - 1, 1.0f
			) + 0.5f
			);

		// save the range of the st vts on this edge
		auto& pr = m_triEdge_stVts_map[orig_e]; // TODO: OPTIMIZE_AWAY_THIS
		m_triEdge_stVtsRange_list[ei] = std::make_pair(m_stVts.size(), n_samples);
		pr = m_triEdge_stVtsRange_list[ei];
		stVtsOnE.clear();

		for (int i = 0; i < n_samples; ++i)
		{
			interp_weight = float(i+1) / (n_samples+1);
			m_stVts.push_back(
				trimesh::mix(this->m_stVts[orig_e[0]], m_stVts[orig_e[1]], interp_weight)
				);
			stVtsOnE.push_back(m_stVts.size() - 1);
			m_stV_triEdge_map.push_back(make_pair(ei, i));
		}
	}

	/*set<TriEdge> visited_origedges;
	for (unsigned fi = 0; fi < m_origG_ptr->faces.size(); ++fi)
	{
		const auto& f = m_origG_ptr->faces[fi];
		connect_st_edges_of_face(f, visited_origedges, this->m_stEdges);
	}
	visited_origedges.clear();*/

	return true;
}

bool SteinerSubdivision::adaptiveMidPointSubdivide(const shared_ptr<MyGraph>& _g_ptr, const double _param)
{
	this->clear();
	m_origG_ptr = _g_ptr;
	this->m_subd_scheme = ADAPTIVE_MIDPOINT;

	// compute average edge length of graph and put _d samples on that avg edge
	/*float avg_len = 0.0f;
	for (unsigned ei = 0; ei < m_origG_ptr->edges.size(); ++ei)
	{
		const auto& e = m_origG_ptr->edges[ei];
		avg_len += trimesh::dist(m_origG_ptr->vts[e[0]], m_origG_ptr->vts[e[1]]);
	}
	this->m_d = avg_len * _param;*/
	//this->m_d = (avg_len / m_origG_ptr->edges.size()) / (1 << _l);
	this->m_d = _param;

	m_stVts = m_origG_ptr->vts;
	m_triEdge_stVtsRange_list.resize(m_origG_ptr->edges.size());

	// generate new steiner vts on each edge
	// and connect each steiner vert to its neighbor on this edge
	float interp_weight;
	vector<unsigned> stVtsOnE;
	int n_samples;
	for (unsigned ei = 0; ei < m_origG_ptr->edges.size(); ++ei)
	{
		// determine how many samples on this edge. At least there will be 1.
		TriEdge orig_e = m_origG_ptr->edges[ei];
		n_samples = (int)(
			std::max(
			trimesh::dist(
			this->getVert(orig_e[0]), this->getVert(orig_e[1])
			) / m_d - 1, 1.0f
			) + 0.5f
			);

		// DIFFERENCE between this and MIDPOINT:
		// find the smallest power of 2 that's greater than n_samples
		// Such value will be the number of intervals. n_samples will take such value minus 1.
		int n = 0;
		int n_intervals = n_samples + 1;
		while ( (n_intervals >> n) > 0 ) n++;
		n_intervals = 1 << (n-1);
		n_samples = n_intervals - 1;

		// save the range of the st vts on this edge
		auto& pr = m_triEdge_stVts_map[orig_e]; // TODO: OPTIMIZE_AWAY_THIS
		m_triEdge_stVtsRange_list[ei] = std::make_pair(m_stVts.size(), n_samples);
		pr = m_triEdge_stVtsRange_list[ei];
		stVtsOnE.clear();

		for (int i = 0; i < n_samples; ++i)
		{
			interp_weight = float(i+1) / (n_samples+1);
			m_stVts.push_back(
				trimesh::mix(this->m_stVts[orig_e[0]], m_stVts[orig_e[1]], interp_weight)
				);
			stVtsOnE.push_back(m_stVts.size() - 1);
			m_stV_triEdge_map.push_back(make_pair(ei, i));
		}

		if (stVtsOnE.empty())
			m_stEdges.push_back(orig_e);
	}

	set<TriEdge> visited_origedges;
	for (unsigned fi = 0; fi < m_origG_ptr->faces.size(); ++fi)
	{
		const auto& f = m_origG_ptr->faces[fi];
		connect_st_edges_of_face(f, visited_origedges, this->m_stEdges);
	}
	visited_origedges.clear();

	build_st_edge_mapping();

	return true;
}

unsigned SteinerSubdivision::sizeOfStVts() const
{
	return sizeOfVts() - m_origG_ptr->vts.size();
}

int SteinerSubdivision::getClosestVertOnEdge(unsigned _u, unsigned _v, unsigned _which) const
{
	unsigned u, v;
	if (_which == 0)
	{
		u = _u;
		v = _v;
	}
	else
	{
		u = _v;
		v = _u;
	}

	const auto& stOnE_range = m_triEdge_stVts_map.find(util::makeEdge(u, v))->second;
	if (stOnE_range.second == 0)
		return v;
	else
		return u < v ? stOnE_range.first : stOnE_range.first+stOnE_range.second-1;
}

int SteinerSubdivision::getResidingFaceIdx(const TriEdge _e) const
{
	assert(isSteinerVert(_e[0]) || isSteinerVert(_e[1]));

	TriEdge e1;
	TriEdge e2;

	unsigned u;
	unsigned v;
	unsigned w;

	if (!isSteinerVert(_e[0]))
	{
		w = _e[0];
		e1 = getResidingEdge(_e[1]);
		u = e1[0]; 
		v = e1[1];
	}
	else if (!isSteinerVert(_e[1]))
	{
		w = _e[1];
		e1 = getResidingEdge(_e[0]);
		u = e1[0]; 
		v = e1[1];
	}
	else // both are st vts
	{
		e1 = getResidingEdge(_e[1]);
		e2 = getResidingEdge(_e[0]);
		u = e1[0];
		v = e1[1];
		w = e2[0] == u || e2[0] == v ? e2[1] : e2[0];
	}

	auto find_it = m_origG_ptr->triFace_faceIdx_map.find(util::makeFace(u, v, w));
	return find_it == m_origG_ptr->triFace_faceIdx_map.end() ? -1 : find_it->second;
}

int SteinerSubdivision::getNbVtsInFace(
	const unsigned _vi, const unsigned _fi, bool _exclude_orig, vector<int>& _nb_vts ) const
{
	vector<unsigned> stVts_onE;
	int nb_v1;
	int nb_v2;
	if (isSteinerVert(_vi))
	{
		unsigned cur_ei;
		int idx_withinE;
		getResidingEdgeIdx(_vi, cur_ei, idx_withinE);
		TriEdge cur_e = m_origG_ptr->edges[cur_ei];
		stVts_onE.clear();
		getStVertIndicesOnTriEdge(cur_ei, stVts_onE);
		/*if (this->m_origG_ptr->isManifoldEdge(cur_e))
		{
			nb_v1 = idx_withinE == 0 ? cur_e[0] : stVts_onE[idx_withinE-1];
			nb_v2 = idx_withinE == stVts_onE.size()-1 ? cur_e[1] : stVts_onE[idx_withinE+1];
			nb_vts.push_back(nb_v1);
			nb_vts.push_back(nb_v2);
		}*/
		
		if ( _exclude_orig ) // do not include adjacent orig vts
		{
			if ( idx_withinE > 0 )
			{
				nb_v1 = stVts_onE[ idx_withinE - 1 ];
				_nb_vts.push_back( nb_v1 );
			}
			if ( idx_withinE < stVts_onE.size() - 1 )
			{
				nb_v2 = stVts_onE[ idx_withinE + 1 ];
				_nb_vts.push_back( stVts_onE[ idx_withinE + 1 ] );
			}
		}
		else // need adjacent orig vts
		{
			nb_v1 = idx_withinE == 0 ? cur_e[ 0 ] : stVts_onE[ idx_withinE - 1 ];
			nb_v2 = idx_withinE == stVts_onE.size() - 1 ? cur_e[ 1 ] : stVts_onE[ idx_withinE + 1 ];
			_nb_vts.push_back( nb_v1 );
			_nb_vts.push_back( nb_v2 );
		}

		unsigned oppo_vi = util::oppositeVert(m_origG_ptr->faces[_fi], cur_e);
		if ( !_exclude_orig ) // need the oppo orig vert as well.
			_nb_vts.push_back(oppo_vi);

		const auto& oppo_e1 = util::makeEdge(oppo_vi, cur_e[0]);
		stVts_onE.clear();
		this->getStVertIndicesOnTriEdge(oppo_e1, stVts_onE);
		_nb_vts.insert(_nb_vts.end(), stVts_onE.begin(), stVts_onE.end());

		const auto& oppo_e2 = util::makeEdge(oppo_vi, cur_e[1]);
		stVts_onE.clear();
		this->getStVertIndicesOnTriEdge(oppo_e2, stVts_onE);
		_nb_vts.insert(_nb_vts.end(), stVts_onE.begin(), stVts_onE.end());
	}
	else // vi is an orig tri vert
	{
		TriEdge oppo_e = util::oppositeEdge(m_origG_ptr->faces[_fi], _vi);
		/*if ( this->m_origG_ptr->isManifoldEdge(util::makeEdge(oppo_e[0], _vi)) )
		{
			nb_v1 = getClosestVertOnEdge(oppo_e[0], _vi, 1);
			nb_vts.push_back(nb_v1);
		}
		if ( this->m_origG_ptr->isManifoldEdge(util::makeEdge(oppo_e[1], _vi)) )
		{
			nb_v2 = getClosestVertOnEdge(oppo_e[1], _vi, 1);
			nb_vts.push_back(nb_v2);
		}*/
		nb_v1 = getClosestVertOnEdge(oppo_e[0], _vi, 1);
		_nb_vts.push_back(nb_v1);
		nb_v2 = getClosestVertOnEdge(oppo_e[1], _vi, 1);
		_nb_vts.push_back(nb_v2);

		stVts_onE.clear();
		getStVertIndicesOnTriEdge(oppo_e, stVts_onE);
		_nb_vts.insert(_nb_vts.end(), stVts_onE.begin(), stVts_onE.end());
	}

	return 1;
}

void SteinerSubdivision::getStEdgesOfEdge(const unsigned _ei, vector<TriEdge>& _st_edges_on_e) const
{
	vector<unsigned> all_vts_on_e;
	getVtsOnTriEdge(_ei, all_vts_on_e);
	for (size_t i = 0; i < all_vts_on_e.size() - 1; ++i)
	{
		_st_edges_on_e.push_back(util::makeEdge(all_vts_on_e[i],all_vts_on_e[i+1]));
	}
}

void SteinerSubdivision::getStEdgesWithinFace(const unsigned _fi, vector<TriEdge>& _st_edges_within_f) const
{
	vector<unsigned> st_vts_on_edges[3]; // st vts on each edge
	int oppo_vi[3]; // orig v oppo. to an edge
	TriEdge e;
	const auto& f = m_origG_ptr->faces[_fi];
	for (unsigned i = 0; i < 3; ++i)
	{
		e = util::makeEdge(f[i], f[(i+1)%3]);
		int ei = m_origG_ptr->getEdgeIdx(e);
		getStVertIndicesOnTriEdge(ei, st_vts_on_edges[i]);
		oppo_vi[i] = f[(i+2)%3];
	}

	for (unsigned i = 0; i < 3; ++i)
	{
		// connect st edges from one edge to the next 
		const auto& stvts_on_cur_e = st_vts_on_edges[i];
		const auto& stvts_on_next_e = st_vts_on_edges[(i+1)%3];
		for (auto cur_it = stvts_on_cur_e.begin(); cur_it != stvts_on_cur_e.end(); ++cur_it)
		{
			for (auto next_it = stvts_on_next_e.begin(); next_it != stvts_on_next_e.end(); ++next_it)
			{
				_st_edges_within_f.push_back( util::makeEdge(*cur_it, *next_it) );
			}
		}
		
		if (CurrentConnectScheme() == CONNECT_TO_ORIG_VTS)
		{
			// connect st edges from one edge to its oppo orig vert
			for (auto cur_it = stvts_on_cur_e.begin(); cur_it != stvts_on_cur_e.end(); ++cur_it)
			{
				_st_edges_within_f.push_back( util::makeEdge(*cur_it, oppo_vi[i]) );
			}
		}
	}
}

void SteinerSubdivision::getStEdgesOfFace(const unsigned _fi, vector<TriEdge>& _st_edges_of_f) const
{
	const auto& f = m_origG_ptr->faces[_fi];
	// append st edges on each edge to output list
	TriEdge e;
	for (unsigned i = 0; i < 3; ++i)
	{
		e = util::makeEdge(f[i], f[(i+1)%3]);
		getStEdgesOfEdge(m_origG_ptr->getEdgeIdx(e), _st_edges_of_f);
	}
	// append st edges lying within the face to output list
	getStEdgesWithinFace(_fi, _st_edges_of_f);
}

// TODO: finish this function.
void SteinerSubdivision::getNbVerts(const unsigned _vi, vector<unsigned>& _nbvts)
{
	vector<unsigned> stVts_onE;
	unsigned nb_v1;
	unsigned nb_v2;
	_nbvts.clear();
	if (isSteinerVert(_vi))
	{
		unsigned residing_ei;
		int idx_withinE;
		getResidingEdgeIdx(_vi, residing_ei, idx_withinE);
		TriEdge residing_e = getOrigTriEdge(residing_ei);
		const auto& _f_indices = m_origG_ptr->getNbFaces( residing_e );

		stVts_onE.clear();
		getStVertIndicesOnTriEdge(residing_ei, stVts_onE);
		nb_v1 = idx_withinE == 0 ? residing_e[0] : stVts_onE[idx_withinE-1];
		nb_v2 = idx_withinE == stVts_onE.size()-1 ? residing_e[1] : stVts_onE[idx_withinE+1];
		_nbvts.push_back(nb_v1);
		_nbvts.push_back(nb_v2);

		for (auto fi_it = _f_indices.begin(); fi_it != _f_indices.end(); ++fi_it)
		{
			unsigned oppo_vi = util::oppositeVert(m_origG_ptr->faces[*fi_it], residing_e);
			_nbvts.push_back(oppo_vi);

			const auto& oppo_e1 = util::makeEdge(oppo_vi, residing_e[0]);
			stVts_onE.clear();
			this->getStVertIndicesOnTriEdge(oppo_e1, stVts_onE);
			_nbvts.insert(_nbvts.end(), stVts_onE.begin(), stVts_onE.end());

			const auto& oppo_e2 = util::makeEdge(oppo_vi, residing_e[1]);
			stVts_onE.clear();
			this->getStVertIndicesOnTriEdge(oppo_e2, stVts_onE);
			_nbvts.insert(_nbvts.end(), stVts_onE.begin(), stVts_onE.end());
		}
	}
	else // vi is an orig tri vert
	{
		const auto& _f_indices = m_origG_ptr->nbFacesOfVert[_vi];
		for (auto fi_it = _f_indices.begin(); fi_it != _f_indices.end(); ++fi_it)
		{
			TriEdge oppo_e = util::oppositeEdge(m_origG_ptr->faces[*fi_it], _vi);
			nb_v1 = getClosestVertOnEdge(oppo_e[0], _vi, 1);
			nb_v2 = getClosestVertOnEdge(oppo_e[1], _vi, 1);
			_nbvts.push_back(nb_v1);
			_nbvts.push_back(nb_v2);

			stVts_onE.clear();
			getStVertIndicesOnTriEdge(oppo_e, stVts_onE);
			_nbvts.insert(_nbvts.end(), stVts_onE.begin(), stVts_onE.end());
		}
	}
}

/***************************/	
/***** private helpers *****/	
/***************************/	
void SteinerSubdivision::clear()
{
	this->m_n_stEdges = 0;
	this->m_stEdges.clear();
	this->m_stVts.clear();
	this->m_triEdge_stVtsRange_list.clear();
	this->m_stV_triEdge_map.clear();

	this->m_nFixed = 0;
	this->m_subd_scheme = FIXED;
	this->m_conn_scheme = CONNECT_TO_ORIG_VTS;
}

void SteinerSubdivision::connect_st_edges_of_face(
	const TriFace _f, 
	set<TriEdge>& _visited_origedges, 
	vector<TriEdge>& _edges
	)
{
	const auto& edges_f = util::edgesFromFace(_f);
	vector<unsigned> stVtsOnE[3];
	getStVertIndicesOnTriEdge(edges_f[0], stVtsOnE[0]);
	getStVertIndicesOnTriEdge(edges_f[1], stVtsOnE[1]);
	getStVertIndicesOnTriEdge(edges_f[2], stVtsOnE[2]);

	// 1. connect st edges along each orig tri edge
	for (unsigned ei = 0; ei < 3; ++ei)
	{
		if ( _visited_origedges.count(edges_f[ei]) )
			continue;

		_visited_origedges.insert(edges_f[ei]);	
		if (!stVtsOnE[ei].empty())
		{
			// v0 to steiner vert 1
			_edges.push_back(util::makeEdge(edges_f[ei][0], stVtsOnE[ei][0]));

			// last steiner vert to v1
			_edges.push_back(util::makeEdge(stVtsOnE[ei].back(), edges_f[ei][1]));

			// steiners in between v0 & v1
			if ( stVtsOnE[ei].size() > 1 )
				for (unsigned st_i = 0; st_i < stVtsOnE[ei].size()-1; ++st_i)
				{
					_edges.push_back(util::makeEdge(stVtsOnE[ei][st_i], stVtsOnE[ei][st_i+1]));
				}
		}
	}

	// 2. connect st edges between each pair of edges
	for (unsigned ei = 0; ei < 3; ++ei)
	{
		const auto& cur_stVts = stVtsOnE[ei];
		const auto& next_stVts = stVtsOnE[(ei+1)%3];

		unsigned oppo_v = util::oppositeVert(_f, edges_f[ei]);
		for (auto u_it = cur_stVts.begin(); u_it != cur_stVts.end(); ++u_it)
		{
			for (auto v_it = next_stVts.begin(); v_it != next_stVts.end(); ++v_it)
			{
				_edges.push_back(util::makeEdge(*u_it, *v_it));
			}
			_edges.push_back(util::makeEdge(*u_it, oppo_v));
		}
	}
}

void SteinerSubdivision::build_st_edge_mapping()
{
	for (unsigned ei = 0; ei < m_stEdges.size(); ++ei)
		m_stEdge_idx_map[m_stEdges[ei]] = ei;
}