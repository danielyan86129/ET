#include <iostream>

#include "graphDiffusion.h"

using std::cout;
using std::endl;

/*
** Public interfaces 
**
*/
GraphDiffuser::GraphDiffuser()
{
	
}

GraphDiffuser::~GraphDiffuser()
{

}

GraphDiffuser::Errcode GraphDiffuser::diffuseOnGraph( 
	const MyGraph* const _g, 
	GraphWeightFunc& _wFunc, 
	const vector<bool>& _is_known, 
	vector<float>& _scalars )
{
	Errcode err_code;
	err_code = setupDiffuseSystem( _g->vts, _g->nbVtsOfVert, _wFunc, _is_known );
	if (err_code == FAILED)
		return FAILED;
	err_code = solve(_scalars);
	return err_code;
}

GraphDiffuser::Errcode GraphDiffuser::setupDiffuseSystem(
	const MyGraph* const _g,
	GraphWeightFunc& _wFunc, /*the weighting functor*/
	const vector<bool>& _is_known /*the indicating list*/
	)
{
	return setupDiffuseSystem(
		_g->vts, _g->nbVtsOfVert, _wFunc, _is_known
		);
}

/// the system: A.x = b, where b = C.f, and f 
/// is only determined in run time (the known scalars)
GraphDiffuser::Errcode GraphDiffuser::setupDiffuseSystem(
	const vector<TriPoint>& _vts,
	const vector<vector<int>>& _vts_adjacency,
	GraphWeightFunc& _wFunc, /*the weighting functor*/
	const vector<bool>& _is_known 
	)
{
	// for debug purpose (remove for optimization)
	this->m_gVts = _vts;
	set<TriEdge> edge_union;
	for (auto orig_vi = 0; orig_vi < _vts_adjacency.size(); ++orig_vi)
	{
		const auto& nbs = _vts_adjacency[orig_vi];
		for (auto vj = 0; vj < nbs.size(); ++vj)
		{
			edge_union.insert( util::makeEdge(orig_vi, nbs[vj]) );
		}
	}
	this->m_gEdges.clear();
	m_gEdges.insert(m_gEdges.end(), edge_union.begin(), edge_union.end());
	edge_union.clear();

	this->m_isKnown = _is_known;
	vector<SpTriplet> A,C;
	A.clear();

	// original vert id -> new id in known/unknown domain (as known and unknowns will be separated)
	vector<unsigned> orig_to_knownOrUnknown_Idx_map(_vts_adjacency.size(), 0);
	// separate nodes with and w/o scalar values
	int unknown_vts_cnt = 0, known_vts_cnt = 0;
	int isolated_cnt = 0;
	m_new_orig_vertIdx_map.clear();
	for ( unsigned orig_vi = 0; orig_vi < _vts_adjacency.size(); ++orig_vi )
	{
		if ( _vts_adjacency[orig_vi].empty() )
		{
			// this node has no incident edges. 
			// don't include it in the graph.
			isolated_cnt ++;
			continue;
		}

		unsigned new_vi;
		if ( !_is_known[orig_vi] )
		{
			new_vi = unknown_vts_cnt;
			unknown_vts_cnt ++;
		}
		else
		{
			new_vi = known_vts_cnt;
			known_vts_cnt ++;
		}
		orig_to_knownOrUnknown_Idx_map[orig_vi] = new_vi;
		m_new_orig_vertIdx_map.push_back(orig_vi);
	}

	cout << "# isolated nodes skipped: "<<isolated_cnt << endl;
	if ( unknown_vts_cnt == 0 ) // there is no unknown vars.
		return ZERO_UNKNOWNS;

	cout << "# knowns/unknowns: "<<known_vts_cnt<<"/"<<unknown_vts_cnt<<endl;

	// construct A & C
	float wij_max = 0.0f;
	float wij_min = numeric_limits<float>::max();
	for ( unsigned new_vi = 0; new_vi < m_new_orig_vertIdx_map.size(); ++new_vi )
	{
		auto orig_vi = m_new_orig_vertIdx_map[new_vi];
		if ( !_is_known[orig_vi] )
		{
			auto unknown_i = orig_to_knownOrUnknown_Idx_map[orig_vi];
			float sum_w = 0.0f; // weight for wii
			const auto& nbs = _vts_adjacency[orig_vi];
			for ( auto it = nbs.begin(); it != nbs.end(); ++it )
			{
				// compute the weight for edge between vi and this nb
				float wij = _wFunc(orig_vi, *it);
				if( ::is_valid(wij) == false )
					cout << "Invalid number! w("<<orig_vi<<","<<*it<<")="<<wij<<", ";
				sum_w += wij;
				wij_min = std::min(wij_min, wij);
				wij_max = std::max(wij_max, wij);

				// put the weight for nb in the proper sparse matrix (A or C)
				if ( !_is_known[*it] )
				{
					auto unknown_j = orig_to_knownOrUnknown_Idx_map[*it];
					A.push_back( SpTriplet(unknown_i, unknown_j, -wij) );
				}
				else
				{
					auto known_j = orig_to_knownOrUnknown_Idx_map[*it];
					C.push_back( SpTriplet(unknown_i, known_j, wij) );
				}
			}

			// put wii into _A
			A.push_back( SpTriplet(unknown_i, unknown_i, sum_w) );
			if (::almost_zero(sum_w, 0.00001f) )
			{
				//cout << "orig_vi, new_vi: "<<orig_vi<<","<<new_vi<<endl;
				cout << "A: zero diagonal identified! ("
					<< unknown_i<<","<<unknown_i<<")="<<sum_w<<endl;
			}
		}
	}
	cout << endl;
	cout <<"wij range: "<<wij_min<<"/"<<wij_max<<endl;

	f.resize(known_vts_cnt);
	A_sm.resize(unknown_vts_cnt, unknown_vts_cnt);
	A_sm.setFromTriplets( A.begin(), A.end() );
	C_sm.resize(unknown_vts_cnt, known_vts_cnt);
	C_sm.setFromTriplets( C.begin(), C.end() );

	//this->outputMatA();

	/* sanity check */
	// A_sm must have nonzero diagonal
	auto diag_A = A_sm.diagonal();
	if ( diag_A.isZero(0.00001) )
	{
		cout << "A_sm diagonal contains zeros!"<<endl;
	}

	this->decomp_of_A.compute(A_sm);
	auto ret = decomp_of_A.info();
	if (ret == Eigen::Success)
	{
		cout << "decomposition success."<<endl;
		return SUCCESS;
	}
	else if (ret == Eigen::NumericalIssue)
	{
		cout << "decomposition numerical issue!"<<endl;
		return FAILED;
	}

	return SUCCESS;
}

GraphDiffuser::Errcode GraphDiffuser::solve(
	vector<float>& _scalar_vals 
	)
{
	// known scalars
	vector<float> known_scalars;
	for (unsigned new_vi = 0; new_vi < m_new_orig_vertIdx_map.size(); new_vi++)
	{
		auto orig_i = m_new_orig_vertIdx_map[new_vi];
		if ( m_isKnown[orig_i] )
		{
			known_scalars.push_back( _scalar_vals[orig_i] );
		}
	}

	if (known_scalars.empty())
	{
		cout << "Cannot perform diffusion: there is no known scalar."<<endl;
		return FAILED;
	}
	if (known_scalars.size() == m_new_orig_vertIdx_map.size())
	{
		cout << "Zero unknowns." << endl;
		return ZERO_UNKNOWNS;
	}

	// solve the diffusion system
	vector<float> unknown_solved;
	solve_diffusion_system(known_scalars, unknown_solved);

	int unknown_i = 0;
	float val_max = numeric_limits<float>::min();
	float val_min = numeric_limits<float>::max();
	for (unsigned new_vi = 0; new_vi < m_new_orig_vertIdx_map.size(); new_vi ++)
	{
		auto orig_i = m_new_orig_vertIdx_map[new_vi];
		if ( !m_isKnown[orig_i] )
		{
			auto val = unknown_solved[unknown_i];
			_scalar_vals[orig_i] = val;
			unknown_i ++;
			val_max = std::max(val_max, val);
			val_min = std::min(val_min, val);
		}
	}

	cout << "unknown scalars (solved) range: "<< val_min<<"/"<<val_max<<endl;
	return SUCCESS;
}

void GraphDiffuser::outputDiffuseSystem( 
	const vector<float>& _scalars 
	)
{
	ofstream o_file("diffuse_on_graph.txt"); 
	std::stringstream out_string;
	out_string << "{"<< endl; // Big {

	// Begin: output graph vts
	out_string << "{ (*Begin: vts*)"<< endl; 
	for (unsigned i = 0; i < m_gVts.size(); ++i)
	{
		const auto & p = m_gVts[i];
		out_string << p[0]<<","<<p[1]<<","<<p[2];
		if (i != m_gVts.size()-1)
		{
			out_string <<",";
		}
	}
	out_string << endl;
	out_string << "}, (*End: vts*)"<< endl; 
	// End: output graph vts

	// Begin: output graph edges
	out_string << "{ (*Begin: graph edges*)"<< endl; 
	for (unsigned i = 0; i < m_gEdges.size(); ++i)
	{
		const auto & e = m_gEdges[i];
		out_string << e[0]+1<<","<<e[1]+1;
		if (i != m_gEdges.size()-1)
		{
			out_string <<",";
		}
	}
	out_string << endl;
	out_string << "}, (*End: graph edges*)"<< endl; 
	// End: output graph faces

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
	for ( int  i = 0; i < m_isKnown.size(); ++i )
	{
		out_string << (m_isKnown[i] ? 1 : 0) << ",";
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

void GraphDiffuser::outputMatA(
	//const vector<float>& _scalars 
	)
{
	ofstream o_file("diffuse_on_graph_matA.txt"); 
	std::stringstream out_string;
	out_string << "{"<< endl; // Big {

	// Begin: output graph vts
	out_string << "{ (*Begin: vts*)"<< endl; 
	for (unsigned i = 0; i < m_gVts.size(); ++i)
	{
		const auto & p = m_gVts[i];
		out_string << p[0]<<","<<p[1]<<","<<p[2];
		if (i != m_gVts.size()-1)
		{
			out_string <<",";
		}
	}
	out_string << endl;
	out_string << "}, (*End: vts*)"<< endl; 
	// End: output graph vts

	// Begin: output graph edges
	out_string << "{ (*Begin: graph edges*)"<< endl; 
	for (unsigned i = 0; i < m_gEdges.size(); ++i)
	{
		const auto & e = m_gEdges[i];
		out_string << e[0]+1<<","<<e[1]+1;
		if (i != m_gEdges.size()-1)
		{
			out_string <<",";
		}
	}
	out_string << endl;
	out_string << "}, (*End: graph edges*)"<< endl; 
	// End: output graph faces

	// Begin: output new vi to orig vi map
	out_string << "{ (*Begin: new_to_orig_vi_map*)"<< endl; 
	for ( int i = 0; i < m_new_orig_vertIdx_map.size(); ++i )
	{
		out_string << m_new_orig_vertIdx_map[i]<<",";
	}
	out_string.seekp(-1, out_string.cur);
	out_string << "}, (*End: new_to_orig_vi_map*)"<< endl; 

	// End: output new vi to orig vi map

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

	// Begin: output unknown flags
	out_string << "{ (*Begin: flags (0: unknown 1: known) *)"<< endl; 
	for ( int  i = 0; i < m_isKnown.size(); ++i )
	{
		out_string << (m_isKnown[i] ? 1 : 0) << ",";
	}
	out_string.seekp(-1, out_string.cur);
	out_string << "} (*End: flags*)" <<endl;
	// End: output unknown flags

	out_string << "}" << endl;// Big }

	o_file << out_string.rdbuf();
	o_file.close();
}

void GraphDiffuser::solve_diffusion_system( 
	const vector<float>& _known, 
	vector<float>& _unknown )
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