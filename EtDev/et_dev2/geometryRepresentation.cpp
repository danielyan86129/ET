#include "geometryRepresentation.h"

#include <fstream>
#include <deque>
#include <set>
#include <queue>
#include <stack>
#include <limits>
#include <iostream>
#include <boost/pool/pool_alloc.hpp> 
#include <boost/heap/fibonacci_heap.hpp> // implementaion for Q during burning
#define _USE_MATH_DEFINES
#include <cmath>
#include <random>
#include <ctime>

#include <QString>

#include <XForm.h> // trimesh's transformation matrix

#include <omp.h> // open-mp

// for geodesic-based MA-prune metric
#include "geodesic/geodesic_algorithm_exact.h"

#include "commonUtils.h"
#include "hybridSkeleton.h"
#include "triangulation.h"

#ifdef PRINT_MEM_USAGE
#include "windows.h"
#include "psapi.h"
#endif // PRINT_MEM_USAGE

void preprocess(shared_ptr<MyMesh> _m, vector<float>& _radii, float _eps)
{
	// -1. first find all faces whose spheres cannot produce intersections
	// then run a clustering on vertices of these faces.
	typedef double ldouble;
	typedef Vec<3,ldouble> DPoint;
	typedef Vec<3,ldouble> DVec;
	TriFace f;
	DPoint a, b, c;
	TriColor color;
	set <unsigned> non_isect_faces_vts;
	// run intersector on all faces to find non-isect faces
	/*for (unsigned fi = 0; fi < _m->faces.size(); ++fi)
	{
	f = _m->faces[fi];
	a = DPoint(_m->vertices[f[0]]);
	b = DPoint(_m->vertices[f[1]]);
	c = DPoint(_m->vertices[f[2]]);
	if ( !_radii.empty() )
	{
	ra = _radii[f[0]];
	rb = _radii[f[1]];
	rc = _radii[f[2]];
	}
	util::IsectType isect_t = 
	util::threeSphereIntersect(a, b, c, ra, rb, rc, vector<DPoint>());
	if (isect_t == util::NONE_ISECT)
	{
	non_isect_faces_vts.insert(f[0]);
	non_isect_faces_vts.insert(f[1]);
	non_isect_faces_vts.insert(f[2]);
	}
	}*/

	_m->need_neighbors();
	vector<unsigned> component;
	queue<unsigned> component_front;
	vector<bool> visited(_m->vertices.size(), false);
	unsigned cluster_cnt = 0;
	unsigned vts_clustered_cnt = 0;

	// cluster within non-intersected faces' vertices
	//for (unsigned vi = 0; vi < _m->vertices.size(); ++vi)
	//{
	//	auto vit = non_isect_faces_vts.find(vi);
	//	if ( visited[vi] || vit == non_isect_faces_vts.end() )
	//		continue;

	//	component.clear();
	//	component_front.push(vi);
	//	visited[vi] = true;

	//	while( !component_front.empty() )
	//	{
	//		unsigned cur_vi = component_front.front();
	//		component_front.pop();
	//		component.push_back(cur_vi);

	//		vector<int>* nbs = &(_m->neighbors[cur_vi]);
	//		for ( unsigned i = 0; i < nbs->size(); ++i )
	//		{
	//			unsigned ni = (*nbs)[i];
	//			if ( 
	//				!visited[ni] && 
	//				non_isect_faces_vts.find(ni) != non_isect_faces_vts.end() &&
	//				trimesh::dist(_m->vertices[ni], _m->vertices[cur_vi]) < _eps
	//				)
	//			{
	//				component_front.push(ni);
	//				visited[ni] = true;
	//			}
	//		}
	//	}

	//	if ( component.size() <= 1 )
	//		continue;

	//	cluster_cnt ++;
	//	vts_clustered_cnt += component.size();

	//	// re-assign coord. of all vts in this component to a representative's
	//		for (unsigned i = 0; i < component.size(); ++i)
	//			_m->vertices[component[i]] = _m->vertices[component[0]];
	//}

	cout << "# clusters clustered: " << cluster_cnt << endl;
	cout << "# vts clustered: " << vts_clustered_cnt << endl;

	// basic replacement: a vert id -> the vert to replace to
	vector<unsigned> replace_map;
	replace_map.resize(_m->vertices.size());
	for (unsigned i = 0; i < _m->vertices.size(); ++i)
	{
		replace_map[i] = i;
	}
	vector<bool> face_2_remove(_m->faces.size(), false);
	unsigned face_remove_cnt = 0;

	// 0.	cluster very close vertices
	//	using flooding: first identify almost-degenerate connected component,
	//	then re-assign coord. of the whole component to the representative vert
	//visited.assign(_m->vertices.size(), false);
	//for (unsigned vi = 0; vi < _m->vertices.size(); ++vi)
	//{
	//	if ( visited[vi] )
	//		continue;

	//	component.clear();
	//	component_front.push(vi);
	//	visited[vi] = true;

	//	while( !component_front.empty() )
	//	{
	//		unsigned cur_vi = component_front.front();
	//		component_front.pop();
	//		component.push_back(cur_vi);

	//		vector<int>* nbs = &(_m->neighbors[cur_vi]);
	//		for ( unsigned i = 0; i < nbs->size(); ++i )
	//		{
	//			unsigned ni = (*nbs)[i];
	//			if ( !visited[ni] && trimesh::dist(_m->vertices[ni], _m->vertices[cur_vi]) < _eps )
	//			{
	//				component_front.push(ni);
	//				visited[ni] = true;
	//			}
	//		}
	//	}

	//	if ( component.empty() )
	//		continue;

	//	// re-assign coord. of all vts in this component to a representative's
	//	for (unsigned i = 0; i < component.size(); ++i)
	//		_m->vertices[component[i]] = _m->vertices[component[0]];
	//}

	// 1.	build basic replacement rules
	//	using flooding: first identify degenerate connected component,
	//	then rename the whole component to the representative vert
	visited.assign(_m->vertices.size(), false);
	for (unsigned vi = 0; vi < _m->vertices.size(); ++vi)
	{
		if ( visited[vi] )
			continue;

		component.clear();
		component_front.push(vi);
		visited[vi] = true;

		while( !component_front.empty() )
		{
			unsigned cur_vi = component_front.front();
			component_front.pop();
			component.push_back(cur_vi);

			vector<int>* nbs = &(_m->neighbors[cur_vi]);
			for ( unsigned i = 0; i < nbs->size(); ++i )
			{
				unsigned ni = (*nbs)[i];
				if ( !visited[ni] && _m->vertices[ni] == _m->vertices[cur_vi] )
				{
					component_front.push(ni);
					visited[ni] = true;
				}
			}
		}

		if ( component.empty() )
			continue;

		// rename all vts in this component to a representative
		for (unsigned i = 0; i < component.size(); ++i)
			replace_map[component[i]] = component[0];
	}

	// remove degenerate faces
	for (unsigned fi = 0; fi < _m->faces.size(); ++fi)
	{
		f = _m->faces[fi];
		if ( util::anySame(_m->vertices[f[0]], _m->vertices[f[1]], _m->vertices[f[2]], -1.0f) )
		{
			face_2_remove[fi] = true;
			face_remove_cnt ++;
		}
	}
	if (face_remove_cnt == 0)
		return;
	trimesh::remove_faces(_m.get(), face_2_remove);

	// 2.   remove unused vertices by going through and renaming the vertices
	//		that orig. vertices are mapped to
	map <unsigned, unsigned> mapped2vts;
	for (unsigned vi = 0; vi < replace_map.size(); ++vi)
	{
		auto mapped_it = mapped2vts.find(replace_map[vi]);
		// renaming mapped-to-vert to its expected position after unused vertices removed
		if ( mapped_it == mapped2vts.end() )
		{
			unsigned sz = mapped2vts.size();
			mapped2vts.insert(mapped_it, make_pair(replace_map[vi], sz));
			replace_map[vi] = sz;
		}
		else
		{
			replace_map[vi] = mapped_it->second;
		}
	}

	// 4.	re-position radii according to the expected new position of each vert. 
	//		also, re-position vertices
	if ( !_radii.empty() )
	{
		vector<float> radii(mapped2vts.size());
		for (auto it = mapped2vts.begin(); it != mapped2vts.end(); ++it)
		{
			radii[it->second] = _radii[it->first];
		}

		// swap
		_radii.swap(radii);
	}
	vector<trimesh::point> vertices(mapped2vts.size());
	for (auto it = mapped2vts.begin(); it != mapped2vts.end(); ++it)
	{
		vertices[it->second] = _m->vertices[it->first];
	}
	_m->vertices.swap(vertices);
	vertices.clear(); 
	vertices.shrink_to_fit();

	// 5.	finally, rename each face's vertices according to replacement map
	for ( unsigned fi = 0; fi < _m->faces.size(); ++fi )
	{
		TriFace* f = &(_m->faces[fi]);
		for ( unsigned i = 0; i < 3; ++i )
		{
			unsigned vi = (*f)[i];
			// rename
			(*f)[i] = replace_map[vi];
		}
	}

	//trimesh::remove_unused_vertices(_m.get());

	// verify
	face_remove_cnt = 0;
	for (unsigned fi = 0; fi < _m->faces.size(); ++fi)
	{
		f = _m->faces[fi];
		if ( util::anySame(_m->vertices[f[0]], _m->vertices[f[1]], _m->vertices[f[2]], -1.0f) )
		{
			face_remove_cnt ++;
		}
	}
	if (face_remove_cnt == 0)
		cout << "Verifying... # faces still need to remove (if not 0 then something is wrong...): " << face_remove_cnt << endl;

}

void find2Closest(shared_ptr<MyMesh> _mMA, shared_ptr<MyMesh> _m3d, vector<TriPoint>& _2closest)
{
	// kd tree built upon m3d's vertices
	float* points = new float[3*_m3d->vertices.size()];
	TriPoint p;
	for (unsigned i = 0; i < _m3d->vertices.size(); ++i)
	{
		p = _m3d->vertices[i];
		points[3*i+0] = p[0];
		points[3*i+1] = p[1];
		points[3*i+2] = p[2];
	}
	trimesh::KDtree kdtree(points, _m3d->vertices.size());
	_2closest.resize(_mMA->faces.size()*2); // 2 points for each face of MA mesh

	// find the 2 closest points in _m3d for the center of each face in _mMA
	TriFace f;
	TriPoint center;
	vector<const float *> knn;
	float cos_30 = std::sqrt(3.0f) / 2.0f;
	unsigned close_cnt = 0;
	for (unsigned fi = 0; fi < _mMA->faces.size(); ++fi)
	{
		f = _mMA->faces[fi];
		center = 1 / 3.0f * ( _mMA->vertices[f[0]] + _mMA->vertices[f[1]] + _mMA->vertices[f[2]] );
		knn.clear();
		kdtree.find_k_closest_to_pt(knn, 2, center, -1.0f);
		// copy 2 nearest points to the return list
		_2closest[2*fi + 0] = TriPoint(knn[0]);
		_2closest[2*fi + 1] = TriPoint(knn[1]);

		/// debug
		// check if 2 closest points are on the same side of the face
		if ( !util::anySame(_mMA->vertices[f[0]], _mMA->vertices[f[1]], _mMA->vertices[f[2]], -1.0f) )
		{
			vec normal = trimesh::normalize( vec(
				(_mMA->vertices[f[0]] - _mMA->vertices[f[1]]) CROSS (_mMA->vertices[f[0]] - _mMA->vertices[f[2]])
				) );
			vec v1 = _2closest[2*fi + 0] - center;
			trimesh::normalize(v1);
			vec v2 = _2closest[2*fi + 1] - center;
			trimesh::normalize(v2);

			if ( (v1 DOT normal) * (v2 DOT normal) > 0 )
			{
				close_cnt ++;
				/*cout << "tri <abc> -----------"
				<<fi<<":" << f[0]<<", "<<f[1]<<", "<<f[2]<< endl;
				cout << "a, b, c: "
				<<_mMA->vertices[f[0]]<<","<<_mMA->vertices[f[1]]<<","<<_mMA->vertices[f[2]]<<endl;
				cout << "closest points on same side:"
				<<_2closest[2*fi + 0]<<", "<<_2closest[2*fi + 1]<<endl;*/
			}
		}
	}
	cout << "# of faces whose 2 closest points on same side: "<<close_cnt<<endl;
}

void SteinerGraph::printBurnTime(vector<size_t>& _vts_indices)
{

	/* print burn time for a few sample points each with different topology */
	set<TopoType> topo_found;
	std::default_random_engine e(1);
	std::uniform_int_distribution<> distr(0, this->m_origG->vts.size() - 1);
	for (size_t i = 0; i < 10; ++i)
	{
		//size_t vi = float(e()) / e.max() * (this->bt2MA_vert.size() - 1);
		size_t vi = distr(e);
		if (this->getTopoType(vi) != NONMANIFOLD_2D)
			i --;
		else
			_vts_indices.push_back(vi);
	}
	for (size_t i = 0; i < 10; ++i)
	{
		//size_t vi = float(e()) / e.max() * (this->bt2MA_vert.size() - 1);
		size_t vi = distr(e);
		if (this->getTopoType(vi) != MANIFOLD_2D)
			i --;
		else
			_vts_indices.push_back(vi);
	}

	std::streamsize ss = std::cout.precision();
	auto flgs = cout.flags();
	cout.setf(ios::left, ios::adjustfield);
	for (size_t i = 0; i < _vts_indices.size(); ++i)
	{
		auto vi = _vts_indices[i];
		cout.precision(numeric_limits<float>::max_digits10);
		cout << "sample vert ";
		cout << setw(8) << vi;
		cout <<" ("<< topo_type_str[this->getTopoType(vi)]<<") burnt time\t";
		cout <<this->bt2MA_vert[vi]<<endl;
	}
	cout.precision(ss);
	cout.setf(flgs);
}

//////////////////////////////////////////////////////////////////////////
/// -----------definition of SteinerGraph-----------
//////////////////////////////////////////////////////////////////////////

SteinerGraph::SteinerGraph(
	const shared_ptr<MyGraph>& _g, 
	const vector<float>& _dist2surf, 
	SteinerSubdivision::SubdivideScheme _scheme, double _steiner_param,
	int _edgeWeight_idx,
	SteinerGraph::BurnScheme _burn_scheme/* = BurnScheme::ORIGINAL_AND_STEINER*/
	)
{
#ifdef PRINT_MEM_USAGE
	PROCESS_MEMORY_COUNTERS pmc;
	GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc));
	SIZE_T memPreUsed, memCurUsed;
	memPreUsed = memCurUsed = pmc.WorkingSetSize;
	cout << "mem. usage (starting SteinerGraph() ) in MB: " << memCurUsed / (1024*1024) << endl;
#endif // PRINT_MEM_USAGE

	this->m_origG = _g;
	/*this->st_vts = orig_g.vts;*/

	init();

	// set general scheme
	m_curBurnScheme = _burn_scheme;

	// remember current edge weight of choice
	this->setEdgeWeightFunc(_edgeWeight_idx);

	// steiner-subdivision 
	cout << "-------performing steiner subdivision--------" << endl;
	this->steinerSubdivide(_steiner_param, _scheme);
	this->computeBT3ForSteinerVts(_dist2surf);
	cout << "-------steiner subdivision done!--------" << endl;

	cout << "# orig triangle vts: " << m_origG->vts.size() << endl;
	cout << "# orig tri faces: " << m_origG->faces.size() << endl;
	cout << "# orig tri edges: " << m_origG->edges.size() << endl;
	cout << "# total st graph vts: " << m_stSubdiv.sizeOfVts() << endl;
	cout << "# total st graph edges: " << m_stSubdiv.sizeOfStEdges() << endl;

#ifdef PRINT_MEM_USAGE
	GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc));
	memCurUsed = pmc.WorkingSetSize;
	cout << "mem. usage: " << memCurUsed / (1024*1024) << endl;
	cout << "mem. usage (by steinerSubdivide() ) in MB: " << (memCurUsed - memPreUsed) / (1024*1024) << endl;
	memPreUsed = memCurUsed;
#endif // PRINT_MEM_USAGE

	cout << "-------building topo graphs--------" << endl;
	this->buildTopoGraphs();
	cout << "-------building topo graphs done!--------" << endl;

#ifdef PRINT_MEM_USAGE
	GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc));
	memCurUsed = pmc.WorkingSetSize;
	cout << "mem. usage: " << memCurUsed / (1024*1024) << endl;
	cout << "mem. usage (by buildTopoGraphs() ) in MB: " << (memCurUsed - memPreUsed) / (1024*1024) << endl;
	memPreUsed = memCurUsed;
#endif // PRINT_MEM_USAGE

	this->printTopoStats();

	cout << "-------burning--------" << endl;
	this->burn();
	cout << "-------burning done!--------" << endl;
	cout << endl;

	cout << "-------generating finner triangulation for MA faces--------" << endl;
	this->triangulate();
	cout << "-------triangulation done!--------" << endl;
	cout << endl;

#ifdef PRINT_MEM_USAGE
	GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc));
	memCurUsed = pmc.WorkingSetSize;
	cout << "mem. usage: " << memCurUsed / (1024*1024) << endl;
	cout << "mem. usage (by burn() ) in MB: " << (memCurUsed - memPreUsed) / (1024*1024) << endl;
	memPreUsed = memCurUsed;
#endif // PRINT_MEM_USAGE

	cout << "----------unburnt vts & their topo graphs----------" << endl;
	unsigned cnt = 10;
	for (unsigned i = 0, c = 0; i < m_stSubdiv.sizeOfVts() && c < cnt; ++i)
	{
		if (!burnt[i] && t_graphs[i] != nullptr)
		{
			cout << "Topo graph for " << i << ": " << endl;
			cout << t_graphs[i];
			++c;
		}
	}

}

SteinerGraph::~SteinerGraph()
{
	for (size_t i = 0; i < t_graphs.size(); ++i)
		del_topo(i);
}

SteinerGraph::BurnScheme SteinerGraph::getBurnScheme() const
{
	return m_curBurnScheme;
}

static TriColor topo_colors[5] = {
	TriColor(1.0f, 1.0f, 0.0f), // ISOLATE
	TriColor(1.0f, 0.0f, 0.0f),	// BOUNDARY-2d
	TriColor(0.0f, 0.0f, 0.0f), // JUNCTION-2d
	TriColor(0.7f, 0.7f, 0.7f),	// MANIFOLD-2d
	TriColor(0.0f, 0.0f, 1.0f)	// NONMANIFOLD-2d
};

const char* topo_type_str[7] = {
	"ISOLATE", "MANIFOLD_1D", "JUNCTION_1D", 
	"BOUNDARY_2D", "JUNCTION_2D", "MANIFOLD_2D", "NONMANIFOLD_2D"
};

SteinerGraph* TopoGraph::stg = NULL;

TopoGraph::TopoGraph()
{
	type = ISOLATE;
};

bool TopoGraph::is_topo_edge_open(TopoGraph::EdgeIdx _i) const
{
	TriEdge te = t_edges[_i].first;

	if (!util::notPureLoop(te))
		return false;

	for (unsigned j = 0; j < 2; ++j)
	{

		const auto nbs = t_nbVtsMap.find(te[(std::size_t)j]);
		assert(nbs != t_nbVtsMap.end());
		if (nbs->second.size() == 1)
			return true;

	}

	return false;
}

bool TopoGraph::has_any_open_topo_edge(void) const
{
	for (EdgeIdx i = 0; i < t_edges.size(); ++i)
		if (is_topo_edge_open(i))
			return true;
	return false;
}

int TopoGraph::getEmptyInstanceIdxOfTEdge(TriEdge _te) const
{
	const auto& find_iter = t_edgeIdxesMap.find(_te);
	if (find_iter == t_edgeIdxesMap.end())
		return -1;

	const auto& te_indices = find_iter->second;
	for (auto idx_it = te_indices.begin(); idx_it != te_indices.end(); ++idx_it)
	{
		if (t_edges[*idx_it].second.empty())
		{
			return *idx_it;
		}
	}

	// no such instance found
	return -1;
}

int TopoGraph::nConnCmpnt() const
{
	unsigned cc_cnt = 0;
	set<unsigned> visited;
	queue<unsigned> q;
	for (unsigned vi = 0; vi < t_vts.size(); ++vi)
	{
		unsigned cur_tv = t_vts[vi];
		if (visited.count(cur_tv))
			continue;

		cc_cnt++;
		visited.insert(cur_tv);
		q.push(cur_tv);

		while(!q.empty())
		{
			cur_tv = q.front();
			q.pop();

			const auto& nbs = t_nbVtsMap.find(cur_tv)->second;
			for (auto it = nbs.begin(); it != nbs.end(); ++it)
			{
				if (visited.count(*it) == 0)
				{
					visited.insert(*it);
					q.push(*it);
				}
			}
		}
	}

	return cc_cnt;
}

ostream& operator << (ostream& _os, const TopoGraph& _tg)
{
	_os << "topo type: " << topo_type_str[_tg.type] << endl;
	_os << "# topo vts: " << _tg.t_vts.size() << endl;
	_os << "# topo edges: " << _tg.t_edges.size() << endl;
	if (_tg.t_edges.size()) 
		_os << "they are: ";
	for (auto te_itor = _tg.t_edges.begin(); te_itor != _tg.t_edges.end(); ++te_itor)
	{
		_os << te_itor->first << ", assoc with faces ";
		for (auto f_it = te_itor->second.begin(); f_it != te_itor->second.end(); ++f_it)
		{
			_os << *f_it<<" ";
		}
		_os << endl;
	}

	auto stg = TopoGraph::stg;
	_os << "mapping: orig face -> topo edge: " << _tg.t_origFace2TopoMap.size() << endl;
	for (auto it = _tg.t_origFace2TopoMap.begin(); it != _tg.t_origFace2TopoMap.end(); ++it)
	{
		auto fi = it->first;
		const auto& f = stg->m_origG->faces[fi];
		_os <<"("<<fi<<" ("
			<<f[0]<<" "<<topo_type_str[stg->getTopoType(f[0])]<<","
			<<f[1]<<" "<<topo_type_str[stg->getTopoType(f[1])]<<","
			<<f[2]<<" "<<topo_type_str[stg->getTopoType(f[2])]<<")"<<")"
			<<" -> ("<< it->second<<")"<<endl;
	}

	_os << endl<< endl;
	return _os;
}

unsigned TopoGraph::size()
{
	unsigned sz = 0;
	// size of topo vts
	sz += t_vts.size() * sizeof(int);
	// size of topo edges
	for (auto te = t_edges.begin(); te != t_edges.end(); ++te)
	{
		sz += sizeof(te->first);
		sz += te->second.capacity() * sizeof(TriEdge) + sizeof(te->second);
	}
	// size of nbVtsMap
	sz += sizeof(t_nbVtsMap);
	for (auto itor = t_nbVtsMap.begin(); itor != t_nbVtsMap.end(); ++itor)
	{
		sz += sizeof(itor->first);
		sz += itor->second.capacity() * sizeof(int) + sizeof(itor->second);
	}
	// size of edgeIdxesMap
	sz += sizeof(t_edgeIdxesMap);
	for (auto itor = t_edgeIdxesMap.begin(); itor != t_edgeIdxesMap.end(); ++itor)
	{
		sz += sizeof(itor->first);
		sz += itor->second.capacity() * sizeof(EdgeIdx) + sizeof(itor->second);
	}
	// size of map<TriEdge, pair<char, int>> t_origEdge2TopoMap
	sz += sizeof(t_origFace2TopoMap);
	for (auto itor = t_origFace2TopoMap.begin(); itor != t_origFace2TopoMap.end(); ++itor)
	{
		sz += sizeof(itor->first);
		sz += sizeof(itor->second);
	}

	return sz;
}

void TopoGraph::reserve()
{
	t_vts.reserve(20);
	t_edges.reserve(20);
}

void TopoGraph::shrinkToFit()
{
	t_vts.shrink_to_fit();
	t_edges.shrink_to_fit();
}

bool TopoGraph::isOpenSheet( int _vi, int _si )
{
	auto tt = stg->getTopoType( _vi );
	return tt == BOUNDARY_2D || tt == JUNCTION_2D && (stg->t_graphs[ _vi ])->is_topo_edge_open( _si );
}

bool TopoGraph::hasAnyOpenSheet( int _vi )
{
	return stg->getTopoType( _vi ) == BOUNDARY_2D || 
		stg->getTopoType( _vi ) == JUNCTION_2D &&
		(stg->t_graphs[ _vi ])->has_any_open_topo_edge();
}

bool TopoGraph::has2dNbhood( int _vi )
{
	auto tt = stg->getTopoType( _vi );
	return tt != ISOLATE && tt != JUNCTION_1D && tt != MANIFOLD_1D;
}

TriColor SteinerGraph::getTopoColor(int _i)
{
	assert(!t_graphs.empty());
	return topo_colors[getTopoType(_i)];
}

TopoType SteinerGraph::getTopoType(int _i) const
{
	assert(!t_type.empty());
	return t_type[_i];
}

bool SteinerGraph::is_leading_topo_edge(unsigned _v, TopoGraph::EdgeIdx _tei) const
{
	assert(!prev_vert.empty());
	return prev_vert[_v][_tei] != -2;
}

int SteinerGraph::getPrevVert(unsigned _v, TopoGraph::EdgeIdx _tei) const
{
	return is_leading_topo_edge(_v, _tei) ? prev_vert[_v][_tei] : prev_vert[_v][prev_tedge[_v][_tei]];
}

int SteinerGraph::getPrevTopoEdge(unsigned _v, TopoGraph::EdgeIdx _tei) const
{
	return is_leading_topo_edge(_v, _tei) ? prev_tedge[_v][_tei] : prev_tedge[_v][prev_tedge[_v][_tei]];
}

/*
**private helpers
*/

void SteinerGraph::steinerSubdivide(double _steiner_param, SteinerSubdivision::SubdivideScheme _scheme)
{
	cout << "steiner scheme: "<<_scheme<<endl;
	switch(_scheme)
	{
	case SteinerSubdivision::FIXED:
		//fixedSubdivide(_nSamples);
		m_stSubdiv.fixedSubdivide(this->m_origG, _steiner_param);
		break;
	case SteinerSubdivision::ADAPTIVE:
		m_stSubdiv.adaptiveSubdivide(this->m_origG, _steiner_param);
		break;
	case SteinerSubdivision::ADAPTIVE_MIDPOINT:
		m_stSubdiv.adaptiveMidPointSubdivide(this->m_origG, _steiner_param);
		break;
	default:
		//fixedSubdivide(_nSamples);
		m_stSubdiv.fixedSubdivide(this->m_origG, _steiner_param);
		break;
	}

}

/*
// subdivide the graph with fixed scheme
void SteinerGraph::fixedSubdivide(int _nSamples)
{
if (_nSamples <= 0)
_nSamples = 1;

// increase space to accommodate new steiner vts
st_vts.reserve(_nSamples*m_origG.edges.size());

// generate new steiner vts on each edge
// and connect each steiner vert to its neighbor on this edge
float interp_weight;
for (unsigned ei = 0; ei < this->m_origG.edges.size(); ++ei)
{
TriEdge orig_e = m_origG.edges[ei];
vector<int>& stVtsOnE = triE_stVts_map[orig_e];
for (int i = 0; i < _nSamples; ++i)
{
interp_weight = float(i+1) / (_nSamples+1);
st_vts.push_back(
trimesh::mix(this->st_vts[orig_e[0]], st_vts[orig_e[1]], interp_weight)
);
stVtsOnE.push_back(st_vts.size() - 1);
}

// v0 to steiner vert 1
st_edges.push_back(util::makeEdge(orig_e[0], stVtsOnE[0]));

// last steiner vert to v1
st_edges.push_back(util::makeEdge(stVtsOnE[stVtsOnE.size() - 1], orig_e[1]));

// steiners in between v0 & v1
if ( stVtsOnE.size() > 1 )
for (unsigned st_i = 0; st_i < stVtsOnE.size()-1; ++st_i)
{
st_edges.push_back(util::makeEdge(stVtsOnE[st_i], stVtsOnE[st_i+1]));
}
}

// now connect each steiner vert to verts on other edges, face by faces:
for (unsigned fi = 0; fi < m_origG.faces.size(); ++fi)
{
TriFace f = m_origG.faces[fi];

// 1. connect to opposite orig vert
for (int i = 0; i < 3; ++i)
{
int oppo_v = f[i]; // orig vert opposit to curr edge
TriEdge cur_e = util::makeEdge( f[(i+1) % 3], f[(i+2) % 3] ); // curr edge
vector<int>& st_vts_onE = triE_stVts_map[cur_e]; // steiner vts on curr edge

for (unsigned st_i = 0; st_i < st_vts_onE.size(); ++st_i)
{
st_edges.push_back( util::makeEdge(oppo_v, st_vts_onE[st_i]) );
}
}

// 2. connect to st_vts on other edge
vector<TriEdge> e_of_f = util::edgesFromFace(f);
for (int i = 0; i < 3; ++i)
{
vector<int>& st_vts_currE = triE_stVts_map[e_of_f[i]];
vector<int>& st_vts_nextE = triE_stVts_map[e_of_f[(i+1) % 3]];
for (unsigned j = 0; j < st_vts_currE.size(); ++j)
for (unsigned k = 0; k < st_vts_nextE.size(); ++k)
{
st_edges.push_back( util::makeEdge(st_vts_currE[j], st_vts_nextE[k]) );
}
}
}
}
*/

void SteinerGraph::computeBT3ForSteinerVts(const vector<float>& _bt3_orig)
{
	radii_vert = _bt3_orig;
	if (radii_vert.empty())
	{
		radii_vert.assign(m_stSubdiv.sizeOfVts(), 0.0f);
		bt3MA_vert = radii_vert;
		return;
	}

	// for each edge, determine its steiner vts dist to surf
	// using c^2 = a^2 + b^2 - 2*a*b*Cos(angle between b, c)
	// if no valid distance available, then linear interpolate the distance.
	bool valid_dist;
	vector<unsigned> stVtsOnE;
	for (unsigned ei = 0; ei < m_origG->edges.size(); ++ei)
	{
		TriEdge e = m_origG->edges[ei];
		unsigned u = e[0], v = e[1];
		float ru = radii_vert[u], rv = radii_vert[v];
		TriPoint p_u = m_stSubdiv.getVert(u);
		TriPoint p_v = m_stSubdiv.getVert(v);
		float d_uv2 = trimesh::dist2(p_u, p_v);
		float d_uv = sqrt(d_uv2);
		float cosine, d_st2, r_st;
		stVtsOnE.clear();
		m_stSubdiv.getStVertIndicesOnTriEdge(ei, stVtsOnE);
		int st_v;
		bool valid_dist;

		if (ru + rv > d_uv && std::abs(ru - rv) < d_uv)
		{
			valid_dist = true;
			cosine = 
				( d_uv2 + pow(ru, 2.0f) - pow(rv, 2.0f) ) / 
				(2.0f * ru * sqrt(d_uv2));
		}
		else
		{
			valid_dist = false;
		}

		d_uv = d_uv < 1e-8 ? 1.0f : d_uv;
		for (unsigned i = 0; i < stVtsOnE.size(); ++i)
		{
			// use linear interpolation
			float a = float(stVtsOnE.size() - i) / (1.0f+stVtsOnE.size());
			//float a = trimesh::dist(m_stSubdiv.getStVert(stVtsOnE[i]), p_v) / d_uv;
			r_st = a * ru + (1.0f - a) * rv;

			// append the dist-to-surf of steiner vert to the list
			radii_vert.push_back(r_st);
		}
	}

	if (this->m_curEdgeWt == SteinerGraph::EW_EUCLIDEAN) // uniform burning needs bt3
	{
		bt3MA_vert = radii_vert;
	}
	else if (this->m_curEdgeWt == SteinerGraph::EW_BT2BT3DIFF) // bt2bt3diff burning needs bt3 to be all zero
	{
		bt3MA_vert.assign(radii_vert.size(), 0.0f);
	}
	else // by default assume uniform burning
	{
		bt3MA_vert = radii_vert;
	}
}

void SteinerGraph::buildTopoGraphs()
{
#ifdef PRINT_MEM_USAGE
	PROCESS_MEMORY_COUNTERS pmc;
	GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc));
	SIZE_T memCurUsed = pmc.WorkingSetSize;
#endif // PRINT_MEM_USAGE

	TopoGraph::stg = this;

	//cout << "entering buildTopoGraphs()..." << endl;
	t_graphs.assign(m_stSubdiv.sizeOfVts(), nullptr);
	t_holds_actual_storage.assign(m_stSubdiv.sizeOfVts(), false);
	t_type.assign(m_stSubdiv.sizeOfVts(), ISOLATE);
	t_n_sector.assign(m_stSubdiv.sizeOfVts(), 0);
	resetTopoStats();

#ifdef PRINT_MEM_USAGE
	GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc));
	cout << "mem. usage ( due to t_graphs.resize() ) in MB: " << (pmc.WorkingSetSize - memCurUsed) / (1024*1024) << endl;
	memCurUsed = pmc.WorkingSetSize;
#endif // PRINT_MEM_USAGE

	typedef vector<int> OneRingVList;
	typedef map<int, vector<int>> NbVMap;
	typedef vector<TriEdge> AssocEdges;
	typedef map<TriEdge, vector<int>> EdgeIDMap;

	vector<int> one_ring_vts; one_ring_vts.reserve(30);
	map<int, vector<int>> nbVts_map;
	TopoGraph::TEdgeList one_ring_edges; one_ring_edges.reserve(30);
	vector<int> faces_assoc; faces_assoc.reserve(30);
	map<TriEdge, vector<int>> edgeIdxes_map; 

	vector<unsigned> st_vts_onE;

	//cout << "starting build process..." << endl;
	// 1. topo graph for orig. vts
	for (unsigned vi = 0; vi < m_origG->vts.size(); ++vi)
	{
		if (vi % 100000 == 0) //debug
			cout << "# orig topo graphs built: " << vi << endl;

		// find 1-ring neighborhood of vi
		one_ring_vts.clear();
		nbVts_map.clear();
		one_ring_edges.clear();
		faces_assoc.clear();
		edgeIdxes_map.clear();
		vector<int>& nbFaces = m_origG->nbFacesOfVert[vi];
		for (unsigned i = 0; i < nbFaces.size(); ++i)
		{
			faces_assoc.clear(); 

			const auto& f = m_origG->faces[nbFaces[i]];
			TriEdge oppo_e = util::oppositeEdge(f, vi);

			// the steiner vert closest to vi on each incident edge
			int st1 = m_stSubdiv.getClosestVertOnEdge(vi, oppo_e[0], 0);
			int st2 = m_stSubdiv.getClosestVertOnEdge(vi, oppo_e[1], 0);

			// that's our 1-ring edge
			// associate cur face with that 1-ring edge
			// then save this <edge, assoc_face> pair
			TriEdge one_r_e = util::makeEdge(st1, st2);
			/*int ret_code = m_origG->getFaceIdx( util::makeFace( m_origG->faces[ nbFaces[ i ] ] ) );
			assert(ret_code >= 0);*/
			faces_assoc.push_back( nbFaces[ i ] );
			one_ring_edges.push_back(make_pair(one_r_e, faces_assoc));

			// index 1-ring edge and store 1-ring vts (test for duplicate)
			edgeIdxes_map[one_r_e].push_back(one_ring_edges.size()-1);
			unsigned pre_size = nbVts_map.size();
			nbVts_map[st1].push_back(st2);
			/*if (pre_size < nbVts_map.size()) 
				one_ring_vts.push_back(st1);
			pre_size = nbVts_map.size();*/
			nbVts_map[st2].push_back(st1);
			/*if (pre_size < nbVts_map.size()) 
				one_ring_vts.push_back(st2);*/
		}

		// merge degree-2 vts in 1-ring
		for ( auto it = nbVts_map.begin(); it != nbVts_map.end(); ++it )
		{
			int nb_v0 = it->first;
			vector<int>& nbVts_v = it->second;
			// skip vert that has non-2-degree
			// or loops to itself
			if (
				nbVts_v.size() != 2 || 
				edgeIdxes_map.find(util::makeEdge(nb_v0, nb_v0)) != edgeIdxes_map.end()
				)
				continue;

			faces_assoc.clear();
			int nb_v1 = nbVts_v[0];
			int nb_v2 = nbVts_v[1];
			vector<int>& nbVts_v1 = nbVts_map[nb_v1];
			vector<int>& nbVts_v2 = nbVts_map[nb_v2];

			// case 1 & 2:
			// will create a loop edge after merging
			if (nb_v1 == nb_v2)
			{
				TriEdge new_edge; // will be filled in case branching
				bool loop_created = false;
				// case (1) merge <v0, v1> to obtain loop edge <v1, v1>
				if (nbVts_v1.size() > 2)
				{
					nbVts_v.clear();
					nbVts_v1.erase(
						std::remove(nbVts_v1.begin(), nbVts_v1.end(), nb_v0), 
						nbVts_v1.end());

					// new vert neighbor relationship
					nbVts_v1.push_back(nb_v1);
					new_edge = TriEdge(nb_v1, nb_v1);
					loop_created = true;
				}
				// case (2) merge <v0, v1> to obtain loop edge <-1, -1>, i.e. loop edge w/o node
				else if (nbVts_v1.size() == 2)
				{
					nbVts_v.clear();
					nbVts_v1.clear();

					new_edge = util::pureLoopEdge();
					loop_created = true;
				}

				if (loop_created)
				{
					// also need to merge assoc.ed faces.
					vector<int>& edge_idxes = edgeIdxes_map[util::makeEdge(nb_v0, nb_v1)];
					mergeAssocFaces(faces_assoc, edge_idxes, one_ring_edges);
					edge_idxes.clear();

					// create new loop edge
					one_ring_edges.push_back(make_pair(new_edge, faces_assoc));
					edgeIdxes_map[new_edge].push_back(one_ring_edges.size()-1);
				}
			}
			else// case (3) no loop edge will be created. just merge
			{
				// delete v0 entry from map and v1 v2's neighborhood
				nbVts_map[nb_v0].clear();
				nbVts_v1.erase(
					std::remove(nbVts_v1.begin(), nbVts_v1.end(), nb_v0), 
					nbVts_v1.end()
					);
				nbVts_v2.erase(
					std::remove(nbVts_v2.begin(), nbVts_v2.end(), nb_v0), 
					nbVts_v2.end()
					);

				// assoc faces with new edge <v1, v2>, they are:
				// faces assoc.ed with topo edges <v0, v1>, <v0, v2>
				vector<int>* edge_idxes = &edgeIdxes_map[util::makeEdge(nb_v0, nb_v1)];
				mergeAssocFaces(faces_assoc, *edge_idxes, one_ring_edges);
				edge_idxes->clear();
				edge_idxes = &edgeIdxes_map[util::makeEdge(nb_v0, nb_v2)];
				mergeAssocFaces(faces_assoc, *edge_idxes, one_ring_edges);
				edge_idxes->clear();

				// after merge new edges & nb vts relationship
				nbVts_v1.push_back(nb_v2);
				nbVts_v2.push_back(nb_v1);
				TriEdge new_edge = util::makeEdge(nb_v1, nb_v2);
				one_ring_edges.push_back(make_pair(new_edge, faces_assoc));
				edgeIdxes_map[new_edge].push_back(one_ring_edges.size()-1);
			}
		}

		TopoGraph* tg = make_topo( vi );
		
		// then, only include the 1-cell 1-ring neighbors when there is no 2d bits
		if ( one_ring_edges.empty() )
		{
			const auto& nbVts = m_origG->nbVtsOfVert[ vi ];
			for ( auto i = 0; i < nbVts.size(); ++i )
			{
				//if ( m_origG->getNbFaces( util::makeEdge( vi, nbVts[ i ] ) ).empty() )
				tg->t_vts.push_back( nbVts[ i ] ); // add the nb v to rep this 1-cell-1-ring neighbor
			}
		}

		// finalize topo vts
		vector<NbVMap::iterator> itors2remove;
		for (NbVMap::iterator itor = nbVts_map.begin(); itor != nbVts_map.end(); ++itor)
		{
			if (itor->second.empty())
			{
				itors2remove.push_back(itor);
				continue;
			}

			tg->t_vts.push_back(itor->first);
		}
		for (auto itor = itors2remove.begin(); itor != itors2remove.end(); ++itor)
			nbVts_map.erase(*itor);

		// save nb vert relationship
		tg->t_nbVtsMap = nbVts_map;

		// finalize topo edges
		for (TopoGraph::TEdgeList::iterator te_itor = one_ring_edges.begin(); 
			te_itor != one_ring_edges.end(); ++te_itor)
		{
			auto itor = edgeIdxes_map.find(te_itor->first);
			if (itor->second.empty())
			{
				continue;
			}

			assert(!te_itor->second.empty());
			tg->t_edges.push_back(*te_itor);
			tg->t_edgeIdxesMap[te_itor->first].push_back(tg->t_edges.size()-1);

			// build mapping: orig face -> its associated topo edge
			auto& map = tg->t_origFace2TopoMap;
			for (auto fi_it = te_itor->second.begin(); 
				fi_it != te_itor->second.end(); ++fi_it)
			{
				assert(map.find(*fi_it) == map.end());
				map[*fi_it] = tg->t_edges.size()-1;
			}
		}

		// determine topo type
		TopoType tt = identify(*tg);
		t_type[vi] = tt;
		updateTopoStats(tt);

		// let's save memory: if the topo graph has only one or zero sector, 
		// then we can do w/o actually storing it. 
		t_n_sector[vi] = tg->t_edges.size();
		if (t_n_sector[vi] <= 1)
		{
			del_topo(vi);
		}
	} // end of building topo graph for each original vertex

	// 2. topo graph for steiner vts
	// 
	vector<unsigned> st_onE; // only st vts on one edge
	vector<unsigned> vts_onE;// including the 2 edge end
	vector<unsigned> st_onE1;
	vector<unsigned> st_onE2;
	vector<unsigned> vts_onOtherE;
	// we only need to store one topo graph for each orig edge
	map<unsigned, TopoGraph*> edge_topograph_map;
	// tentatively distribute 1-ring edges & vts into each steiner vert's topo graph
	for (unsigned fi = 0; fi < m_origG->faces.size(); ++fi)
	{
		TriFace f = m_origG->faces[fi];
		vector<TriEdge> edges_of_f = util::edgesFromFace(f);

		for (unsigned ei = 0; ei < 3; ++ei)
		{
			TriEdge e = edges_of_f[ei];

			/// Begin: optimization
			// fetch from map if topo graph has been allocated for it already
			TopoGraph* tg;
			auto e_idx = m_origG->getEdgeIdx(e);
			auto find_it = edge_topograph_map.find(e_idx);
			if (find_it != edge_topograph_map.end())
			{
				tg = find_it->second;
			}
			/// End: optimization

			st_onE.clear();
			m_stSubdiv.getStVertIndicesOnTriEdge(e, st_onE);
			vts_onE.resize(2+st_onE.size());
			vts_onE.front() = e[0];
			vts_onE.back() = e[1];
			std::copy(st_onE.begin(), st_onE.end(), vts_onE.begin()+1);

			/// Begin: optimization
			// only create one topo graph for each edge: 
			// let's pick the first st vert on the edge as the rep. vert
			// and create a topo graph for it
			// all other st vts on this edge should refer to its topo graph
			int v = vts_onE[1];
			int v1 = e[0];
			int v2 = e[1];
			tg = make_topo(v);
			edge_topograph_map[e_idx] = tg;
			// new edge
			faces_assoc.clear();
			TriEdge one_r_e = util::makeEdge(v1, v2);
			faces_assoc.push_back(fi);
			tg->t_edges.push_back(make_pair(one_r_e, faces_assoc));
			tg->t_edgeIdxesMap[one_r_e].push_back(tg->t_edges.size());
			// new topo vts (test duplicate)
			int pre_size = tg->t_nbVtsMap.size();
			tg->t_nbVtsMap[v1].push_back(v2);
			if (tg->t_nbVtsMap.size() != pre_size)
				tg->t_vts.push_back(v1);
			pre_size = tg->t_nbVtsMap.size();
			tg->t_nbVtsMap[v2].push_back(v1);
			if (tg->t_nbVtsMap.size() != pre_size)
				tg->t_vts.push_back(v2);
			// now make every other st vts refer to this topo graph
			for (unsigned i = 2; i < vts_onE.size()-1; ++i)
			{
				make_topo(vts_onE[i], tg);
			}
			/// End: optimization

			/// the unoptimized version
			// new topo edge and faces associated to it
			/*for (unsigned i = 1; i < vts_onE.size()-1; ++i)
			{
				int v = vts_onE[i];
				int v1 = vts_onE[i-1];
				int v2 = vts_onE[i+1];
				TopoGraph* tg = make_topo(v);

				// new edge
				faces_assoc.clear();
				TriEdge one_r_e = util::makeEdge(v1, v2);
				faces_assoc.push_back(fi);
				tg->t_edges.push_back(make_pair(one_r_e, faces_assoc));
				tg->t_edgeIdxesMap[one_r_e].push_back(tg->t_edges.size());
				// new topo vts (test duplicate)
				int pre_size = tg->t_nbVtsMap.size();
				tg->t_nbVtsMap[v1].push_back(v2);
				if (tg->t_nbVtsMap.size() != pre_size)
					tg->t_vts.push_back(v1);
				pre_size = tg->t_nbVtsMap.size();
				tg->t_nbVtsMap[v2].push_back(v1);
				if (tg->t_nbVtsMap.size() != pre_size)
					tg->t_vts.push_back(v2);
			}*/
		}
	}
	// for each steiner [merge those 1-ring steiners of degree 2, and] finalize topo graph
	vector<int> edges_idxes;
	edges_idxes.push_back(0);
	edges_idxes.push_back(1);
	for (unsigned vi = m_stSubdiv.getStartOfStVts(); vi < m_stSubdiv.getEndOfStVts(); ++vi)
	{
		if ((vi - m_stSubdiv.getStartOfStVts()) % 100000 == 0) //debug
			cout << "# steiner topo graphs built: " << (vi - m_stSubdiv.getStartOfStVts()) << endl;

		/// Begin: optimization
		// if there is no actual topo graph allocated for vi, 
		// then vi will refer to the rep vert's topo graph
		if (!t_holds_actual_storage[vi])
		{
			continue;
		}
		// otherwise, vi is a rep vert for its edge
		// we need to go ahead and finalize its topo graph
		/// End: optimization
		
		//TopoGraph* tg = make_topo(vi);
		TopoGraph* tg = t_graphs[vi];

		if (tg->t_edges.size() > 2)
		{
			int stop = 1;
		}
		faces_assoc.clear();
		if (tg->t_edges.size() == 2)
		{
			// merge the only 2 topo edges of vi into a loop edge <-1, -1>
			mergeAssocFaces(faces_assoc, edges_idxes, tg->t_edges);

			// then finalize
			tg->t_vts.clear();
			tg->t_edges.clear();
			tg->t_edges.push_back(make_pair(util::pureLoopEdge(), faces_assoc));
			tg->t_edgeIdxesMap.clear();
			tg->t_edgeIdxesMap[util::pureLoopEdge()].push_back(tg->t_edges.size()-1);
			tg->t_nbVtsMap.clear();
			auto& map = tg->t_origFace2TopoMap;
			for (auto fi_it = faces_assoc.begin(); fi_it != faces_assoc.end(); ++fi_it)
			{
				assert(map.find(*fi_it) == map.end());
				map[*fi_it] = 0;
			}
		}
		else
		{
			// no need to merge. simply finalize.
			// tg.t_vts: computed
			// tg.t_edges: computed
			// tg.t_nbVtsMap: computed
			// tg.t_edgeIdxesMap:
			tg->t_edgeIdxesMap.clear();
			for (unsigned tei = 0; tei < tg->t_edges.size(); ++tei)
			{
				TopoGraph::TEdge& te = tg->t_edges[tei];
				tg->t_edgeIdxesMap[te.first].push_back(tei);
			}
			for (unsigned tei = 0; tei < tg->t_edges.size(); ++tei)
			{
				TopoGraph::TEdge& te = tg->t_edges[tei];
				for (auto fi_it = te.second.begin(); fi_it != te.second.end(); ++fi_it)
				{
					tg->t_origFace2TopoMap[*fi_it] = tei;
				}
			}
		}
		// determine topo type
		TopoType tt = identify(*tg);
		t_type[vi] = tt;
		t_n_sector[vi] = tg->t_edges.size();
		/// Begin: optimization
		// need to set the dependent st vts' # sector and type to be the same of the rep.
		st_onE.clear();
		auto e_idx = m_stSubdiv.getResidingEdgeIdx(vi);
		m_stSubdiv.getStVertIndicesOnTriEdge(e_idx, st_onE);
		for (size_t i = 0; i < st_onE.size(); ++i)
		{
			auto st_vi = st_onE[i];
			// assign valid type & # sectors to *dependent* st vt
			t_n_sector[st_vi] = t_n_sector[vi];
			t_type[st_vi] = tt;
			updateTopoStats(tt);
		}
		/// End: optimization

		// let's save memory here
		if (t_n_sector[vi] == 1)
		{
			del_topo(vi);

			/// Begin: optimization
			// need to also delete dependent vts that refer to this vert's topo graph
			for (size_t i = 0; i < st_onE.size(); ++i)
			{
				auto st_vi = st_onE[i];
				del_topo(st_vi);
			}
			/// End: optimization
		}
	}

	// compact memory usage
	for (auto tg = t_graphs.begin(); tg != t_graphs.end(); ++tg)
		if (*tg)
			(*tg)->shrinkToFit();

#ifdef PRINT_MEM_USAGE
	GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc));
	cout << "mem. usage ( by all topo graphs ) in MB: " << (pmc.WorkingSetSize - memCurUsed) / (1024*1024) << endl;
	memCurUsed = pmc.WorkingSetSize;
#endif // PRINT_MEM_USAGE
}

void SteinerGraph::identifyAll()
{
	resetTopoStats();
	for (unsigned i = 0; i < t_graphs.size(); ++i)
	{
		TopoType tt = identify(*t_graphs[i]);
		t_type[i] = tt;
		updateTopoStats(tt);
	}
}

TopoType SteinerGraph::identify(const TopoGraph& _tg)
{
	if ( _tg.t_edges.empty() )
		return _tg.t_vts.empty() ? ISOLATE : 
		_tg.t_vts.size() == 2 ? MANIFOLD_1D : JUNCTION_1D;

	const map<TriEdge, vector<TopoGraph::EdgeIdx>>& edgeIdxesMap = _tg.t_edgeIdxesMap;
	bool has_cycle = false;
	bool multi = false;

	// manifold: one edge w/o vts
	if (_tg.t_edges.size() == 1 && _tg.t_vts.empty())
		return MANIFOLD_2D;

	// more than one edge
	// has cycle?
	if (edgeIdxesMap.find(util::pureLoopEdge()) != edgeIdxesMap.end())
	{	
		has_cycle = true;
	}
	else
	{
		// spanning tree to determine
		map<int, bool> visited;
		map<int, bool>::iterator vis_itor;
		deque<int> q;

		// more than one iteration if more than one connected components
		for (unsigned vi = 0; vi < _tg.t_vts.size(); ++vi)
		{
			int v = _tg.t_vts[vi];
			if (visited[v])
				continue;

			q.push_back(v);
			while( !q.empty() )
			{
				int cur = q.front();
				q.pop_front();

				vis_itor = visited.find(cur);
				if (
					(vis_itor != visited.end() && vis_itor->second == true) || 
					(edgeIdxesMap.find(util::makeEdge(cur, cur)) != edgeIdxesMap.end())
					)
				{
					has_cycle = true;
					goto LABEL_END_CYCLE_DETECTION;
				}

				visited[cur] = true;
				// add unvisited neighbor vts
				const vector<int>& nbVts = _tg.t_nbVtsMap.find(cur)->second;
				for (unsigned ni = 0; ni < nbVts.size(); ++ni)
				{
					if (visited.find(nbVts[ni]) == visited.end())
						q.push_back(nbVts[ni]);
				}
			}
		}
		has_cycle = false;
	}
LABEL_END_CYCLE_DETECTION:

	for (unsigned tei = 0; tei < _tg.t_edges.size(); ++tei)
	{
		TriEdge e = _tg.t_edges[tei].first;
		multi = edgeIdxesMap.find(e)->second.size() > 1 && util::notLoop(e) && _tg.t_vts.size() == 2;
		if (multi)
			break;
	}

	if (!has_cycle)
		return BOUNDARY_2D;
	if (multi)
		return NONMANIFOLD_2D;
	// lastly
	return JUNCTION_2D;
}

void SteinerGraph::getAssocFaces(int _v, TopoGraph::EdgeIdx _tei, vector<int>& _faces)
{
	if (simpleTopoAt(_v))
	{
		// in this simple case an explicit topo graph is not available
		// but we can easily find the assoc. faces based whether the vertex is 
		// an original or a steiner vertex
		if (_tei != 0)
			return;
		if (m_stSubdiv.isSteinerVert(_v))
		{
			auto containing_edge = m_stSubdiv.getResidingEdge(_v);
			_faces = m_origG->getNbFaces( containing_edge );
		}
		else
		{
			_faces = m_origG->nbFacesOfVert[_v];
		}
	}
	else
	{
		// in this case we need to use the previously stored topo graph
		_faces = t_graphs[_v]->t_edges[_tei].second;
	}
}

void SteinerGraph::resetTopoStats()
{
	t_cnt.clear();
	t_cnt.resize(7);
}

void SteinerGraph::updateTopoStats(TopoType _tt)
{
	++t_cnt[_tt];
}

void SteinerGraph::printTopoStats()
{
	for (unsigned i = 0; i < t_cnt.size(); ++i)
	{
		cout << topo_type_str[i] << ": " << t_cnt[i] << endl;
	}
}

TopoGraph* SteinerGraph::make_topo(int _vi, TopoGraph* _tg_ptr)
{
	t_graphs[_vi] = _tg_ptr;
	t_holds_actual_storage[_vi] = false;
	return t_graphs[_vi];
}

TopoGraph* SteinerGraph::make_topo(int _vi)
{
	if (!t_graphs[_vi])
		t_graphs[_vi] = new TopoGraph;
	t_holds_actual_storage[_vi] = true;
	return t_graphs[_vi];
}
void SteinerGraph::del_topo(int _vi)
{
	if (t_holds_actual_storage[_vi])
	{
		if (t_graphs[_vi])
		{
			delete t_graphs[_vi];
			t_graphs[_vi] = nullptr;
			t_holds_actual_storage[_vi] = false;
		}
	}
	else
	{
		t_graphs[_vi] = nullptr;
	}
}

/*
void SteinerGraph::augmentTopoGraphs()
{
vector<int> nbVts_xu;
vector<int> nbVts_xv;
typedef set<TriEdge, EdgeCompare> EdgeSet;
EdgeSet tedges_xu, tedges_xv;
vector<int> nbVts_u, nbVts_v;
vector<int> common_vts;

// if the edge is not non-manifold, 
// then there is a unique face that contains this edge
// we don't need to store anything;
// else, we map topo edges using the edges' neighbor faces
auto map_topo_for_edge = [&] (const TriEdge& _e)
{
int u = _e[0], v = _e[1];
TopoGraph& gu = t_graphs[u];
TopoGraph& gv = t_graphs[v];

auto& tedge_assoc_u = gu.t_origEdge2TopoMap.find(_e)->second;
auto& tedge_assoc_v = gv.t_origEdge2TopoMap.find(_e)->second;

if (tedge_assoc_u.first == 'e') // case 1. 
{
// unique mapping if assoc.ed with topo edge
/ *st_topoEdgeMap[TEdgeTernary(u, v, tedge_assoc_u.second)] = 
tedge_assoc_u.second;
st_topoEdgeMap[TEdgeTernary(v, u, tedge_assoc_v.second)] = 
tedge_assoc_v.second;* /
}
else // case 2
{
TriEdge orig_e = _e;
if (m_stSubdiv.isSteinerVert(u))
{
orig_e = m_stSubdiv.getResidingEdge(u);
}
else if (m_stSubdiv.isSteinerVert(v))
{
orig_e = m_stSubdiv.getResidingEdge(v);
}

/ * for each face, find the topo edge at each end
and then map them to each other
* /
const auto& nb_faces = m_origG->nbFacesOfEdge.find(orig_e)->second;
// the topo edges at u & v that cover a nb face
TopoGraph::EdgeIdx tei_u, tei_v; 
for (auto it = nb_faces.begin(); it != nb_faces.end(); ++it)
{
const auto& nb_f = m_origG->faces[*it];
bool success = getAssocTEdgeOfVertForFace(u, nb_f, tei_u);
success = getAssocTEdgeOfVertForFace(v, nb_f, tei_v);

TriEdge te_u = gu.t_edges[tei_u].first;
if (v == te_u[0])
gu.t_map[tei_u][0] = tei_v;
else
gu.t_map[tei_u][1] = tei_v;
TriEdge te_v = gv.t_edges[tei_v].first;
if (u == te_v[0])
gv.t_map[tei_v][0] = tei_u;
else
gv.t_map[tei_v][1] = tei_u;
}
}
};

// compute for st edges
for (unsigned ei = 0; ei < m_stSubdiv.sizeOfStEdges(); ++ei)
{
map_topo_for_edge(m_stSubdiv.getStEdge(ei));
}
}*/

// find the topo edge of _v that _fi mapped to
// return: 
// -1 if face is not mapped to any topo edge of _v 
// -2 if topo graphs not constructed yet;
int SteinerGraph::mapTopo(int _v, unsigned _fi) const
{
	if (t_graphs.empty())
		return -2;

	if (simpleTopoAt(_v))
	{
		// in this simple topo case only one sector is available locally
		// so just return index 0
		return 0;
	}

	// else, it's a more complex topo graph with multiple sectors
	TopoGraph* gv = t_graphs[_v];

	// return the topo edge index that associated with face _fi at _v
	auto o2t_itor = gv->t_origFace2TopoMap.find(_fi);
	if(o2t_itor == gv->t_origFace2TopoMap.end())
	{
		return -1;
	}
	return o2t_itor->second;
}

bool SteinerGraph::simpleTopoAt(int _vi) const
{
	return t_n_sector[_vi] == 1;
}

void SteinerGraph::setEdgeWeightFunc(int _edgeWeight_idx)
{
	switch ( _edgeWeight_idx )
	{
	case 0:
		this->m_curEdgeWt = EW_EUCLIDEAN;
		this->edgeWeight = &SteinerGraph::euclidean_edge_weight;
		break;
	case 1:
		this->m_curEdgeWt = EW_BT2BT3DIFF;
		this->edgeWeight = &SteinerGraph::bt2bt3diff_edge_weight;
		break;
	default:
		cerr << "Error: Unrecognized edge weight index "<<_edgeWeight_idx
			<<". euclidean edge weight is used."<<endl;
		this->edgeWeight = &SteinerGraph::euclidean_edge_weight;
		break;
	}
}

void SteinerGraph::burn()
{
	// start timing
#ifdef PROFILE_SPEED
	auto t_start = clock();
#endif

#ifdef _DEBUG_BURN
	unsigned iter_limit = 20;
	unsigned cnt_limit = 0;

	set<int> vts_debug;
	/*vts_debug.insert(13262);
	vts_debug.insert(172864);*/
	/*vts_debug.insert(668);
	vts_debug.insert(373787);*/

	//for (unsigned vi = 0; vi < m_stSubdiv.sizeOfVts(); ++vi)
	//{
	//	const auto& tg = t_graphs[vi];
	//	if (tg.type == JUNCTION)
	//		vts_debug.insert(vi);
	//}

	const char * debug_burn_file = "debug_burn.txt";
	std::ostream& log_os = std::ofstream(debug_burn_file);
	if (!log_os.good())
	{
		cout << "ERROR: debug burn failed. couldn't open " << debug_burn_file << endl;
	}
#endif // _DEBUG_BURN

	final.clear();
	final.resize(m_stSubdiv.sizeOfVts());
	min_tedge.clear();
	min_tedge.resize(m_stSubdiv.sizeOfVts());
	bt2MA_vert_per_sheet.clear();
	bt2MA_vert_per_sheet.resize(m_stSubdiv.sizeOfVts());
	bt2MA_vert.clear();
	bt2MA_vert.resize( m_stSubdiv.sizeOfVts(), infiniteBurnDist() );
	prev_vert.clear();
	prev_vert.resize(m_stSubdiv.sizeOfVts());
	prev_tedge.clear();
	prev_tedge.resize(m_stSubdiv.sizeOfVts());
	burnt.clear();
	burnt.resize(m_stSubdiv.sizeOfVts());
	bt2bt3diffMA_vert.clear();
	bt2bt3diffMA_vert.resize( m_stSubdiv.sizeOfVts() );

	typedef pair<int, float> IntFloatPair;
	struct IntFloatPairComp
	{
		bool operator () (const IntFloatPair& _a, const IntFloatPair& _b) const
		{
			return _a.second > _b.second;
		}
	};
	//priority_queue<IntFloatPair, vector<IntFloatPair>, IntFloatPairComp> q;
	typedef boost::heap::fibonacci_heap<IntFloatPair, boost::heap::compare<IntFloatPairComp> > fib_heap;
	fib_heap q;
	vector<fib_heap::handle_type> hdls; // store handle for each vertex
	hdls.resize(m_stSubdiv.sizeOfVts());
	vector<bool> is_hdl_valid(hdls.size(), false);

	// put boundary vertices/junction vertices (with open topo edge correctly set) to q
	// NOTE: if BurnScheme::STEINER_ONLY, then only add boundary st vts to q.
	assignBoundaryStartDist(bt3MA_vert);
	std::cout << "Done assigning start dist to boundary sheets." << endl;
	for (unsigned vi = 0; vi < m_stSubdiv.sizeOfVts(); ++vi)
	{
		const auto& tg = t_graphs[vi];
		auto has_open_sheet = TopoGraph::hasAnyOpenSheet( vi );
#ifdef _DEBUG_BURN
		if (vts_debug.count(vi))
		{
			log_os << vi<<"'s topo graph:"<<endl;
			if (!simpleTopoAt(vi))
			{
				log_os << tg<<endl;
				log_os << "has open edge? " << has_open_sheet << endl;
			}
			else
			{
				log_os << "simply has one topo edge." << endl;
				log_os << "has open edge? " << false << endl;
			}
		}
#endif

		auto tt = getTopoType(vi);
		if (tt == BOUNDARY_2D || tt == JUNCTION_2D && has_open_sheet )
		{
			if ( 
				getBurnScheme() == ORIGINAL_AND_STEINER ||
				getBurnScheme() == STEINER_ONLY && m_stSubdiv.isSteinerVert( vi ) )
			{
				hdls[ vi ] = q.push( IntFloatPair( vi, bt2MA_vert[ vi ] ) );
				is_hdl_valid[ vi ] = true;
			}

#ifdef _DEBUG_BURN
			if (vts_debug.count(vi))
			{
				log_os << vi<<" added to q with burn dist "<<bt2MA_vert[vi]<<endl;
			}
#endif
		}
	}

#ifdef _DEBUG_BURN
	for (unsigned vi = 0; vi < m_stSubdiv.sizeOfVts(); ++vi)
	{
		if ( vts_debug.count(vi) )
		{
			log_os << "\n initial states of " <<vi<<" --------:"<<endl;
			log_os << "burnDist:"<<bt2MA_vert[vi]<<endl;
			log_os << "dist:";
			std::for_each(bt2MA_vert_per_sheet[vi].begin(), bt2MA_vert_per_sheet[vi].end(), 
				[&log_os](float _d){log_os<<_d<<" ";});
			log_os << endl;
			log_os << *t_graphs[vi] << endl;
			log_os << endl;
		}
	}
#endif

	// start the burning iteration
	vector< std::pair<TopoGraph::EdgeIdx,int> > fnIdxes;
	// save temp tuples: <tri edges, the face they are on, the topo edge that covers the face>
	vector< std::tuple<vector<int>, unsigned, unsigned> > assoc_info;
	vector<int> assoc_faces;

	// debug
	unsigned cntInQueue_burnt = 0;
	unsigned cntInQueue_burnDistInf = 0;
	unsigned cntInQueue_minTedgeFinal = 0;
	unsigned cntInQueue_burnDistInconsistent = 0;

	try {
		while (!q.empty())
		{
			if (q.size() % 100000 == 0 /*1*/) //debug
				cout << "burn q size: " << q.size() << endl;

			IntFloatPair cur = q.top();
			q.pop();
			int curV = cur.first;
			float curD = cur.second;
			is_hdl_valid[curV] = false;

			if ( burnt[curV] ) //debug
				cntInQueue_burnt ++;
			else if ( bt2MA_vert[curV] == numeric_limits<float>::max() )//debug
				cntInQueue_burnDistInf ++;
			else if ( final[curV][min_tedge[curV]] )//debug
				cntInQueue_minTedgeFinal ++;
			else if ( curD != bt2MA_vert[curV] )
				cntInQueue_burnDistInconsistent ++;

#ifdef _DEBUG_BURN
			if (vts_debug.count(curV)/* || 
									 t_graphs[curV].type == BOUNDARY && curD > 0.0f*/
									 ) // debug
			{
				const auto& tg = t_graphs[curV];
				log_os << "\n new iter-------"<<endl;
				log_os << "curV&D: " << curV <<" "<<curD<< ". " << topo_type_str[getTopoType(curV)] << endl;
				log_os << "true D: " << bt2MA_vert[curV] << endl;
				log_os << "burnt or not: " << burnt[curV] << endl;
				log_os << "checking if skipping curV ..." << endl;
			}
#endif

			// skip if this copy of curV is out-of-date
			if (burnt[curV] || 
				bt2MA_vert[curV] == numeric_limits<float>::max() || 
				curD != bt2MA_vert[curV])
			{
				continue;
			}
			//curD = burnDist_vert[curV];

#ifdef _DEBUG_BURN
			if (vts_debug.count(curV)) // debug
			{
				auto& tg = t_graphs[curV];
				log_os << "curV&D: " << curV <<" "<<curD<< ". " << topo_type_str[getTopoType(curV)] << endl;
				log_os << "min topo edge: " << (int)min_tedge[curV] << endl;
				log_os << "that edge finalized? : " << final[curV][min_tedge[curV]] << endl;
			}
#endif

			const TopoGraph* tg = t_graphs[curV];

			// finalize the burnt sheet and open edges
			TopoGraph::EdgeIdx curTEi = min_tedge[curV];
			fnIdxes.clear();
			if (simpleTopoAt(curV))
			{
				// in this case there is only one sector. simply finalize it.
				fnIdxes.push_back(make_pair(curTEi, -1));
				final[curV][0] = true;
			}
			else
			{
				openTopo(*tg, curTEi, final[curV], fnIdxes);
			}
			assert(!fnIdxes.empty());

#ifdef _DEBUG_BURN
			if (vts_debug.count(curV)) // debug
			{
				log_os << "after openning-------" << endl;
				//log_os << "min topo edge: " << (int)curTEi << endl;
				log_os << "prev of min topo edge: " << prev_vert[curV][curTEi] << endl;
				log_os << "prev sheet of min topo edge: " << (int)prev_tedge[curV][curTEi] << endl;
				log_os << "edges finalized: ";
				for (auto itor = fnIdxes.begin(); itor != fnIdxes.end(); ++itor)
					log_os << (int)(itor->first) << ", ";
				log_os << endl;
			}
#endif // _DEBUG_BURN

			// update v.dist(f), v.prev(f) for just finalized edges
			for (unsigned i = 0; i < fnIdxes.size(); ++i)
			{
				// the fn pair indicating <the topo edge f, and whether it's a regular expose sheet>
				const auto& fn_tuple = fnIdxes[i]; 
				TopoGraph::EdgeIdx f = fn_tuple.first;
				bt2MA_vert_per_sheet[curV][f] = curD;
				//if (f != curTEi)
				if ( fn_tuple.second >= 0 )
				{
					// this topo edge is finalized by another topo edge 
					// in the same topo graph
					prev_vert[curV][f] = -2; 
					prev_tedge[curV][f] = (TopoGraph::EdgeIdx)(fn_tuple.second);
				}

#ifdef _DEBUG_BURN
				if (vts_debug.count(curV))//debug
				{
					log_os << curV <<"'s dist["<<(int)f<<"] updated! "<<bt2MA_vert_per_sheet[curV][f]<< endl;
					log_os << curV <<"'s prev_vert["<<(int)f<<"] updated! "<<prev_vert[curV][f]<< endl;
					log_os << curV <<"'s prev_tedge["<<(int)f<<"] updated! "<<(int)prev_tedge[curV][f]<< endl;
					int stop = true;
				}
#endif // _DEBUG_BURN
			}

			// is curV just charred? burnt?
			burnt[curV] = true;
			float minD_nonFinal = numeric_limits<float>::max();
			for (TopoGraph::EdgeIdx ei = 0; ei < t_n_sector[curV]; ++ei)
			{
				if (final[curV][ei])
					continue;
				// still there are edges not final. 
				// update distF & min_tedges:
				// give it higher priority if this edge is a unburnt boundary (open edge)
				burnt[curV] = false;
				if (bt2MA_vert_per_sheet[curV][ei] <= minD_nonFinal)
				{
					minD_nonFinal = bt2MA_vert_per_sheet[curV][ei];
					bt2MA_vert[curV] = bt2MA_vert_per_sheet[curV][ei];
					min_tedge[curV] = ei;

					//if (tg.type == JUNCTION && tg.isTopoEdgeOpen(ei))
					//	break; // preempt other non-boundary edge
				}
			}

#ifdef _DEBUG_BURN
			if (vts_debug.count(curV))//debug
			{
				log_os << "is "<<curV<<" burnt? "<<burnt[curV]<<endl;
			}
#endif // _DEBUG_BURN

			// add curV back in q if just charred. don't propagate
			if (!burnt[curV]/* && minD_nonFinal < numeric_limits<float>::max()*/)
			{
				hdls[curV] = q.push(IntFloatPair(curV, bt2MA_vert[curV])); // TODO: debug this logic...
				is_hdl_valid[curV] = true;
			}
			else
			{
				// propagate from curV to all its nbs
				// prepare all neighbor edges from the assoc faces
				assoc_info.clear();
				TopoGraph::EdgeIdx tei;
				for (unsigned i = 0; i < t_n_sector[curV]; ++i)
				{
					tei = i;
					assoc_faces.clear();
					getAssocFaces(curV, tei, assoc_faces);
					for (auto fi_it = assoc_faces.begin(); fi_it != assoc_faces.end(); ++fi_it)
					{
						assoc_info.push_back(std::make_tuple(vector<int>(), *fi_it, tei));
						m_stSubdiv.getNbVtsInFace( curV, *fi_it, getBurnScheme()==STEINER_ONLY, 
							std::get<0>(assoc_info.back()) );
					}
				}

				// propagate!
				for (unsigned i = 0; i < assoc_info.size(); ++i)
				{
					TopoGraph::EdgeIdx tei_curV = std::get<2>(assoc_info[i]);
					unsigned fi_curV = std::get<1>(assoc_info[i]);
					const auto& nb_vts = std::get<0>(assoc_info[i]);

					for (unsigned ii = 0; ii < nb_vts.size(); ++ii)
					{
						// find topo sheet fi_u that fi_curV is mapped to
						int u = nb_vts[ii];
						int ret_code = mapTopo(u, fi_curV);
						assert(ret_code >= 0);
						TopoGraph::EdgeIdx fi_u = (TopoGraph::EdgeIdx)ret_code;

						// update f' burn distance if not final
#ifdef _DEBUG_BURN
						if ( vts_debug.count(u)||vts_debug.count(curV) )//debug
						{
							log_os << "is final? "<<u<<", topo edge "<<(int)fi_u<<"... "<<final[u][fi_u]<<endl;
						}
#endif

						if ( !final[u][fi_u] )
						{
							// update f' burn distance if bt2 later than bt3
							/*float protect_d = 
							bt3MA_vert.empty() ? 0.0f : bt3MA_vert[u];*/
							
							float tempd = 
								curD + (this->*edgeWeight)(curV, u);
							//curD + trimesh::dist(m_stSubdiv.getStVert(e[0]), m_stSubdiv.getStVert(e[1]));
#ifdef _DEBUG_BURN
							if ( vts_debug.count(u)||vts_debug.count(curV) )//debug
							{
								float protect_d = bt3MA_vert[ u ];
								log_os << "trying to update "<<u<<" "<<topo_type_str[getTopoType(u)]<<", topo edge "<<(int)fi_u
									<<" from "<<curV<<" "<<topo_type_str[getTopoType(curV)]<<" (curD "<<curD<<"), "
									<<"via "<<(int)tei_curV<<" (new dist "<<tempd<<", old dist "<<bt2MA_vert_per_sheet[u][fi_u]<<") ...\n";

								log_os << "will suspend 2d burning (new dist < bt3)? "<< (tempd < protect_d) <<endl;
							}
#endif

							tempd = update_with_protection(u, tempd); //std::max( tempd, protect_d );

							if (tempd < bt2MA_vert_per_sheet[u][fi_u])
							{
								bt2MA_vert_per_sheet[u][fi_u] = tempd;
								prev_vert[u][fi_u] = curV;
								prev_tedge[u][fi_u] = tei_curV;

#ifdef _DEBUG_BURN
								if (vts_debug.count(u)||vts_debug.count(curV))//debug
								{
									log_os << u << "'s dist["<<(int)fi_u<<"] updated! "<< bt2MA_vert_per_sheet[u][fi_u] << ". "
										<< "prev_vert["<<(int)fi_u<<"] updated! "<< prev_vert[u][fi_u] <<". "
										<< "prev_tedge["<<(int)fi_u<<"] updated! "<< (int)prev_tedge[u][fi_u] <<endl;
									int stop = true;
								}
#endif // _DEBUG_BURN

								if (tempd < bt2MA_vert[u])
								{
									bt2MA_vert[u] = tempd;
									min_tedge[u] = fi_u;

									// add this neighbor u to q since not burnt, or update its copy in the queue
									if ( !is_hdl_valid[u] )
									{
										hdls[u] = q.push(IntFloatPair(u, bt2MA_vert[u]));
										is_hdl_valid[u] = true;
									}
									else
									{
										(*hdls[u]).second = bt2MA_vert[u];
										q.increase(hdls[u]);
									}

#ifdef _DEBUG_BURN
									if (vts_debug.count(u)||vts_debug.count(curV))//debug
									{
										log_os << u <<"'s distF updated! "<<bt2MA_vert[u]<< ". "
											<< "min edge updated! " << (int)min_tedge[u]<< endl;
										int stop = true;
									}
#endif // _DEBUG_BURN
								}
							}
							else
							{
#ifdef _DEBUG_BURN
								if (vts_debug.count(u)||vts_debug.count(curV))//debug
								{
									log_os<<"new dist >= old dist. no update.\n";
								}
#endif // _DEBUG_BURN
							}
						}
					}
				} // propagate
			} // end propagating to neighbors

#ifdef _DEBUG_BURN
			if (vts_debug.count(curV)) // debug
			{
				burnt[curV] ? log_os << curV << " burnt. " << endl
					: log_os << curV << " charred. inserted back."<<endl;
				log_os <<"distF "<<bt2MA_vert[curV] << ". "
					<< "min edge "<<(int)min_tedge[curV]<< endl; 
				log_os << "prev_vert/prev_tedge: "<<endl;
				for (TopoGraph::EdgeIdx tei = 0; tei < t_n_sector[curV]; ++tei)
				{
					log_os <<prev_vert[curV][tei]<<"/"<<(unsigned)prev_tedge[curV][tei]<<" ";
				}
				log_os << endl;
			}

			//dynamic_cast<ofstream*>(log_os)->close();
#endif // _DEBUG_BURN
		} // burning iteration

		if ( getBurnScheme() == STEINER_ONLY )
		{
			cout << "recovering bt for orig vts from st vts... " << endl;
			recover_bt_for_orig_vts();
			cout << "recovering bt done." << endl << endl;
		}

		// reset burn time to those vertices that have unburned pieces 
		// (e.g. due to being part of closed pockets)
		for ( auto i = 0; i < m_stSubdiv.sizeOfVts(); ++i )
		{
			if ( !burnt[ i ] )
			{
				bt2MA_vert[ i ] = infiniteBurnDist();
			}
		}

		// compute difference between dist2surf & burn dist for orig. verts
		float min_burn, max_burn, min_d2s, max_d2s, min_diff, max_diff;
		min_burn = min_d2s = min_diff = numeric_limits<float>::max();
		max_burn = max_d2s = max_diff = numeric_limits<float>::min();
		for (unsigned i = 0; i < m_stSubdiv.sizeOfVts(); ++i)
		{
			float d2s = bt3MA_vert.empty() ? 0.0f : bt3MA_vert[i];
			float d_bt2 = bt2MA_vert[ i ];
			if ( d_bt2 < infiniteBurnDist() )
			{
				bt2bt3diffMA_vert[ i ] = std::max( ( d_bt2 - d2s ), 0.0f );
				min_burn = std::min( bt2MA_vert[ i ], min_burn );
				max_burn = std::max( bt2MA_vert[ i ], max_burn );
				min_diff = std::min( bt2bt3diffMA_vert[ i ], min_diff );
				max_diff = std::max( bt2bt3diffMA_vert[ i ], max_diff );
			}
			else
			{
				bt2bt3diffMA_vert[ i ] = infiniteBurnDist();
			}
			min_d2s = std::min(d2s, min_d2s);
			max_d2s = std::max(d2s, max_d2s);
		}
		cout << "min/max burn: " << min_burn<<", "<<max_burn<<endl; // debug
		cout << "min/max d2s: " << min_d2s<<", "<<max_d2s<<endl; // debug
		cout << "min/max diff: " << min_diff<<", "<<max_diff<<endl; // debug

		// update stats
		int unburnt_cnt = 0; // debug
		resetVtsBurnStats();
		for (unsigned vi = 0; vi < m_origG->vts.size()/*st_vts.size()*/; ++vi)
		{
			if (!burnt[vi])
			{
				unburnt_cnt++;
				//cout << "vert not burnt. vid: " << vi << endl;// debug
				//cout << "topo type: " << topo_type_str[t_graphs[vi].type] << endl;// debug
				continue;
			}
			if ( burnt[ vi ] && bt2MA_vert[ vi ] == infiniteBurnDist() )
				cout << "burnt but infinite distF! vid: " << vi << endl;
			min_burnDistOrigVts = std::min(min_burnDistOrigVts, bt2MA_vert[vi]);
			max_burnDistOrigVts = std::max(max_burnDistOrigVts, bt2MA_vert[vi]);
		}
		min_burnDistAllVts = min_burnDistOrigVts;
		max_burnDistAllVts = max_burnDistOrigVts;
		for (unsigned vi = m_stSubdiv.getStartOfStVts(); vi < m_stSubdiv.getEndOfStVts(); ++vi)
		{
			if (!burnt[vi])
			{
				unburnt_cnt++;
				//cout << "vert not burnt. vid: " << vi << endl;// debug
				//cout << "topo type: " << topo_type_str[t_graphs[vi].type] << endl;// debug
				continue;
			}
			if (burnt[vi] && bt2MA_vert[vi] == numeric_limits<float>::max())
				cout << "burnt but infinite distF! vid: " << vi << endl;
			min_burnDistAllVts = std::min(min_burnDistAllVts, bt2MA_vert[vi]);
			max_burnDistAllVts = std::max(max_burnDistAllVts, bt2MA_vert[vi]);
		}

		if (unburnt_cnt > 0 )
		{
			cout << "------------WARNING: unburned vertices------------" << endl;
			cout << "# un-burnt: " << unburnt_cnt << endl << endl; //debug
		}

		cout << "burnt but still in queue count: " << cntInQueue_burnt << endl; //debug
		cout << "burn dist still inf count: " << cntInQueue_burnDistInf << endl; //debug
		cout << "min_tEdge final but still in queue count: " << cntInQueue_minTedgeFinal << endl; //debug
		cout << "burn dist inconsistent in queue count:" << cntInQueue_burnDistInconsistent << endl; //debug
	}
	catch(const std::exception& se)
	{
		std::cerr <<
			"General error: " <<
			se.what() << std::endl;
	}

	// end timing
#ifdef PROFILE_SPEED
	auto t_duration = clock() - t_start;
	cout << "burn() took " << t_duration * 1000.0f / CLOCKS_PER_SEC << " ms."<<endl;
#endif // PROFILE_SPEED

	//// post-checks
	//auto vi_debug = 7149;
	//cout << "sheet-id/pre-vid for " << vi_debug << ": ";
	//for ( auto tei = 0; tei < t_graphs[vi_debug]->t_edges.size(); ++tei )
	//{
	//	cout << "(" << tei << "," << getPrevVert( vi_debug, tei ) << ") ";
	//}
	//cout << endl;

} // SteinerGraph::burn()

void SteinerGraph::triangulate()
{
	// order vts on a face (in ccw or cw)
	vector<unsigned> vts_2_traverse;
	// order vts based on their bt2 value
	// tuple: bt2, vi, relative position in vts2traverse
	vector<std::tuple<float, unsigned, unsigned> > vert_bt2_pairs;
	// pair: vi, and relative pos in vts2traverse
	vector< std::pair<unsigned,unsigned> > vts_sorted_by_bt2;
	// the triangulator
	Triangulation triangulator;

	vector<TriFace> finner_tris_of_f;
	m_MA_finer_tris.clear();

	// triangulate each face
	const auto& faces = m_origG->faces;
	TriEdge e;
	vector<unsigned> stVts_e;
	for (size_t fi = 0; fi < faces.size(); ++fi)
	{
		finner_tris_of_f.clear();
		vts_2_traverse.clear();
		const auto& f = faces[fi];

		// order vts on f consistently
		for (size_t i = 0; i < 3; ++i)
		{
			vts_2_traverse.push_back(f[i]);
			e = util::makeEdge(f[i], f[(i+1)%3]);
			stVts_e.clear();
			m_stSubdiv.getStVertIndicesOnTriEdge(e, stVts_e);
			if ( e[0] != f[i] )
			{// need to reverse st-vts-on-edge
				vts_2_traverse.insert(vts_2_traverse.end(), stVts_e.rbegin(), stVts_e.rend());
			}
			else
			{
				vts_2_traverse.insert(vts_2_traverse.end(), stVts_e.begin(), stVts_e.end());
			}
		}

		// order vts based on their bt2
		vert_bt2_pairs.clear();
		vts_sorted_by_bt2.clear();
		for ( size_t i = 0; i < vts_2_traverse.size(); ++i )
		{
			auto vi = vts_2_traverse[i];
			auto tei = mapTopo(vi, fi);
			assert( tei >= 0 );
			auto bt2val = this->bt2MA_vert_per_sheet[vi][tei];
			vert_bt2_pairs.push_back(std::make_tuple(bt2val, vi, i));
		}
		std::sort( vert_bt2_pairs.begin(), vert_bt2_pairs.end(), 
			[] (std::tuple<float,unsigned,unsigned> _a, std::tuple<float,unsigned,unsigned> _b) {
				return std::get<0>(_a) < std::get<0>(_b);
		} );

		for (auto it = vert_bt2_pairs.begin(); it != vert_bt2_pairs.end(); ++it)
			vts_sorted_by_bt2.push_back( std::make_pair(std::get<1>(*it), std::get<2>(*it)) );
		triangulator.triangulateChordBased( vts_2_traverse, vts_sorted_by_bt2, finner_tris_of_f );
		for (auto fine_it = finner_tris_of_f.begin(); fine_it != finner_tris_of_f.end(); ++fine_it)
		{
			m_MA_finer_tris.push_back( std::make_pair(*fine_it, fi) );
		}
	}
}

//void SteinerGraph::recordDualwrtStEdge(unsigned _fi, TriEdge _st_e, unsigned _e_dual_vi, unsigned _poly_dual_vi)
//{
//	assert(!m_stEdgeIdx_dual_map.empty());
//
//	// input: _st_e is a steiner edge that contains a dual point
//	// i.e. _dual_vi, which is connected to a poly face dual, _polly_dual_vi
//	auto eid = m_stSubdiv.getStEdgeIndex(util::makeEdge(_st_e[0],_st_e[1]));
//	if (!(eid >= 0))
//	{
//		cout << "Fatal error: invalid edge id!" << endl;
//		exit(-1);
//	}
//	m_stEdgeIdx_dual_map[eid] = _e_dual_vi;
//
//	auto eDual_face_pr = make_pair(_e_dual_vi, _fi);
//	//assert(m_eDualFaceIdPr_polyDual_map.find(eDual_face_pr) == m_eDualFaceIdPr_polyDual_map.end());
//	m_eDualFaceIdPr_polyDual_map[eDual_face_pr] = _poly_dual_vi;
//}

void SteinerGraph::burnMedialCurveNetwork(
	bool _start_time,
	bool _stop_at_junction,
	bool _protect_bt2,
	bool _under_estimate_dist)
{
	// related states of burning
	bt1_medialCurve.resize( dual_vts.size() );
	bt1_medialCurve.assign( dual_vts.size(), std::numeric_limits<float>::max() );
	burnPrev_medialCurve.resize( dual_vts.size() );
	burnPrev_medialCurve.assign( dual_vts.size(), -1 );
	burnNext_medialCurve.resize( dual_vts.size() );
	burnNext_medialCurve.assign( dual_vts.size(), -1 );

	// construct neighbor-ship
	vector< vector<unsigned> > nbVts(dual_vts.size());
	for (unsigned dei = 0; dei < dual_edges.size(); ++dei)
	{
		TriEdge e = dual_edges[dei];
		nbVts[e[0]].push_back(e[1]);
		nbVts[e[1]].push_back(e[0]);
	}

	// count initial degrees of each vert
	vector<unsigned> remain_degs(dual_vts.size());
	for (unsigned vi = 0; vi < dual_vts.size(); ++vi)
	{
		remain_degs[vi] = nbVts[vi].size(); 
	}
	typedef std::pair<unsigned, float> uint_float_pair;
	struct uint_float_comp
	{
		bool operator () (const uint_float_pair& _a, const uint_float_pair& _b)
		{
			return _a.second > _b.second;
		}
	};
	
	std::priority_queue< uint_float_pair,vector<uint_float_pair>,uint_float_comp > q; // containing <firefront vert, its dist>
	
	// init q and start time of boundary nodes
	vector<unsigned int> st_vts_on_e;
	for (unsigned vi = 0; vi < dual_vts.size(); ++vi)
	{
		//if (!_under_estimate_dist)
		{
			/* BEGIN: orig. burning code */
			if ( nbVts[ vi ].size() == 1 )
			{
				bool is_on_pocket = false;
				//if ( isEdgeDual( vi ) >= 0 )
				//{
				//	auto e = m_stSubdiv.getOrigTriEdge( isEdgeDual( vi ) );
				//	// it's ok to only inspect the first st vert on edge
				//	is_on_pocket = !burnt[ m_stSubdiv.getClosestVertOnEdge( e[ 0 ], e[ 1 ], 0 ) ]; 
				//}
				if ( !is_on_pocket )
				{
					q.push( std::make_pair( vi, bt2_MC[ vi ] ) );
					if ( _start_time )
					{
						bt1_medialCurve[ vi ] = bt2_MC[ vi ];
						burnPrev_medialCurve[ vi ] = vi;
					}
					else
					{
						bt1_medialCurve[ vi ] = 0.0f;
						burnPrev_medialCurve[ vi ] = vi;
					}
				}
			}
			/* END: orig. burning code */
		}
		//else
		//{
		//	/* BEGIN: modified burning with under-estimating burn distance */
		//	if ( nbVts[vi].size() == 1 )
		//	{
		//		if (this->is_face_dual[vi] >= 0)
		//		{
		//			// use its next ( an edge dual ) as starting bndry point
		//			unsigned e_dual = nbVts[vi][0];
		//			burnNext_medialCurve[vi] = e_dual;
		//			-- remain_degs[vi];
		//			-- remain_degs[e_dual];

		//			if (remain_degs[e_dual] == 1)
		//			{
		//				if (_protect_bt2)
		//				{
		//					bt1_medialCurve[e_dual] = std::max(
		//						bt2_MC[vi] + trimesh::dist(dual_vts[vi], dual_vts[e_dual]), 
		//						bt2_MC[e_dual]
		//					);
		//				}
		//				else
		//				{
		//					bt1_medialCurve[e_dual] = 
		//						bt2_MC[vi] + trimesh::dist(dual_vts[vi], dual_vts[e_dual]);
		//				}
		//				burnPrev_medialCurve[e_dual] = vi;
		//				q.push(std::make_pair(e_dual, bt1_medialCurve[e_dual]));
		//			}
		//		}
		//		else
		//		{
		//			q.push(std::make_pair(vi, bt2_MC[vi]));
		//		}
		//	}
		//	/* END: modified burning with under-estimating burn distance */
		//}
	}
	// burning iterations start
	while ( !q.empty() )
	{
		auto cur_item = q.top();
		q.pop();
		unsigned cur_v = cur_item.first;
		float cur_d = cur_item.second;
		if (remain_degs[cur_v] == 0) // skip if burnt
			continue;

		assert(remain_degs[cur_v] == 1);

		/* modified burning: must be an edge dual */
		//if (_under_estimate_dist)
		//	assert( this->is_face_dual[cur_v] == -1 );

		const auto& nbs = nbVts[cur_v];
		int last_nb = -1;
		for (auto nb_it = nbs.begin(); nb_it != nbs.end(); ++nb_it)
		{
			if (_stop_at_junction && remain_degs[*nb_it] > 2)
				continue;
			if (remain_degs[*nb_it] > 0)
			{
				last_nb = (int)(*nb_it);
				break;
			}
		}

		burnNext_medialCurve[cur_v] = last_nb;
		if (last_nb == -1)
			continue;

		-- remain_degs[cur_v];
		-- remain_degs[last_nb];

		//if (!_under_estimate_dist)
		{
			// orig burning algorithm
			if ( remain_degs[last_nb] == 1 )
			{
				if (_protect_bt2)
					bt1_medialCurve[last_nb] = std::max(
					bt1_medialCurve[cur_v] + trimesh::dist(dual_vts[cur_v], dual_vts[last_nb]),
					bt2_MC[last_nb]);
				else
					bt1_medialCurve[last_nb] = 
					bt1_medialCurve[cur_v] + trimesh::dist(dual_vts[cur_v], dual_vts[last_nb]);
				burnPrev_medialCurve[last_nb] = cur_v;
				q.push( std::make_pair(last_nb, bt1_medialCurve[last_nb]) );
			}
		}
		//else
		//{
		//	/* BEGIN: modified burning with under-estimating burn distance */
		//	// last_nb must be a face dual
		//	assert(this->is_face_dual[last_nb] >= 0);

		//	if (remain_degs[last_nb] == 1)
		//	{
		//		// find the last nb of last_nb, 
		//		// i.e. the last edge dual of the face corresponding to the face dual last_nb
		//		const auto& nbs = nbVts[last_nb];
		//		int last_eDual = -1;
		//		for (auto nb_it = nbs.begin(); nb_it != nbs.end(); ++nb_it)
		//		{
		//			if (_stop_at_junction && remain_degs[*nb_it] > 2)
		//				continue;
		//			if (remain_degs[*nb_it] > 0)
		//			{
		//				last_eDual = (int)(*nb_it);
		//				break;
		//			}
		//		}
		//		this->burnNext_medialCurve[last_nb] = last_eDual;
		//		if (last_eDual == -1)
		//			continue;

		//		// now last edge dual located. 
		//		-- remain_degs[last_eDual];
		//		-- remain_degs[last_nb];

		//		// compute the BT1 output by the face corresponding to last_nb
		//		float max_BT1 = -1.0f;
		//		for (auto eDual_it = nbs.begin(); eDual_it != nbs.end(); ++eDual_it)
		//		{
		//			if (*eDual_it == last_eDual)
		//				continue;
		//			max_BT1 = std::max(max_BT1, 
		//				bt1_medialCurve[*eDual_it] + trimesh::dist(dual_vts[*eDual_it], dual_vts[last_eDual])
		//				);
		//		}

		//		last_eDual;
		//		auto findIt_lastEdgeDual = bt2_MC_perFace.find(trimesh::ivec2(last_eDual, is_face_dual[last_nb]));
		//		assert(findIt_lastEdgeDual != bt2_MC_perFace.end());
		//		float bt2_lastEdgeDual = findIt_lastEdgeDual->second;
		//		auto findIt_lastNeighbor = bt2_MC_perFace.find(trimesh::ivec2(last_nb, is_face_dual[last_nb]));
		//		float bt2_lastNeighbor = findIt_lastNeighbor->second;
		//		float max_BT2 = std::max(bt2_lastNeighbor, bt2_lastEdgeDual);
		//		// change the pos of the face dual to where the last edge dual is
		//		// so that the underestimated burn dist is realized
		//		dual_vts[last_nb] = dual_vts[last_eDual];

		//		findIt_lastNeighbor->second = bt2_lastEdgeDual;
		//		bt2_MC[last_nb] = findIt_lastNeighbor->second; //bt2_MC[last_eDual];
		//		bt3_MC[last_nb] = bt3_MC[last_eDual];
		//		if (_protect_bt2)
		//			bt1_medialCurve[last_nb] = std::max(max_BT1, max_BT2);
		//		else
		//			bt1_medialCurve[last_nb] = max_BT1;

		//		if (remain_degs[last_eDual] == 1)
		//		{
		//			if (_protect_bt2)
		//				bt1_medialCurve[last_eDual] = std::max(max_BT1, max_BT2);
		//			else
		//				bt1_medialCurve[last_eDual] = max_BT1;
		//			this->burnPrev_medialCurve[last_eDual] = last_nb;
		//			q.push(std::make_pair(last_eDual, max_BT1));
		//		}
		//	}
		//	/* END: modified burning with under-estimating burn distance */
		//}
	} // end of burn iterations

	if ( _under_estimate_dist )
	{
		// if under_estimate is specified, that means face duals are re-located 
		// during this burning. Need to perform another burning on the modified dual structure.
		// Note: to prevent recursion, we DO NOT specify under_estimate.
		
		//this->burnMedialCurveNetwork(true, false, _protect_bt2, false);
		relocate_face_dual();
		this->burnMedialCurveNetwork(true, false, _protect_bt2, !_under_estimate_dist);
	}

	unsigned cnt_remainBoundary = 0;
	unsigned cnt_unburntVts = 0;
	for (unsigned vi = 0; vi < dual_vts.size(); ++vi)
	{
		if (remain_degs[vi] == 1)
			cnt_remainBoundary ++;
		if (remain_degs[vi] > 0)
			cnt_unburntVts ++;
	}
	cout << "# remaining network boundary points: " << cnt_remainBoundary << endl;
	cout << "# unburnt network points: " << cnt_unburntVts << endl;
}

void SteinerGraph::relocate_face_dual()
{
	// for each poly/face dual, find the last burned-away edge dual
	// Note: skip the face dual that doesn't have a next-burn-away dual
	// i.e. this face dual is the last burn-away dual in its poly region
	for (unsigned i = 0; i < dual_vts.size(); ++i)
	{
		if ( !(m_is_face_dual[i] >= 0) )
			continue;

		auto next_dual = burnNext_medialCurve[i];
		if (next_dual == -1)
			continue;

		// since poly dual will be located, need to adjust bt2 bt3 accordingly
		auto findIt_lastEdgeDual = m_bt2_MC_perFace.find(trimesh::ivec2(next_dual, m_is_face_dual[i]));
		assert(findIt_lastEdgeDual != m_bt2_MC_perFace.end());
		float bt2_lastEdgeDual = findIt_lastEdgeDual->second;
		auto findIt_lastNeighbor = m_bt2_MC_perFace.find(trimesh::ivec2(i, m_is_face_dual[i]));
		float bt2_lastNeighbor = findIt_lastNeighbor->second;
		float max_BT2 = std::max(bt2_lastNeighbor, bt2_lastEdgeDual);
		findIt_lastNeighbor->second = bt2_lastEdgeDual;
		bt2_MC[i] = findIt_lastNeighbor->second; //bt2_MC[last_eDual];
		bt3_MC[i] = bt3_MC[next_dual];

		// relocate poly dual
		dual_vts[i] = dual_vts[next_dual];
	}
}

// build neighbor-ship given the vts of edges of a line network
void SteinerGraph::build_neighborship(
	const vector<TriPoint>& _vts, const vector<TriEdge>& _edges, 
	vector<vector<unsigned>>& _neighbors, vector<int>& _nbs_cnt)
{
	// build up neighbor-ship
	_neighbors.resize(_vts.size());
	_nbs_cnt.resize(_vts.size(), 0);
	for (unsigned ei = 0; ei < _edges.size(); ++ei)
	{
		const auto& e = _edges[ei];
		_neighbors[e[0]].push_back(e[1]);
		_neighbors[e[1]].push_back(e[0]);
		_nbs_cnt[e[0]] ++;
		_nbs_cnt[e[1]] ++;
	}
}

void SteinerGraph::perturbVts(vector<TriPoint>& _pert_vts, float _pert_amount)
{
	_pert_vts.resize(m_stSubdiv.sizeOfVts());

	// first perturb the orig triangle vts
	TriPoint v;
	trimesh::vec perturb_vec;
	for (unsigned vi = 0; vi < m_origG->vts.size(); ++vi)
	{
		v = m_stSubdiv.getVert(vi);

		perturb_vec = trimesh::normalize( 
			trimesh::vec(rand() % 100000, rand() % 100000, rand() % 100000) );
		perturb_vec = _pert_amount * perturb_vec;
		v += perturb_vec;

		_pert_vts[vi] = v;
	}

	// second reconstruct the steiner vts on each edge that's perturbed
	TriEdge e;
	vector<unsigned> st_vts_onE;
	for (unsigned ei = 0; ei < m_origG->edges.size(); ++ei)
	{
		e = m_origG->edges[ei];
		st_vts_onE.clear();
		m_stSubdiv.getStVertIndicesOnTriEdge(ei, st_vts_onE);

		unsigned n_stVtsOnE = st_vts_onE.size();
		TriPoint new_stV;
		// replace old steiner with new steiner vts
		for (int i = 0; i < n_stVtsOnE; ++i)
		{
			float interp_weight = float(i+1) / (n_stVtsOnE+1);
			new_stV = trimesh::mix(_pert_vts[e[0]], _pert_vts[e[1]], interp_weight);
			_pert_vts[st_vts_onE[i]] = new_stV;
		}
	}
}

void SteinerGraph::debug_dualline_cycles( 
	const unsigned _vts_size, 
	const vector<TriEdge>& _dual_edges, vector<int>& _from_which_face_dual_edge,
	bool _output_smallest_cycle /* = false */, 
	std::ostream* _os /* = std::cout */ )
{
	/* detect all cycle bases in this graph */
	cout << "detect all cycle bases in this graph." << endl;
	typedef list<vector<unsigned> > CycleList;
	CycleList cycles; 	// store actual cycle bases here

	find_cycle_base(_vts_size, _dual_edges, cycles);

	cout << cycles.size() <<" cycles found." << endl;

	if (_os == nullptr)
		return;

	cout << "output cycle base to ostream..." << endl;
	/* now a cycle base has been found. print it to _os */
	map<TriEdge, unsigned> dualE_face_map;
	for (unsigned ei = 0; ei < _dual_edges.size(); ++ei)
		dualE_face_map[_dual_edges[ei]] = _from_which_face_dual_edge[ei];

	unsigned cnt_cycle = 0;
	set<unsigned> face_set;
	for (auto c_it = cycles.begin(); c_it != cycles.end(); ++c_it)
	{
		*_os << "------dual cycle "<<cnt_cycle++<<"------:"<<endl;
		for (unsigned vi = 0; vi < c_it->size(); ++vi)
			*_os << (*c_it)[vi] << " ";
		*_os << endl;

		face_set.clear();
		for (unsigned vi = 0; vi < c_it->size()-1; ++vi)
		{
			TriEdge e = util::makeEdge((*c_it)[vi], (*c_it)[vi+1]);
			auto dualE_face_iter = dualE_face_map.find(e);
			assert(dualE_face_iter != dualE_face_map.end());
			face_set.insert(dualE_face_iter->second);
		}
		for (auto f_it = face_set.begin(); f_it != face_set.end(); ++f_it)
			*_os << *f_it << " ";
		*_os << endl;
	}

	// clean up
	cycles.clear();
	face_set.clear();
	dualE_face_map.clear();
}

void SteinerGraph::debug_burnline_cycles( 
	const vector<TriPoint>& _burn_vts, 
	const vector<std::pair<TriEdge, int> >& _burntE_assocF_list, 
	bool _output_smallest_cycle /* = false */, 
	std::ostream* _os /* = std::cout */ )
{
	map<TriEdge, vector<int> > burntE_assocF_map;
	vector<TriEdge> burntE_list;
	for (auto it = _burntE_assocF_list.begin(); it != _burntE_assocF_list.end(); ++it)
	{
		TriEdge e = it->first;
		e = util::makeEdge(e[0], e[1]);
		unsigned pre_size = burntE_assocF_map.size();
		burntE_assocF_map[e].push_back(it->second);
		if (burntE_assocF_map.size() > pre_size)
			burntE_list.push_back(e);
	}

	list<vector<unsigned> > cycles;
	find_cycle_base(_burn_vts.size(), burntE_list, cycles);
	cout << "find_cycle_base() finished." << endl;

	if (_os != nullptr)
	{
		/* now a cycle base has been found. print it to _os */
		unsigned cnt_cycle = 0;
		for (auto c_it = cycles.begin(); c_it != cycles.end(); ++c_it)
		{
			*_os << "------burn cycle "<<cnt_cycle++<<"------:"<<endl;
			for (unsigned vi = 0; vi < c_it->size(); ++vi)
				*_os << (*c_it)[vi] << " ";
			*_os << endl;

			for (unsigned vi = 0; vi < c_it->size()-1; ++vi)
			{
				TriEdge e = util::makeEdge((*c_it)[vi], (*c_it)[vi+1]);
				auto burntE_assocF_iter = burntE_assocF_map.find(e);
				assert(burntE_assocF_iter != burntE_assocF_map.end());
				for (auto f_it = burntE_assocF_iter->second.begin(); 
					f_it != burntE_assocF_iter->second.end(); ++f_it)
					*_os << *f_it << " ";
			}
			*_os << endl;
		}
	}

	// clean-up
	cycles.clear();
	burntE_assocF_map.clear();
	burntE_list.clear();
	
} // debug_burnline_cycles

void SteinerGraph::debug_dualline_conn_comp( 
	const unsigned _vts_size, const vector<TriEdge>& _dual_edges, 
	vector<int>& _from_which_face_dual_edge,
	bool _output_smallest_cc /* = false */, 
	std::ostream* _os /* = std::cout */ )
{
	cout << "debugging dual c.c. ..."<<endl;

	/* find all c.c.s */
	vector<vector<TriEdge> > cc_eList; // conn. cmpnt.s edge list
	vector<vector<unsigned> > cc_vList; // conn. compts vts list
	find_conn_cmpnts(_vts_size, _dual_edges, cc_eList, cc_vList);

	cout << "done! "<< cc_eList.size()<<" components found."<<endl;

	// sort c.c. according to the cardinality
	vector<unsigned> cc_order(cc_eList.size());
	for (unsigned o = 0; o < cc_order.size(); ++o)
		cc_order[o] = cc_eList[o].size();

	std::sort(cc_order.begin(), cc_order.end());

	if (_os == nullptr)
		return;

	/* output to _os */
	// create a mapping: dual edge -> the face that it's from
	map<TriEdge, unsigned> dualE_face_map;
	for (unsigned ei = 0; ei < _dual_edges.size(); ++ei)
		dualE_face_map[_dual_edges[ei]] = _from_which_face_dual_edge[ei];

	unsigned cnt_cc = 0;
	set<unsigned> face_set;
	for (unsigned o = 0; o < cc_order.size(); ++o)
	{
		const auto& cc_vts = cc_vList[o];
		const auto& cc_edges = cc_eList[o];

		*_os << "------cc "<<cnt_cc++<<"------:"<<endl;
		if (cnt_cc == 1 && cc_edges.size() > 1000)
		{
			*_os << "skip. major component too big." << endl;
			continue;
		}

		for (auto v_it = cc_vts.begin(); v_it != cc_vts.end(); ++v_it)
		{
			*_os << *v_it << " ";
		}
		*_os << endl;
		for (auto e_it = cc_edges.begin(); e_it != cc_edges.end(); ++e_it)
			*_os << *e_it << " ";
		*_os << endl;
		// the faces from which this c.c. comes from
		face_set.clear();
		for (auto e_it = cc_edges.begin(); e_it != cc_edges.end(); ++e_it)
		{
			auto dualE_face_iter = dualE_face_map.find(*e_it);
			assert(dualE_face_iter != dualE_face_map.end());
			face_set.insert(dualE_face_iter->second);
		}
		for (auto f_it = face_set.begin(); f_it != face_set.end(); ++f_it)
			*_os << *f_it << " ";

		*_os << endl;
	}

	cout << cc_eList.size()<<" components output to file."<<endl;

	// clean-up
	cc_eList.clear(); cc_eList.shrink_to_fit();
	cc_vList.clear(); cc_vList.shrink_to_fit();
	dualE_face_map.clear();
	face_set.clear();
} // debug_dualline_conn_comp

void SteinerGraph::find_cycle_base( 
	const unsigned _network_vts_size,
	const vector<TriEdge>& _network_edges, 
	list<vector<unsigned> >& _cycles )
{
	/* build adjacency list for the line network */
	vector<vector<unsigned> > adj_list_dual(_network_vts_size);
	for (unsigned ei = 0; ei < _network_edges.size(); ++ei)
	{
		TriEdge e = _network_edges[ei];
		adj_list_dual[e[0]].push_back(e[1]);
		adj_list_dual[e[1]].push_back(e[0]);
	}

	/* detect all cycle bases in this graph */

	/* 1. find a spanning tree & the edges that are not part of the tree 
	Note: there is a tree for each connected component of the graph 
	*/
	// during DFS, whether the dual vert is visited or not 
	vector<bool> visited(_network_vts_size, false);
	// prev of each dual vert in the spanning tree
	vector<unsigned> prev_in_tree(_network_vts_size);
	// edges of the spanning tree
	set<TriEdge> edges_in_tree;
	// the stack used in DFS
	std::stack<unsigned> s;
	// perform DFS for each connected component
	for (unsigned vi = 0; vi < _network_vts_size; ++vi)
	{
		if (visited[vi])
			continue;
		// else, found a new connected component
		// start DFS with vi as the root

		s.push(vi);
		visited[vi] = true;
		prev_in_tree[vi] = vi;
		while ( !s.empty() )
		{
			unsigned cur_parent = s.top();
			const auto& nbs = adj_list_dual[cur_parent];
			unsigned nbi = 0;

			if (cur_parent == 43 || cur_parent == 52) //debug
				int stop = 1;

			for (; nbi < nbs.size(); ++nbi)
			{
				unsigned nb = nbs[nbi];
				if ( !visited[nb] )
				{
					s.push(nb);
					visited[nb] = true;
					prev_in_tree[nb] = cur_parent;
					TriEdge tree_e = util::makeEdge(cur_parent, nb);
					edges_in_tree.insert(tree_e);
					break;
				}
			}
			if ( nbi == nbs.size() ) // all children of cur parent have been visited
			{
				s.pop();
			}
		}
	}

	/* 2. detect cycles by keeping track of what happens when each outside edge is added */
	auto find_path = [&prev_in_tree] (unsigned _u, unsigned _v, vector<unsigned>& _path) 
	{
		_path.clear();

		/* find the lowest common ancestor (LCA) */
		// find path root -> _u
		vector<unsigned> u_path;
		u_path.push_back(_u);
		while (prev_in_tree[_u] != _u)
		{
			u_path.push_back(prev_in_tree[_u]);
			_u = prev_in_tree[_u];
		}
		// find path root -> _v;
		vector<unsigned> v_path;
		v_path.push_back(_v);
		while (prev_in_tree[_v] != _v)
		{
			v_path.push_back(prev_in_tree[_v]);
			_v = prev_in_tree[_v];
		}
		// LCA is the last common node between u/v_path
		unsigned min_len = std::min(u_path.size(), v_path.size());
		auto u_rit = u_path.rbegin();
		auto v_rit = v_path.rbegin();

		if (*u_rit != *v_rit)
			return;

		for (; u_rit != u_path.rend() && v_rit != v_path.rend(); ++u_rit, ++v_rit)
		{
			if ( *u_rit != *v_rit )
			{
				break;
			}
		}

		/* now u_rit & v_rit both point to LCA. make path: u -> LCA -> v */
		if (u_rit == u_path.rend()) // u & v is on same branch and u_path is shorter
		{
			v_rit--;
			_path.insert(_path.end(), v_rit, v_path.rend());
		}
		else if (v_rit == v_path.rend()) // u & v is on same branch and v_path is shorter
		{
			u_rit--;
			_path.insert(_path.end(), u_rit, u_path.rend());
		}
		else // u & v on different branch. LCA is different from both of u v
		{
			_path.insert(_path.end(), u_path.begin(), u_rit.base()+1);
			_path.insert(_path.end(), v_rit, v_path.rend());
		}
	};
	// for each edge outside spanning tree, there is a new base cycle to add
	vector<unsigned> new_cycle;
	for (unsigned ei = 0; ei < _network_edges.size(); ++ei)
	{
		TriEdge e = _network_edges[ei];
		if ( edges_in_tree.count(e) == 0 )
		{
			// a new base cycle
			new_cycle.clear();
			find_path(e[0], e[1], new_cycle); // now new_cycle contains path e[0]->e[1]
			if (new_cycle.empty())
				continue;

			if (new_cycle.back() == e[1])
				new_cycle.push_back(e[0]);
			else
				new_cycle.push_back(e[1]);

			_cycles.push_back(new_cycle);
		}
	}

	// clean up
	visited.clear(); visited.shrink_to_fit();
	prev_in_tree.clear(); prev_in_tree.shrink_to_fit();
	edges_in_tree.clear();
	new_cycle.clear(); new_cycle.shrink_to_fit();
	adj_list_dual.clear(); adj_list_dual.shrink_to_fit();
} // find_cycle_base

void SteinerGraph::find_conn_cmpnts(
	const unsigned _vts_size, 
	const vector<vector<int>>& _adj_list, 
	vector<vector<TriEdge> >& _cc_eList, 
	vector<vector<unsigned> >& _cc_vList
	)
{
	/* actually identify all conn. cmpnt.s */
	vector<bool> visited(_vts_size, false);
	std::queue<unsigned> q;
	vector<TriEdge> cc_edges;
	set<unsigned> cc_vts;
	for (unsigned vi = 0; vi < _vts_size; ++vi)
	{
		if (visited[vi])
			continue;

		// new c.c. starts from vi
		cc_edges.clear();
		cc_vts.clear();
		visited[vi] = true;
		q.push(vi);
		cc_vts.insert(vi);
		// start tracing new c.c.
		while (!q.empty())
		{
			unsigned cur = q.front();
			q.pop();

			const auto& nbs = _adj_list[cur];
			for (auto nb_it = nbs.begin(); nb_it != nbs.end(); ++nb_it)
			{
				if (!visited[*nb_it])
				{
					q.push(*nb_it);
					visited[*nb_it] = true;
					cc_edges.push_back(util::makeEdge(cur, *nb_it));
					cc_vts.insert(*nb_it);
				}
			}
		} // tracing new c.c.

		if (!cc_vts.empty())
		{
			_cc_eList.push_back(cc_edges);
			_cc_vList.push_back(vector<unsigned>(cc_vts.begin(), cc_vts.end()));
		}
	}

	// clean-up
	visited.clear(); visited.shrink_to_fit();
	cc_edges.clear(); cc_edges.shrink_to_fit();
	cc_vts.clear();
}

void SteinerGraph::find_conn_cmpnts( 
	const unsigned _vts_size, 
	const vector<TriEdge>& _network_edges, 
	vector<vector<TriEdge> >& _cc_eList, 
	vector<vector<unsigned> >& _cc_vList 
	)
{
	/* build adjacency list for the line network */
	vector<vector<int> > adj_list(_vts_size);
	for (unsigned ei = 0; ei < _network_edges.size(); ++ei)
	{
		TriEdge e = _network_edges[ei];
		adj_list[e[0]].push_back(e[1]);
		adj_list[e[1]].push_back(e[0]);
	}

	find_conn_cmpnts(_vts_size, adj_list, _cc_eList, _cc_vList);

	adj_list.clear(); adj_list.shrink_to_fit();
} // find_conn_cmpnts

bool SteinerGraph::outputToWenping(const char* _file_name, const vector<float>& _radii)
{
	ofstream ofile(_file_name);
	if (!ofile.is_open())
	{
		cout << "Error: couldn't open file "<< _file_name<<endl;
		return false;
	}

	// output # vts, edges, faces
	ofile << m_origG->vts.size() << " "
		<< m_origG->edges.size() << " "
		<< m_origG->faces.size() << endl;

	// output vts, one on each line
	for (unsigned vi = 0; vi < m_origG->vts.size(); ++vi)
	{
		auto v = m_origG->vts[vi];
		auto r = _radii[vi];

		ofile << "v " << setiosflags(ios::fixed) 
			<< v[0] <<" "<<v[1]<<" "<<v[2]<<" "
			<< r << endl;
	}

	// output edges, one on each line
	for (unsigned ei = 0; ei < m_origG->edges.size(); ++ei)
	{
		auto e = m_origG->edges[ei];

		ofile << "e "
			<< e[0]<<" "<<e[1]<<endl;
	}

	// output faces, one on each face
	for (unsigned fi = 0; fi < m_origG->faces.size(); ++fi)
	{
		auto f = m_origG->faces[fi];

		ofile << "f "
			<<f[0]<<" "<<f[1]<<" "<<f[2]<<endl;
	}

	ofile.close();

	return true;
}

// export bt to the specified file
void SteinerGraph::exportPerVertexET(std::string _file_name)
{
	// immediately return if bt is not available yet.
	if (bt2bt3diffMA_vert.size() < m_origG->vts.size())
		return;

	ofstream out_file(_file_name);
	// encode little-endianess.
	//out_file << (is_little_endian() ? 1 : 0) << endl;
	// num of original vts
	out_file << m_origG->vts.size() << endl;
	std::streamsize ss = out_file.precision();
	auto flgs = out_file.flags();
	out_file.precision(numeric_limits<float>::max_digits10);
	for (size_t i = 0; i < m_origG->vts.size(); ++i)
	{
		// vertex-id its-burntime
		out_file << i <<" " << bt2bt3diffMA_vert[i]<<endl;
	}
	out_file.precision(ss);
	out_file.setf(flgs);
	out_file.close();
}
void SteinerGraph::exportPerSectorET(std::string _file_name)
{
	// immediately return if bt is not available yet.
	if (bt2MA_vert_per_sheet.size() < m_origG->vts.size())
		return;

	ofstream out_file(_file_name);
	// encode little-endianess.
	//out_file << (is_little_endian() ? 1 : 0) << endl;
	// num of original vertices and num of original faces
	out_file << m_origG->vts.size() <<" "<< m_origG->faces.size() << endl;
	
	// for each vert: the bt on each sector, 
	// and for each sector the associated faces' ids
	std::streamsize ss = out_file.precision();
	auto flgs = out_file.flags();
	out_file.precision(numeric_limits<float>::max_digits10);
	vector<int> assoc_faces;
	for (size_t vi = 0; vi < m_origG->vts.size(); ++vi)
	{
		// num of sector
		out_file << t_n_sector[vi] <<endl;
		// ET on each sector
		auto r = this->bt3MA_vert[vi];
		for (size_t si = 0; si < t_n_sector[si]; ++si)
		{
			out_file << bt2MA_vert_per_sheet[vi][si] - r << " ";
		}
		out_file<<endl;
		// faces associated with each sector
		for (size_t si = 0; si < t_n_sector[vi]; ++si)
		{
			assoc_faces.clear();
			this->getAssocFaces(vi, si, assoc_faces);
			for (auto it = assoc_faces.begin(); it != assoc_faces.end(); ++it)
			{
				out_file << *it << " ";
			}
			out_file << endl;
		}
	}
	out_file.precision(ss);
	out_file.setf(flgs);
	out_file.close();
}

void SteinerGraph::getBurnedStEdgesFromBurnTree(
	vector<std::pair<TriEdge, int> >& _burntStE_assocF_list, 
	vector< std::pair<size_t,list<TriEdge>> >* _face_with_crossed_paths /*= nullptr*/)
{
	/* get all burnt st edges that are part of any burn path tree */
	typedef std::pair<TriEdge, int> edge_face_tuple;
	set<edge_face_tuple> burntE_assocF_set;
	auto make_edge_face_tuple = [] (unsigned _u, unsigned _v, int _f) 
	{
		return make_pair(TriEdge(_u, _v), _f);
		//return trimesh::ivec3(_u, _v, _f);
	};

	/* trace burn path trees. associate each burn path edge with face burnt by it. */
	vector<bool> visited(m_stSubdiv.sizeOfVts(), false);
	queue<int> out_branches;
	set<unsigned> pre_vts_set;
	map<int, vector<TopoGraph::EdgeIdx> > out_branch_pair_map;
	set<TopoGraph::EdgeIdx> v_tedges_visited_set;
	set<TopoGraph::EdgeIdx> u_tedges_visited_set;
	
	/*unsigned v_debug = 19;
	TriEdge e_debug = util::makeEdge(668, 373787);*/

	// collect burnt st edge and its assoc face into a list
	vector<TriEdge> stEdges_of_f;
	for (unsigned fi = 0; fi < m_origG->faces.size(); ++fi)
	{
		int from = -1, to = -1;
		stEdges_of_f.clear();
		m_stSubdiv.getStEdgesOfFace(fi, stEdges_of_f);
		for (size_t i = 0; i < stEdges_of_f.size(); ++i)
		{
			const auto& st_e = stEdges_of_f[i];
			if ( getBurnScheme() == STEINER_ONLY && 
				( !m_stSubdiv.isSteinerVert( st_e[ 0 ] ) || !m_stSubdiv.isSteinerVert( st_e[ 1 ] ) )
				)
				continue;
			if ( isBurnPath( st_e, fi, from, to ) )
			{
				burntE_assocF_set.insert(make_edge_face_tuple(from, to, fi));
			}
		}
	}
	cout << "initial step of collecting burnt edges -------------: " << endl;
	cout << "# burnt edges initially collected: " << burntE_assocF_set.size() << endl;

	/* fix any possible burn path intersections in a triangle face */
	bool reverse[] = {false, false, false};
	TriEdge edges[3];
	vector<unsigned> stVts_onE;//vts on cur edge
	vector<int> stVts_other; // vts on other edge(s)
	vector<int> vts2traverse; // traverse vts in this order later
	int nVts_eachEdge[3]; // # of st vts on each edge
	vector<int> vOrder_eID_map; // st vert -> edge id
	list<TriEdge> potential_crossing_paths; // paths within cur face that could cross each other
	auto crossing_paths_cpy = potential_crossing_paths;
	// crossings will be fixed later

	// update the burn path topology, where the path lies in face:
	// _new_to's prev <- _new_from, also other topo graph related info
	auto update_pathTopo_within_face = [&] (unsigned _to, unsigned _new_from, unsigned _old_from, unsigned _fi) 
	{
		auto& to_preTedges = this->prev_tedge[_to];
		auto& to_preVts = this->prev_vert[_to];

		// change _to's prev vert to _new_from
		TopoGraph::EdgeIdx to_te = (TopoGraph::EdgeIdx)this->mapTopo(_to, _fi);
		TopoGraph::EdgeIdx newFrom_te = (TopoGraph::EdgeIdx)this->mapTopo(_new_from, _fi);
		TopoGraph::EdgeIdx oldFrom_te = to_preTedges[to_te];
		for (TopoGraph::EdgeIdx te = 0; te < to_preVts.size(); ++te)
		{
			if (to_preVts[te] == _old_from)
			{
				to_preVts[te] = _new_from;
				to_preTedges[te] = newFrom_te;
			}
		}

		// debug
#ifdef _DEBUG_BURN
		unsigned to_debug = 8695, new_from_debug = _new_from, old_from_debug = _old_from;
		if (_to == to_debug && _new_from == new_from_debug && _old_from == old_from_debug)
		{
			cout << "updating path: to "<<to_debug<<
				" now from "<<new_from_debug<<" before from "<<old_from_debug<<endl;
			for (TopoGraph::EdgeIdx te = 0; te < this->prev_vert[_to].size(); ++te)
			{
				cout << "pre_vert, pre_tedge = "<<this->prev_vert[_to][te]<<
					", "<<(unsigned)(this->prev_tedge[_to][te])<<endl;
			}
			cout << endl;
		}
#endif
	};

	set<unsigned> f_debug;
	//f_debug.insert( 7276 );
	auto fixBurnPathIntersections = [&] (unsigned _fi) -> bool 
	{
		bool fixed = true;
		int intersect_cnt = 0;
		for (auto cur_path_it = potential_crossing_paths.begin(); 
			cur_path_it != potential_crossing_paths.end(); )
		{
			bool switched = false;
			// check if cur path intersect with other path
			for (auto other_path_it = potential_crossing_paths.begin(); 
				other_path_it != potential_crossing_paths.end(); )
			{
				// switch if intersected
				bool intersected = ((*cur_path_it)[0] - (*other_path_it)[0]) * 
					((*cur_path_it)[1] - (*other_path_it)[0]) * 
					((*cur_path_it)[0] - (*other_path_it)[1]) * 
					((*cur_path_it)[1] - (*other_path_it)[1]) < 0;
				if ( !intersected )
				{
					// move to the next "other" path
					++ other_path_it;
				}
				else
				{
					intersect_cnt ++;
					if (f_debug.count(_fi))
					{
						std::cout << "intersection detected!-----" << endl;
					}
					// actually switch: 
					// - find out direction of each path & switch them;
					// - append switched paths to the list;
					// - remove old paths;
					unsigned cur_from = (*cur_path_it)[0];
					unsigned cur_to = (*cur_path_it)[1];
					unsigned other_from = (*other_path_it)[0];
					unsigned other_to = (*other_path_it)[1];

					TriEdge temp;
					/*bool cur_nonManifold = isOnNonManifoldCurve(
						vts2traverse[cur_from], vts2traverse[cur_to], temp);*/
					int cur_assocF = _fi; //cur_nonManifold ? _fi : -1;
					/*bool other_nonManifold = isOnNonManifoldCurve(
						vts2traverse[other_from], vts2traverse[other_to], temp);*/
					int other_asssocF = _fi; //other_nonManifold ? _fi : -1;

					if (f_debug.count(_fi))
					{
						cout << "Done: assoc face set."<< endl;
					}

					TriEdge new_edges[2] = {
						TriEdge(other_from, cur_to), 
						TriEdge(cur_from, other_to)
					};
					bool cur_cross_pair_fixed = false;
					// only need to change one path to resolve intersection
					// we are about to try the first new edge. 
					// if it can't be used, then try the second edge.
					for (unsigned ei = 0; ei < 2; ++ei)
					{
						if (cur_cross_pair_fixed)
							break;

						TriEdge switched_e = new_edges[ei]; // <new_from, old_to>
						unsigned old_from = new_edges[(ei+1)%2][0];
						int eid_1 = vOrder_eID_map[switched_e[0]];
						int eid_2 = vOrder_eID_map[switched_e[1]];

						// about to add s_e. check if s_e overlaps with a tri edge.
						int overlap_ei = -1;
						if ( eid_1 >= 0 && (eid_1 == eid_2) ) // completely inside an orig. edge
						{
							overlap_ei = eid_1;
						}
						else if ( eid_1 < 0 && eid_2 < 0 ) // completely overlaps an orig. edge
						{
							if (eid_1 == -1 && eid_2 == -2 || eid_1 == -2 && eid_2 == -1)
								overlap_ei = 0;
							else if (eid_1 == -2 && eid_2 == -3 || eid_1 == -3 && eid_2 == -2)
								overlap_ei = 1;
							else
								overlap_ei = 2;
						}
						else if ( (eid_1 < 0 ) && (	// one end is a tri orig. v, the other in adjacent edge
							(eid_1 == -1 && (eid_2 == 0 || eid_2 == 2) ) || 
							(eid_1 == -2 && (eid_2 == 0 || eid_2 == 1) ) || 
							(eid_1 == -3 && (eid_2 == 1 || eid_2 == 2) ) ) )
						{
							overlap_ei = eid_2;
						}
						else if ( (eid_2 < 0 ) && (				
							(eid_2 == -1 && (eid_1 == 0 || eid_1 == 2) ) || 
							(eid_2 == -2 && (eid_1 == 0 || eid_1 == 1) ) || 
							(eid_2 == -3 && (eid_1 == 1 || eid_1 == 2) ) ) )
						{
							overlap_ei = eid_1;
						}
						TriEdge st_e = switched_e;
						if ( overlap_ei >= 0 ) // there is overlap !
						{
							// adjust s_e
							if (f_debug.count(_fi))
								cout << "overlap edges: " << overlap_ei << endl;

							// it's the immediate st edge on the orig. tri edge
							// that should be added
							if (overlap_ei == 2) // wrap the 0th vt to the last vert
							{
								if (switched_e[0] == 0)
									switched_e[0] += vts2traverse.size();
								else if (switched_e[1] == 0)
									switched_e[1] += vts2traverse.size();
							}
							if (switched_e[0] > switched_e[1])
							{
								st_e[0] = (switched_e[1] + 1) % vts2traverse.size();
								st_e[1] = switched_e[1];
							}
							else 
							{
								st_e[0] = switched_e[1] - 1;
								st_e[1] = switched_e[1] % vts2traverse.size();
							}
						} // handling overlap
						
						if (f_debug.count(_fi))
						{
							cout << st_e <<" about to be added. checking pre-existence/loop... "<< endl;
						}

						// before add the st edge, check pre-existence / loop
						TriEdge temp;
						int temp_from = -1, temp_to = -1;
						/*bool non_manifoldness = isOnNonManifoldCurve(
							vts2traverse[st_e[0]], vts2traverse[st_e[1]], temp);*/
						int assoc_f = _fi; /*non_manifoldness ? _fi : -1;*/
						auto e_to_add = util::makeEdge(vts2traverse[st_e[1]], vts2traverse[st_e[0]]);
						bool already_exist = isBurnPath(e_to_add, assoc_f, temp_from, temp_to);
						if (already_exist)
						{
							if (temp_from == st_e[0] && temp_to == st_e[1])
							{
								cout << "Fatal logic error: to-add switch path already exists! Quiting..." << endl;
								exit(-1);
							}
							else // the reverse path exists. Loop!
							{
								cur_cross_pair_fixed = false;
								continue;
							}
						}

						if (f_debug.count(_fi))
						{
							cout << st_e <<" safe to add."<< endl;
						}

						// now safely add switch path
						burntE_assocF_set.insert( 
							make_edge_face_tuple(
							vts2traverse[st_e[0]], vts2traverse[st_e[1]], assoc_f) );
						potential_crossing_paths.push_back(st_e);
						cur_cross_pair_fixed = true;
						switched = true;
						if (f_debug.count(_fi))
						{
							cout << st_e <<" added"<< endl;
						}
						// don't forget to update burn path topology
						update_pathTopo_within_face(
							vts2traverse[st_e[1]], vts2traverse[st_e[0]], vts2traverse[old_from], _fi);
						// remove old edges from burnt edge set and output list
						if (ei == 0)
						{
							if (f_debug.count(_fi))
							{
								cout << "deleting <"<<cur_from<<","<<cur_to<<">"<<endl;
							}
							burntE_assocF_set.erase( 
								burntE_assocF_set.find(
								make_edge_face_tuple(vts2traverse[cur_from], vts2traverse[cur_to], cur_assocF)
								) );
							cur_path_it = potential_crossing_paths.erase(cur_path_it);
						}
						else
						{
							if (f_debug.count(_fi))
							{
								cout << "deleting <"<<other_from<<","<<other_to<<">"<<endl;
							}
							burntE_assocF_set.erase( 
								burntE_assocF_set.find(
								make_edge_face_tuple( vts2traverse[other_from], vts2traverse[other_to], other_asssocF )
								) );
							if (cur_path_it == other_path_it)
							{
								cur_path_it = potential_crossing_paths.erase(other_path_it);
								other_path_it = cur_path_it;
							}
							else
								other_path_it = potential_crossing_paths.erase(other_path_it);
						}
					}
					if (!cur_cross_pair_fixed)
					{
						cout << "Logic error: cannot fix intersected paths!" << endl;
						fixed = false;
					}

				//	switched = true;
				//	break;
				} // intersected
				
			} // for each (other) path, test to see if intersected
			
			if (!switched)
			{
				++cur_path_it;
			}
			/*if ( !switched )
				++ cur_path_it;
			else
				cur_fixed = true;*/
		}
		return intersect_cnt > 0 && fixed;
	}; // fixBurnPathIntersections()

	unsigned cnt_fixed_faces = 0;
	pair<unsigned, unsigned> range_eachEdge[3];
	TriEdge e;
	for (unsigned fi = 0; fi < m_origG->faces.size(); ++fi)
	{
		if (fi % 20000 == 0 /*1*/) // debug
			cout << "# of faces where crossings checked: " << fi << endl;

		TriFace f = m_origG->faces[fi];

		//
		// TODO: what follows should be abstracted into a separate function.
		//

		// whether we need to reverse the steiner points order
		// for each edge
		for (unsigned ei = 0; ei < 3; ei ++)
		{
			TriEdge e = util::makeEdge(f[ei], f[(ei+1)%3]);
			edges[ei] = e;
			if (e[0] != f[ei])
				reverse[ei] = true;
		}

		// save the traversal order as we traverse each st&orig vert
		vts2traverse.clear();
		for (unsigned ei = 0; ei < 3; ei ++)
		{
			TriEdge e = edges[ei];
			stVts_onE.clear();
			m_stSubdiv.getStVertIndicesOnTriEdge(e, stVts_onE);

			nVts_eachEdge[ei] = stVts_onE.size();

			if (reverse[ei])
			{
				std::reverse(stVts_onE.begin(), stVts_onE.end());
				stVts_onE.insert(stVts_onE.begin(), e[1]);
			}
			else
			{
				stVts_onE.insert(stVts_onE.begin(), e[0]);
			}
			std::copy(stVts_onE.begin(), stVts_onE.end(), back_inserter(vts2traverse));

		}

		// build mapping: vert order -> the edge the vert's on
		// don't forget: hold-place for triangle's 3 vertices
		vOrder_eID_map.resize(vts2traverse.size());
		vector<int>::iterator fill_itor = vOrder_eID_map.begin();
		*fill_itor = -1; // first tri vert
		fill_itor += 1;
		std::fill_n(fill_itor, nVts_eachEdge[0], 0);
		fill_itor += nVts_eachEdge[0];
		*fill_itor = -2; // second tri vert
		fill_itor += 1;
		std::fill_n(fill_itor, nVts_eachEdge[1], 1);
		fill_itor += nVts_eachEdge[1];
		*fill_itor = -3; // third tri vert
		fill_itor += 1;
		std::fill_n(fill_itor, nVts_eachEdge[2], 2);

		// we collect (burnt) steiner edges lying within a face as potentially crossing paths
		// note: ends of edge use index into the vts2traverse list to reflect a certain order.
		unsigned start_v = 1;
		unsigned end_v = 1 + nVts_eachEdge[0];
		unsigned next_start_v, next_end_v;
		potential_crossing_paths.clear();
		// where do st vts on each edge start & end?
		range_eachEdge[0].first = start_v;
		range_eachEdge[0].second = end_v;
		range_eachEdge[1].first = range_eachEdge[0].second+1;
		range_eachEdge[1].second = range_eachEdge[0].second+1+nVts_eachEdge[1];
		range_eachEdge[2].first = (range_eachEdge[1].second+1) % vts2traverse.size();
		range_eachEdge[2].second = (range_eachEdge[1].second+1+nVts_eachEdge[2]) % vts2traverse.size();
		// now form st edges lying in face
		for (unsigned ei = 0; ei < 3; ei ++)
		{
			// update the start / end index into vts2traverse for next edge
			const auto& s_e_pr = range_eachEdge[ei];
			const auto& s_e_pr_next = range_eachEdge[(ei+1)%3];
			start_v = s_e_pr.first;
			end_v = s_e_pr.second;
			next_start_v = s_e_pr_next.first;
			next_end_v = s_e_pr_next.second;
			unsigned oppo_v_order = next_end_v;

			// connect vts on adjacent two edges: 
			// e1 to e2 and the oppo v of e1, 
			// e2 to e3 and the oppo v of e2,
			// e3 to e1 and the oppo v of e3
			int temp_from = -1, temp_to = -1;
			for (unsigned order = start_v; order != end_v; order = (order+1)%vts2traverse.size() )
			{
				for (unsigned other_order = next_start_v; 
					other_order != next_end_v; other_order = (other_order+1)%vts2traverse.size() )
				{
					e = TriEdge(vts2traverse[order], vts2traverse[other_order]);
					if ( getBurnScheme() == STEINER_ONLY &&
						( !m_stSubdiv.isSteinerVert( e[ 0 ] ) || !m_stSubdiv.isSteinerVert( e[ 1 ] ) )
						)
						continue;
					if ( isBurnPath(e, fi, temp_from, temp_to) )
					{
						if (temp_from == e[0])
							potential_crossing_paths.push_back(TriEdge(order, other_order));
						else
							potential_crossing_paths.push_back(TriEdge(other_order, order));
					}
				}
				if ( getBurnScheme() == ORIGINAL_AND_STEINER)
				{
					e = TriEdge( vts2traverse[ order ], vts2traverse[ oppo_v_order ] );
					if ( isBurnPath( e, fi, temp_from, temp_to ) )
					{
						//potential_crossing_paths.push_back(TriEdge(temp_from, temp_to));
						if ( temp_from == e[ 0 ] )
							potential_crossing_paths.push_back( TriEdge( order, oppo_v_order ) );
						else
							potential_crossing_paths.push_back( TriEdge( oppo_v_order, order ) );
					}
				}
			}
		}

		/*unsigned start_v = 1;
		unsigned end_v = 1 + nVts_eachEdge[0];
		potential_crossing_paths.clear();
		for (unsigned ei = 0; ei < 3; ei ++)
		{
		// vts on cur edge
		stVts_onE.resize(end_v - start_v);
		for (unsigned i = 0; i < stVts_onE.size(); ++i)
		stVts_onE[i] = i + start_v; // index into vts2traverse


		// vts of the other side
		stVts_other.clear();
		if (ei == 0) // vts on other edges 1, 2, including the opposite vert
		{
		stVts_other.resize(vts2traverse.size() - end_v - 1);
		for (unsigned i = 0; i < stVts_other.size(); ++i)
		stVts_other[i] = i + end_v + 1; // index into vts2traverse
		}
		else if (ei == 1) // vts on other edges: 2
		{
		stVts_other.resize(vts2traverse.size() - end_v/ * - 1 + 1* /);
		for (unsigned i = 0; i < stVts_other.size(); ++i)
		stVts_other[i] = i + end_v + 1;// index into vts2traverse
		stVts_other.back() = 0;// index into vts2traverse to opposite vert of edge 1
		}

		// update the start / end index into vts2traverse for next edge
		start_v = end_v;
		end_v = end_v + 1 + nVts_eachEdge[(ei + 1) % 3];

		for (unsigned i = 0; i < stVts_other.size(); ++i)
		{
		for (unsigned j = 0; j < stVts_onE.size(); ++j)
		{
		int order_other = stVts_other[i];
		int order = stVts_onE[j];
		int v = vts2traverse[order];
		int v_other = vts2traverse[order_other];
		TriEdge e = TriEdge(vts2traverse[order], vts2traverse[order_other]);
		bool reversed = false;
		if (burntE_assocF_set.count(make_edge_face_tuple(e[0], e[1], -1)) > 0 )
		{
		potential_crossing_paths.push_back(TriEdge(order, order_other));
		}
		else if (burntE_assocF_set.count(make_edge_face_tuple(e[1], e[0], -1)) > 0)
		{
		potential_crossing_paths.push_back(TriEdge(order_other, order));
		}
		}
		}
		}*/

		if (f_debug.count(fi))
		{
			cout << "----fixing face " << fi <<"----"<< endl;
			cout << "before switching... " << endl;
			cout << f <<": "
				<<m_stSubdiv.getVert(f[0])<<","
				<<m_stSubdiv.getVert(f[1])<<","
				<<m_stSubdiv.getVert(f[2])<< endl;

			cout << "traverse order: ";
			std::for_each( vts2traverse.begin(), vts2traverse.end(), 
				[] (unsigned _vi) { cout << _vi << " "; } );
			cout << endl;

			cout << "in face graph edges: ";
			std::for_each( potential_crossing_paths.begin(), potential_crossing_paths.end(),
				[] (const TriEdge& _e) { cout << _e<<" "; });
			cout << endl;
		}

		// debug begin
		crossing_paths_cpy = potential_crossing_paths;
		// debug end

		// now potential_crossing_paths is ready. check and fix any crossings.
		if (fixBurnPathIntersections(fi))
		{
			//cout << "fixed: "<<fi<<endl;
			cnt_fixed_faces ++;

			// debug begin
			if (_face_with_crossed_paths)
			{
				for (auto it = crossing_paths_cpy.begin(); it != crossing_paths_cpy.end();
					it ++)
				{
					(*it)[0] = vts2traverse[(*it)[0]];
					(*it)[1] = vts2traverse[(*it)[1]];
				}
				_face_with_crossed_paths->push_back(make_pair(fi, crossing_paths_cpy));
			}
			// debug end
		}

		// debug begin
		if (f_debug.count(fi))
		{
			cout << "after switching... " << endl;

			cout << "in face graph edges: ";
			std::for_each( potential_crossing_paths.begin(), potential_crossing_paths.end(),
				[] (const TriEdge& _e) { cout << _e<<" "; });
			cout << endl;
		}
		// debug end

	} // fixed all problematic faces

	cout << "# faces whose burn paths crossed&fixed: " << cnt_fixed_faces <<endl;

	_burntStE_assocF_list.assign(burntE_assocF_set.begin(), burntE_assocF_set.end());
	burntE_assocF_set.clear();

	// debug begin
	/*for (auto it  = _burntStE_assocF_list.begin(); it != _burntStE_assocF_list.end(); ++it)
	{
		if ( util::makeEdge(it->first[0],it->first[1]) == e_debug )
			cout << it->first << " associated with face "<<it->second<<endl;
	}*/
	// debug end

	/* output burntEdge association to intermediate file for later inspection */
	//const char* file_name = "burnPath_face_assoc.txt";
	//ofstream out_file(file_name);
	//if (!out_file.good())
	//	cout << "Error: couldn't open "<< file_name<< endl;
	//else
	//{
	//	for (auto it = _burntStE_assocF_list.begin(); it != _burntStE_assocF_list.end(); ++it)
	//	{
	//		out_file << it->first<<"->"<<it->second<<endl;
	//	}
	//	out_file.close();
	//}

	///* a few checks */
	//cout << endl << "test: does each vert has valid number of prev branches? " << endl;
	//set<int> bad_vts;
	//typedef pair<unsigned, TopoGraph::EdgeIdx> pre_v_te_pair;
	//set<pre_v_te_pair> pre_v_te_pair_set;
	//set<unsigned> vts_debug;
	//vts_debug.insert(15699);
	//vts_debug.insert(15700);
	//for (unsigned vi = 0; vi < m_stSubdiv.sizeOfVts(); ++vi)
	//{
	//	const auto& tg = t_graphs[vi];
	//	pre_v_te_pair_set.clear();
	//	for (TopoGraph::EdgeIdx tei = 0; tei < tg.t_edges.size(); ++tei)
	//	{
	//		if (!is_leading_topo_edge(vi, tei))
	//			continue;
	//		pre_v_te_pair_set.insert( std::make_pair(prev_vert[vi][tei], prev_tedge[vi][tei]) );
	//	}

	//	if ( !( (tg.type == BOUNDARY && pre_v_te_pair_set.size() == 1 /*&& *pre_vts_set.begin() == vi */)
	//		|| (tg.type == MANIFOLD && pre_v_te_pair_set.size() == 1)
	//		|| (tg.type == NONMANIFOLD && pre_v_te_pair_set.size() <= t_graphs[vi].t_edges.size() - 1)
	//		|| (tg.type == JUNCTION)
	//		) )
	//	{
	//		bad_vts.insert(vi);
	//	}

	//	/* print out going branches of bad vts */
	//	if (vts_debug.count(vi) || bad_vts.count(vi))
	//	{
	//		cout << "prev <vert, topo edge> of " << vi <<" "<<topo_type_str[tg.type]<<" : ";
	//		for (TopoGraph::EdgeIdx tei = 0; tei < prev_vert[vi].size(); ++tei)
	//		{
	//			int pre = getPrevVert(vi, tei);
	//			unsigned te = getPrevTopoEdge(vi, tei);
	//			if (!is_leading_topo_edge(vi, tei))
	//				cout << "tei "<<(int)tei << " depends on "<< prev_tedge[vi][tei]<<". ";
	//			cout <<"<"<<pre<<" "<<topo_type_str[t_graphs[pre].type]<<", "
	//				<<(unsigned)te<<"> "; 
	//		}
	//		cout << endl;
	//	}
	//}
	//cout << "# bad vts: " << bad_vts.size() << endl;
	//unsigned bad_vts_manifold = 0;
	//unsigned bad_vts_nonmanifold = 0;
	//unsigned bad_vts_bndry = 0;
	//for (auto v_it = bad_vts.begin(); v_it != bad_vts.end(); ++v_it)
	//{
	//	auto type = t_graphs[*v_it].type;
	//	//cout << *v_it<<", "<<topo_type_str[t_graphs[*v_it].type] << endl;
	//	if (type == MANIFOLD)
	//		bad_vts_manifold ++;
	//	else if (type == NONMANIFOLD)
	//		bad_vts_nonmanifold ++;
	//	else if (type == BOUNDARY)
	//		bad_vts_bndry ++;
	//}
	//cout << "# bad manifold vts: " << bad_vts_manifold << endl;
	//cout << "# bad non-manifold vts: " << bad_vts_nonmanifold << endl;
	//cout << "# bad boundary vts: " << bad_vts_bndry << endl;
	//cout << endl;
} // end of SteinerGraph::getBurnedStEdgesFromBurnTree()

void SteinerGraph::output_in_mathematica_format(
	const vector< std::pair<size_t,list<TriEdge>> >& _face_with_crossed_paths
	)
{
	/* 0. create output file */
	ofstream outfile("faces_w_crossed_paths.txt");
	int elem_cnt = 0;
	if ( !outfile )
	{
		outfile << "Error: testing failed!" << endl;
		exit(-1);
	}

	// biggest {
	outfile << '{' << endl;

	/*1. output a list of all vertices of stg*/
	// output list of centroids
	outfile << "(* BEGIN: vertices list *)"<<endl;
	outfile << "{ ";
	const auto& all_vts = m_stSubdiv.getAllVts();
	for (int i = 0; i < all_vts.size(); ++i)
	{
		auto& c = all_vts[i];
		outfile << "{"<< c[0]<<","<<c[1]<<","<<c[2]<<"}";
		if (i + 1 < all_vts.size())
		{
			outfile << ", ";
		}
	}
	outfile << " }, " << endl;
	outfile << "(* END: vertices list *)"<<endl;

	/* 2. output each face with its associated crossed paths */
	auto last_pair_it = _face_with_crossed_paths.end(); last_pair_it--;
	outfile << "(* BEGIN: face-crossing pair list *)"<<endl;
	outfile << "{ ";
	for ( auto it = _face_with_crossed_paths.begin(); it != _face_with_crossed_paths.end(); ++it )
	{
		// begin current face-crossing pair
		outfile<<"{";

		// output face
		auto& f = m_origG->faces[it->first];
		outfile << "{"<<f[0]+1<<","<<f[1]+1<<","<<f[2]+1<<"},";
		// output crossed paths
		const auto& paths = it->second;
		auto last_p_it = paths.end(); last_p_it--;
		outfile<<"{";
		for ( auto p_it = paths.begin(); p_it != paths.end(); ++p_it )
		{
			outfile<<"{"<<(*p_it)[0]+1<<","<<(*p_it)[1]+1<<"}";
			if ( p_it != last_p_it )
				outfile << ", ";
		}
		outfile<<"}";

		// end of current pair
		outfile<<"}";
		if (it != last_pair_it)
		{
			outfile << ", ";
		}
	}
	outfile << " }" << endl;
	outfile << "(* END: face-crossing pair list *)"<<endl;

	// biggest }
	outfile << '}' << endl;
	outfile.close();

	cout << "faces w. crossed paths: output done." << endl;
}

/// -------------private helpers-----------------

void SteinerGraph::init()
{
	// init burn stats to invalid values
	resetVtsBurnStats();
	max_diffDistOrigFaces = numeric_limits<float>::min();
	min_diffDistOrigFaces = numeric_limits<float>::max();

	// initially all faces are active (e.g. not pruned)
	m_faceActiveFlag.assign(m_origG->faces.size(), true);
}

void SteinerGraph::resetVtsBurnStats()
{
	// burn stats to invalid values
	min_burnDistOrigVts = numeric_limits<float>::max();
	max_burnDistOrigVts = -1.0f;
	min_burnDistAllVts = numeric_limits<float>::max();
	max_burnDistAllVts = -1.0f;
}

void SteinerGraph::mergeAssocFaces( 
	vector<int>& _assoc_faces, 
	const vector<int>& _idxes, 
	const TopoGraph::TEdgeList& _tedges )
{
	for (vector<int>::const_iterator itor = _idxes.begin(); 
		itor != _idxes.end(); 
		++itor)
	{
		const auto& faces_assoc_temp = _tedges[*itor].second;
		_assoc_faces.insert(_assoc_faces.end(), 
			faces_assoc_temp.begin(), faces_assoc_temp.end());
	}
}

void SteinerGraph::openTopo( 
	const TopoGraph& _tg, 
	const TopoGraph::EdgeIdx _te, 
	vector<bool>& _final, 
	vector< std::pair<TopoGraph::EdgeIdx,int> >& _fnIdxes )
{
	// if _te is a boundary sheet, chances are
	// there might be boundary sheets still not finalized. finalize all of them.
	// For each vertex, this is done at most once during the entire algorithm 
	queue< std::pair<TopoGraph::EdgeIdx, int> > oq;
	if ( _tg.is_topo_edge_open(_te) ) 
	{
		for (TopoGraph::EdgeIdx te_id = 0; te_id < _tg.t_edges.size(); ++te_id)
		{
			if ( _tg.is_topo_edge_open(te_id) && !_final[te_id] )
			{
				// add to queue cuz there might be more edges opened up 
				// b.c. this edge is finalized
				oq.push(make_pair(te_id, -1));
			}
		}
	}
	else
	{
		oq.push( make_pair(_te, -1) );
	}

	auto t_nbVts = _tg.t_nbVtsMap;
	auto t_edgeIdxMap = _tg.t_edgeIdxesMap;
	set<int> nbVts_set;
	vector<TopoGraph::EdgeIdx> nbNonFinal;

	while(!oq.empty())
	{
		auto fn_pair = oq.front();
		TopoGraph::EdgeIdx tei = fn_pair.first;
		int exposed_by = fn_pair.second;
		TriEdge te = _tg.t_edges[tei].first;
		oq.pop();
		if ( _final[tei] )
			continue;

		// tei will be finalized
		_fnIdxes.push_back( make_pair(tei, exposed_by) );
		_final[tei] = true;

		if (!util::notPureLoop(te)) // a pure loop edge will not open any other edge
			continue;

		for (int i = 0; i < 2; ++i)
		{
			int v = te[i];
			auto nbs = t_nbVts.find(v)->second;
			if (nbs.size() == 1) // no more nbs for this end of the edge
				continue;
			// at least one non-final neighbor after finalization. could be open.
			nbVts_set.clear();// delete duplicate
			for (unsigned ni = 0; ni < nbs.size(); ++ni)
			{
				nbVts_set.insert(nbs[ni]);
			}
			// collect all non-final neighbor edges indices
			nbNonFinal.clear();
			for (auto itor = nbVts_set.begin(); itor != nbVts_set.end(); ++itor)
			{
				auto nbEdgeIdxes = (t_edgeIdxMap.find(util::makeEdge(v, *itor)))->second;
				for (auto eid = nbEdgeIdxes.begin(); eid != nbEdgeIdxes.end(); ++eid)
				{
					if (_final[*eid] == false)
						nbNonFinal.push_back(*eid);
				}
			}
			// the edge is open if it has only one non-final neighbor u, 
			// and this edge is not a boundary/open sheet (the unchanged nbVtsMap[u] != 1)
			// add this open edge to the queue
			if (nbNonFinal.size() == 1 && !_tg.is_topo_edge_open(nbNonFinal[0]) &&
				util::notLoop( _tg.t_edges[nbNonFinal[0]].first ))
			{
				oq.push( make_pair(nbNonFinal[0], exposed_by == -1 ? tei : exposed_by) );
			}
		}
	} // end of while queue not empty ...
} // SteinerGraph::openTopo()

void SteinerGraph::assignBoundaryStartDist(
	const vector<float>& _dist2surf)
{
	set<unsigned> debug_vts;

	for (unsigned vi = 0; vi < m_stSubdiv.sizeOfVts(); ++vi)
	{
        //if (/*vi % 10000 == 0*/ 1)
        if(vi == 233)
        {
            cout << "vi=" << vi << endl; // debug 233
        }

		auto n_sector = t_n_sector[vi];
		auto tt = getTopoType(vi);
		bool nbhood_has_2d_bits = TopoGraph::has2dNbhood( vi );
		final[ vi ].assign( n_sector, false );
		min_tedge[vi] = 0;
		burnt[vi] = nbhood_has_2d_bits ? false : true;
		bt2MA_vert_per_sheet[ vi ].assign( n_sector, infiniteBurnDist() );
		bt2MA_vert[ vi ] = nbhood_has_2d_bits ? infiniteBurnDist() : _dist2surf[ vi ];
		prev_vert[vi].assign(n_sector, -1);
		prev_tedge[vi].assign(n_sector, 255);
		if ( !nbhood_has_2d_bits )
			continue;

		//  for each open topo edge
		//  treat it like a boundary sheet
		TopoGraph::EdgeIdx tei;
		for ( tei = 0; tei < n_sector; ++tei )
		{
			if ( TopoGraph::isOpenSheet( vi, tei ) )
			{
				bt2MA_vert_per_sheet[ vi ][ tei ] = _dist2surf.empty() ? 0.0f : _dist2surf[ vi ];
				bt2MA_vert[ vi ] = _dist2surf.empty() ? 0.0f : _dist2surf[ vi ];
				prev_vert[ vi ][ tei ] = vi;
				prev_tedge[ vi ][ tei ] = tei;
				min_tedge[ vi ] = tei;
			}
		}
	}
} // SteinerGraph::assignBoundaryStartDist()

float SteinerGraph::update_with_protection( int _vi, float _bt ) const
{
	return max( _bt, bt3MA_vert[ _vi ] );
}

void SteinerGraph::recover_bt_for_orig_vts()
{
	vector<int> nb_vts_on_s;
	vector<int> assoc_faces;
	for ( auto vi = 0; vi < m_origG->vts.size(); ++vi )
	{
		// skip if this is a boundary vertex
		auto tt = getTopoType( vi );
		if ( !TopoGraph::has2dNbhood( vi ) || tt == BOUNDARY_2D )
		{
			burnt[ vi ] = true;
			continue;
		}
		// initialize burn time value
		bt2MA_vert[ vi ] = 0.0f;
		auto cur_v = m_origG->vts[ vi ];
		// determine bt2 on each sheet,
		// the final burn time of the vert is the max among all sheets
		for ( TopoGraph::EdgeIdx si = 0; si < t_n_sector[ vi ]; ++si )
		{
			// skip if this is a boundary sheet
			if ( TopoGraph::isOpenSheet(vi, si) )
				continue;
			float bt_on_s = infiniteBurnDist();
			int pre_vert, pre_s;
			assoc_faces.clear();
			getAssocFaces( vi, si, assoc_faces );
			for ( auto i = 0; i < assoc_faces.size(); ++i )
			{
				auto fi_on_s = assoc_faces[ i ];
				nb_vts_on_s.clear();
				m_stSubdiv.getNbVtsInFace( vi, fi_on_s, true, nb_vts_on_s );

				// compute candidate bt2 from each nb vert
				for ( auto i = 0; i < nb_vts_on_s.size(); ++i )
				{
					auto nb_v = nb_vts_on_s[ i ];
					auto nb_s = mapTopo( nb_v, fi_on_s );
					auto cand_bt2 = update_with_protection(
						vi,
						bt2MA_vert[ nb_v ] + trimesh::dist( cur_v, m_stSubdiv.getVert( nb_v ) )
					);
					cand_bt2 < bt_on_s ? ( bt_on_s = cand_bt2, pre_vert = nb_v, pre_s = nb_s ) : 0;
				}
			}
			
			if ( bt_on_s == infiniteBurnDist() )
			{
				bt2MA_vert[ vi ] = infiniteBurnDist();
				prev_vert[ vi ][ si ] = -1;
				prev_tedge[ vi ][ si ] = 0;
			}
			else
			{
				bt2MA_vert_per_sheet[ vi ][ si ] = bt_on_s;
				bt2MA_vert[ vi ] = std::max( bt2MA_vert[ vi ], bt_on_s );
				prev_vert[ vi ][ si ] = pre_vert;
				prev_tedge[ vi ][ si ] = pre_s;
			}
		}
		burnt[ vi ] = bt2MA_vert[ vi ] < infiniteBurnDist();
	}
}

void SteinerGraph::findBoundaryBurnDir(
	vector<trimesh::vec>& _bndry_dirs,
	map<unsigned, unsigned>& _vIdx_burnDirIdx_map)
{
	_bndry_dirs.clear();
	_bndry_dirs.reserve(t_cnt[BOUNDARY_2D]);

	// 1. set up boundary steiner vertices' direction.
	vector<unsigned> stVts_on_e;
	for (auto e_it = this->m_origG->edges.begin(); e_it != this->m_origG->edges.end(); ++e_it)
	{
		stVts_on_e.clear();
		m_stSubdiv.getStVertIndicesOnTriEdge(e_it-m_origG->edges.begin(), stVts_on_e);

		if (stVts_on_e.empty())
			continue;

		unsigned st_v = stVts_on_e[0];

		if (getTopoType(st_v) != BOUNDARY_2D)
			continue;

		// else, the steiner vts on this edge are all BOUNDARY vts
		// compute the burn dir as the flipped normal of the edge: face_normal X edge_dir
		int fi = this->m_origG->getNbFaces( e_it - m_origG->edges.begin() )[ 0 ];
		const auto& f = m_origG->faces[fi];
		const auto& p1 = m_stSubdiv.getVert(f[0]);
		const auto& p2 = m_stSubdiv.getVert(f[1]);
		const auto& p3 = m_stSubdiv.getVert(f[2]);
		const auto& face_nml = (p1-p2).cross(p1-p3);

		const auto& e_u = m_stSubdiv.getVert((*e_it)[0]);
		const auto& e_v = m_stSubdiv.getVert((*e_it)[1]);
		const auto& edge_dir = e_u-e_v;

		auto burn_dir = trimesh::normalize(face_nml.cross(edge_dir));
		unsigned oppo_vi = util::oppositeVert(f, *e_it);
		const auto& oppo_v = m_stSubdiv.getVert(oppo_vi);

		// we want the burn dir to point inward the tri face
		if ( trimesh::normalize(trimesh::vec(oppo_v-m_stSubdiv.getVert(st_v))).dot(burn_dir) < 0.0f )
			burn_dir = -burn_dir;

		// assign this dir to all steiner vts on this edge 
		_bndry_dirs.push_back(burn_dir);
		for (auto st_it = stVts_on_e.begin(); st_it != stVts_on_e.end(); ++st_it)
			_vIdx_burnDirIdx_map[*st_it] = _bndry_dirs.size() - 1;
	}

	// 2. set up boundary orig. tri vts using the neighbor edges
	for (unsigned vi = 0; vi < m_origG->vts.size(); ++vi)
	{
		if ( !TopoGraph::hasAnyOpenSheet( vi ) )
			continue;

		// find nb st vts, cuz we are going to use the fact that they are all boundary
		auto nb_vts = m_origG->nbVtsOfVert[(int)vi];
		trimesh::vec cur_dir(0.0f, 0.0f, 0.0f);
		for (auto nbV_it = nb_vts.begin(); nbV_it != nb_vts.end(); ++nbV_it)
		{
			if ( !TopoGraph::hasAnyOpenSheet( *nbV_it ) )
				continue;

			// then, there is a unique neighbor face.
			// use the normal of that face to compute dir.
			TriEdge e = util::makeEdge(vi, *nbV_it);
			auto nb_fs = this->m_origG->getNbFaces( e );
			if ( nb_fs.empty() )
				continue;
			int fi = nb_fs[ 0 ];
			const auto& f = m_origG->faces[fi];
			const auto& p1 = m_stSubdiv.getVert(f[0]);
			const auto& p2 = m_stSubdiv.getVert(f[1]);
			const auto& p3 = m_stSubdiv.getVert(f[2]);
			const auto& face_nml = (p1-p2).cross(p1-p3);

			const auto& e_u = m_stSubdiv.getVert(e[0]);
			const auto& e_v = m_stSubdiv.getVert(e[1]);
			const auto& edge_dir = e_u-e_v;

			auto burn_dir = trimesh::normalize(face_nml.cross(edge_dir));
			unsigned oppo_vi = util::oppositeVert(f, e);
			const auto& oppo_v = m_stSubdiv.getVert(oppo_vi);

			// we want the burn dir to point inward the tri face
			if ( 
				trimesh::normalize(trimesh::vec(oppo_v-m_stSubdiv.getVert(vi))).dot(burn_dir) 
				< 0.0f )
				burn_dir = -burn_dir;

			cur_dir += burn_dir;
		}
		cur_dir /= nb_vts.size();
		_bndry_dirs.push_back(cur_dir);
		_vIdx_burnDirIdx_map[vi] = _bndry_dirs.size() - 1;
	}
} // SteinerGraph::assignBoundaryBurnDir()

void SteinerGraph::dualize_poly_subdivision_simple(
	unsigned _fi, 
	const std::vector<vector<unsigned>>& _poly_subdivision,
	const map<unsigned, std::pair<trimesh::vec, float>>& _vert_burnDirDistPair_map,
	const map<TriEdge, vector<int> >& _burntEdge_assocFaces_map, 
	map<TriEdge, movingSum_edgeDual_struct>& edge_dual_map,
	DualOpt _dual_method,
	bool _verbose)
{
	vector<std::pair<unsigned, float>> dualV_dist_pair_curPolyFace;
	for (unsigned i = 0; i < _poly_subdivision.size(); ++i)
	{
		dualV_dist_pair_curPolyFace.clear();
		const auto& cur_poly = _poly_subdivision[i];

		for (unsigned ii = 0; ii < cur_poly.size(); ++ii)
		{
			TriEdge e = util::makeEdge(cur_poly[ii], cur_poly[(ii+1)%cur_poly.size()]);

			bool skip_this_edge = is_edge_dualizable(e[0], e[1], _fi, _burntEdge_assocFaces_map, _verbose);
			if (skip_this_edge)
				continue;
			if (_verbose)
				cout<<"dual vert "<<dualV_dist_pair_curPolyFace.back().first<<endl; // debug

			TriPoint dual_v;
			unsigned dualV_idx;
			const auto& dir_dist1 = _vert_burnDirDistPair_map.find(e[0])->second;
			const auto& dir_dist2 = _vert_burnDirDistPair_map.find(e[1])->second;
			float dual_bt2;
			if (_dual_method == SIMPLE_DUAL)
			{
				dual_v = trimesh::mix(m_stSubdiv.getVert(e[0]), m_stSubdiv.getVert(e[1]), 0.5f);
				dual_bt2 = (dir_dist1.second + dir_dist2.second) / 2.0f;
			}
			else // supposedly POINT_SOURCE_DUAL
			{
				float l1 = trimesh::dist(m_stSubdiv.getVert(e[0]), m_stSubdiv.getVert(e[1]));
				float l2 = dir_dist2.second - dir_dist1.second;

				dual_v = trimesh::mix(m_stSubdiv.getVert(e[0]), m_stSubdiv.getVert(e[1]), (l1+l2)*0.5f/l1);
				dual_bt2 = (l1+l2) * 0.5f + dir_dist1.second;
			}

			auto it = edge_dual_map.find(e);
			if (it == edge_dual_map.end())
			{
				dual_vts.push_back( dual_v );
				dualV_idx = dual_vts.size() - 1;

				bt2_MC.push_back( dual_bt2 );

				auto& tpl = edge_dual_map[e];
				tpl.dual_v_id = dualV_idx;
				tpl.bt2 = dual_bt2;
				tpl.sum = dual_v;
				tpl.count = 1;
			}
			else
			{
				dualV_idx = it->second.dual_v_id;

				auto& tpl = edge_dual_map[e];
				tpl.bt2 += dual_bt2;
				tpl.sum += dual_vts[dualV_idx];
				tpl.count += 1;
			}
			dualV_dist_pair_curPolyFace.push_back(std::make_pair(dualV_idx, dual_bt2));
		}

		TriPoint center_pos(0.0f, 0.0f, 0.0f);
		float center_dist = 0.0f;
		std::for_each( dualV_dist_pair_curPolyFace.begin(), dualV_dist_pair_curPolyFace.end(), 
			[&center_pos, &center_dist, this] (std::pair<unsigned, float> _pr) {
				center_pos = center_pos + this->dual_vts[_pr.first];
				center_dist = center_dist + _pr.second;
		} );
		center_pos = vec(center_pos) / (float)dualV_dist_pair_curPolyFace.size();
		center_dist /= (float)dualV_dist_pair_curPolyFace.size();

		dual_vts.push_back(center_pos);
		this->bt2_MC.push_back(center_dist);
		unsigned center_id = dual_vts.size()-1;
		std::for_each( dualV_dist_pair_curPolyFace.begin(), dualV_dist_pair_curPolyFace.end(), 
			[_fi, center_id, this] (std::pair<unsigned, float> _pr) {
				this->dual_edges.push_back( util::makeEdge(center_id, _pr.first) );
		} );
	}
} // SteinerGraph::dualize_poly_subdivision_simple()

void SteinerGraph::dualize_poly_subdivision(
	const std::vector<unsigned>& _vts, unsigned _fi, bool _f_part_of_pocket,
	const std::vector<vector<unsigned>>& _poly_subdivision,
	const map<unsigned, std::pair<trimesh::vec, float>>& _vert_burnDirDistPair_map,
	const map<TriEdge, vector<int> >& _burntEdge_assocFaces_map,
	map<TriEdge, movingSum_edgeDual_struct>& _edge_dual_map,
	DualOpt _edge_dual_opt/* = WEIGHT_CENTER_CLOSENESS_DUAL*/,
	DualOpt _poly_dual_opt/* = WEIGHT_CENTER_CLOSENESS_DUAL*/,
	float _w, float _eps,
	bool _verbose,
	ostream* _os )
{
	if(_verbose) 
		std::cout << "face " << _fi << std::endl;
	vector<std::pair<unsigned, float>> poly_duals;
	vector<unsigned> edge_dual_indices;
	vector<float> edge_dual_bt2;
	vector<float> edge_dual_bt3;
	vector<bool> edge_no_dual;
	auto f = this->m_origG->faces[_fi];
	TriEdge e;

	/* used by PSEUDO_INV_DUAL, WEIGHT_CENTER_CLOSENESS_DUAL */
	TriPoint p1, p2, p3;
	p1 = m_stSubdiv.getVert(f[0]);
	p2 = m_stSubdiv.getVert(f[1]);
	p3 = m_stSubdiv.getVert(f[2]);
	vector<float> weights(3, 0.0f);
	// stuff used to prepare for poly face dual
	Eigen::Matrix<double, 3, 2> P; // constant across all poly faces in the same subdivision
	P << p1[0]-p3[0], p2[0]-p3[0],
		p1[1]-p3[1], p2[1]-p3[1],
		p1[2]-p3[2], p2[2]-p3[2];
	Eigen::Vector3d p3_cpy(p3[0], p3[1], p3[2]);
	Eigen::Vector3d qi;
	Eigen::Matrix3d	A = Eigen::Matrix3d::Zero();
	Eigen::Vector3d Ai;
	double Ci;
	Eigen::Matrix<double, 3, 1> B = Eigen::Matrix<double, 3, 1>::Zero();
	Eigen::Vector3d dir;
	// intermediate result and poly face dual
	Eigen::Vector3d sol; // {s, t, dist of dual}
	TriPoint poly_dual_pos;
	float poly_dual_bt2;
	float poly_dual_bt3;
	// stuff used to prepare for poly edge dual
	TriPoint p1_, p2_;
	Eigen::Matrix<double, 3, 1> P_;
	Eigen::Vector3d p2_cpy;
	Eigen::Matrix2d A_;
	Eigen::Vector2d Ai_;
	double Ci_;
	Eigen::Matrix<double, 2, 1> B_;
	// intermediate result and poly edge dual
	Eigen::Vector2d sol_; // {s, t, dist of dual}

	/* used by WEIGHT_DUAL */
	vector<float> angle_weights;

	for (auto poly_it = _poly_subdivision.begin(); poly_it != _poly_subdivision.end(); ++poly_it)
	{
		unsigned qis[2]; // indices of two ends of the edge
		edge_dual_indices.clear();
		edge_dual_indices.resize(poly_it->size(), numeric_limits<unsigned>::max());
		edge_dual_bt2.clear();
		edge_dual_bt2.resize( poly_it->size(), numeric_limits<float>::max() );
		edge_dual_bt3.clear();
		edge_dual_bt3.resize( poly_it->size(), numeric_limits<unsigned>::max() );
		edge_no_dual.clear();
		edge_no_dual.resize( poly_it->size(), true );

		for ( unsigned poly_ei = 0; poly_ei < poly_it->size(); ++poly_ei )
		{
			qis[ 0 ] = *( poly_it->begin() + poly_ei );
			qis[ 1 ] = *( poly_it->begin() + ( poly_ei + 1 ) % poly_it->size() );
			TriEdge e = util::makeEdge( qis[ 0 ], qis[ 1 ] );

			// skip if cur edge is not dualizable.
			if ( !is_edge_dualizable( qis[ 0 ], qis[ 1 ], _fi, _burntEdge_assocFaces_map, _verbose ) )
				continue;
			edge_no_dual[ poly_ei ] = false;

			TriPoint edge_dual_pos;
			unsigned dualV_idx;
			// topo edge for the two ends of a poly edge
			int te1, te2;
			float dual_bt2;
			float dual_bt3;
			pair<vec, float> dir_dist1, dir_dist2;
			if ( !_f_part_of_pocket )
			{
				dir_dist1 = _vert_burnDirDistPair_map.find( qis[ 0 ] )->second;
				dir_dist2 = _vert_burnDirDistPair_map.find( qis[ 1 ] )->second;
			}

			/* 1. dualize edges */
			if ( _f_part_of_pocket ) // handle pocket edge specially
			{
				edge_dual_pos = trimesh::mix( m_stSubdiv.getVert( e[ 0 ] ), m_stSubdiv.getVert( e[ 1 ] ), 0.5f );
				dual_bt2 = infiniteBurnDist();
				dual_bt3 = ( bt3MA_vert[ e[ 0 ] ] + bt3MA_vert[ e[ 1 ] ] ) * 0.5f;
			}
			else if ( _edge_dual_opt == SIMPLE_DUAL )
			{
				edge_dual_pos = trimesh::mix( m_stSubdiv.getVert( e[ 0 ] ), m_stSubdiv.getVert( e[ 1 ] ), 0.5f );
				dual_bt2 = ( dir_dist1.second + dir_dist2.second ) * 0.5f;
				dual_bt3 = ( bt3MA_vert[ e[ 0 ] ] + bt3MA_vert[ e[ 1 ] ] ) * 0.5f;
			}
			else if ( _edge_dual_opt == PSEUDO_INV_DUAL || _edge_dual_opt == WEIGHT_CENTER_CLOSENESS_DUAL )
			{
				p1_ = m_stSubdiv.getVert( qis[ 0 ] );
				p2_ = m_stSubdiv.getVert( qis[ 1 ] );
				p2_cpy << p2_[ 0 ], p2_[ 1 ], p2_[ 2 ];
				p2_ = p1_ - p2_;
				P_ << p2_[ 0 ], p2_[ 1 ], p2_[ 2 ];
				A_ = Eigen::Matrix2d::Zero();
				B_ = Eigen::Matrix<double, 2, 1>::Zero();

				for ( unsigned j = 0; j < 2; ++j )
				{
					const auto& dir_dist_pair = _vert_burnDirDistPair_map.find( qis[ j ] )->second;
					const auto& burn_dir = dir_dist_pair.first;
					double di = dir_dist_pair.second;
					const auto& qi_ = m_stSubdiv.getVert( qis[ j ] );
					qi << qi_[ 0 ], qi_[ 1 ], qi_[ 2 ];
					dir << burn_dir[ 0 ], burn_dir[ 1 ], burn_dir[ 2 ];
					Ai_ << dir.dot( P_ ), -1.0;
					A_ += Ai_ * Ai_.transpose();
					Ci_ = ( p2_cpy - qi ).dot( dir ) + di;
					B_ += -1 * Ci_ * Ai_;
				}

				// the center of the cur edge
				TriPoint center( 0.0f, 0.0f, 0.0f );
				center = ( m_stSubdiv.getVert( qis[ 0 ] ) + m_stSubdiv.getVert( qis[ 1 ] ) ) / 2.0f;

				if ( _edge_dual_opt == PSEUDO_INV_DUAL )
				{
					float center_d = 0.0f;
					for ( unsigned j = 0; j < 2; ++j )
					{
						const auto& dir_dist_pair = _vert_burnDirDistPair_map.find( qis[ j ] )->second;
						const auto& burn_dir = dir_dist_pair.first;
						double di = dir_dist_pair.second;

						center_d += ( ( center - m_stSubdiv.getVert( qis[ j ] ) ).dot( burn_dir ) + di );
					}
					center_d /= 2.0f;

					Eigen::Vector2d c;
					/*c << center[0], center[1], center[2];*/
					c << 0.5f, center_d;

					linearSolveSVD<2>( A_, B_, c, sol_ );
				}
				else // weight closeness to center method
				{
					A_( 0, 0 ) += 2 * _w*P_.dot( P_ );
					Eigen::Vector3d c_;
					c_ << center[ 0 ], center[ 1 ], center[ 2 ];
					B_( 0, 0 ) -= 2 * _w*( p2_cpy - c_ ).dot( P_ );

					sol_ = A_.fullPivLu().solve( B_ );
				}

				// linear solve is done. construct poly dual point from the solution
				weights[ 0 ] = (float)sol_[ 0 ];
				weights[ 1 ] = 1 - weights[ 0 ];
				weights[ 2 ] = 0.0f;
				snap_to_boundary( weights );
				edge_dual_pos = weights[ 0 ] * m_stSubdiv.getVert( qis[ 0 ] ) +
					weights[ 1 ] * m_stSubdiv.getVert( qis[ 1 ] );
				dual_bt2 = sol_[ 1 ];
				dual_bt3 = weights[ 0 ] * bt3MA_vert[ qis[ 0 ] ] + weights[ 1 ] * bt3MA_vert[ qis[ 1 ] ];
			}
			else // new WEIGHT_DUAL
			{
				float f1_v1 = dir_dist1.second;
				float f1_v2 = f1_v1 +
					( this->m_stSubdiv.getVert( qis[ 1 ] ) - this->m_stSubdiv.getVert( qis[ 0 ] ) ).dot( dir_dist1.first );
				float f2_v2 = dir_dist2.second;
				float f2_v1 = f2_v2 +
					( this->m_stSubdiv.getVert( qis[ 0 ] ) - this->m_stSubdiv.getVert( qis[ 1 ] ) ).dot( dir_dist2.first );
				float b = std::max( f2_v1 - f1_v1, _eps );
				float a = std::max( f1_v2 - f2_v2, _eps );

				edge_dual_pos = ( a*this->m_stSubdiv.getVert( qis[ 0 ] ) + b*this->m_stSubdiv.getVert( qis[ 1 ] ) ) / ( a + b );
				te1 = this->mapTopo( qis[ 0 ], _fi );
				te2 = this->mapTopo( qis[ 1 ], _fi );
				dual_bt2 = ( a*bt2MA_vert_per_sheet[ qis[ 0 ] ][ te1 ] + b*bt2MA_vert_per_sheet[ qis[ 1 ] ][ te2 ] ) / ( a + b );
				dual_bt3 = ( a*bt3MA_vert[ qis[ 0 ] ] + b*bt3MA_vert[ qis[ 1 ] ] ) / ( a + b );
			}

			// now save the just created point in the edge_dual_map 
			// & the edge_duals list of cur poly
			auto it = _edge_dual_map.find( e );
			if ( it == _edge_dual_map.end() )
			{
				dual_vts.push_back( edge_dual_pos );
				m_is_face_dual.push_back( -1 ); // this is not a face dual
				dualV_idx = dual_vts.size() - 1;

				bt2_MC.push_back( dual_bt2 );
				bt3_MC.push_back( dual_bt3 );

				auto& tpl = _edge_dual_map[ e ];
				tpl.dual_v_id = dualV_idx;
				tpl.bt2 = dual_bt2;
				tpl.sum = edge_dual_pos;
				tpl.count = 1;
			}
			else
			{
				dualV_idx = it->second.dual_v_id;

				auto& tpl = _edge_dual_map[ e ];
				//tpl.dist_sum += dual_dist;
				tpl.sum += dual_vts[ dualV_idx ];
				tpl.count += 1;

				// store the largest bt2/bt3 at the intermediate structure for the dual vert
				// which will be finalized later
				tpl.bt2 = std::max( tpl.bt2, dual_bt2 );
				bt3_MC[ dualV_idx ] = std::max( bt3_MC[ dualV_idx ], dual_bt3 );
			}
			edge_dual_indices[ poly_ei ] = dualV_idx;
			edge_dual_bt2[ poly_ei ] = dual_bt2;
			edge_dual_bt3[ poly_ei ] = dual_bt3;

			// save the {<edge dual id, face id> -> the bt2 w.r.t. that face} relationship
			m_bt2_MC_perFace[ trimesh::ivec2( dualV_idx, _fi ) ] = dual_bt2;
		}// dual of edge

		int cnt_edge_duals = 0;
		for ( auto it = edge_no_dual.begin(); it != edge_no_dual.end(); ++it )
		{
			if ( *it == false )
				cnt_edge_duals++;
		}

		/* dualize poly region */
		TriPoint center_pos(0.0f, 0.0f, 0.0f);
		float center_bt2 = 0.0f;
		float center_bt3 = 0.0f;

		if ( _f_part_of_pocket ) // handle f that's part of pocket
		{
			for ( auto q_it = poly_it->begin(); q_it != poly_it->end(); ++q_it )
			{
				TriPoint q = m_stSubdiv.getVert( *q_it );
				center_pos = center_pos + q;
				center_bt3 += bt3MA_vert[ *q_it ];
			}
			poly_dual_pos = vec( center_pos ) / (float)poly_it->size();
			poly_dual_bt3 = center_bt3 / (float)poly_it->size();
			poly_dual_bt2 = infiniteBurnDist();
		}
		else // normal dualization of face
		{
			if ( cnt_edge_duals == 0 )
			{
				// fall-back plan to compute face dual and related info
				for ( auto q_it = poly_it->begin(); q_it != poly_it->end(); ++q_it )
				{
					TriPoint q = m_stSubdiv.getVert( *q_it );
					center_pos = center_pos + q;
					center_bt3 += bt3MA_vert[ *q_it ];
					auto si = mapTopo( *q_it, _fi );
					center_bt2 += bt2MA_vert_per_sheet[ *q_it ][ si ];
				}
				poly_dual_pos = vec( center_pos ) / (float)poly_it->size();
				poly_dual_bt3 = center_bt3 / (float)poly_it->size();
				poly_dual_bt2 = center_bt2 / (float)poly_it->size();
			}
			else
			{
				if ( _poly_dual_opt == SIMPLE_DUAL )
				{
					for ( unsigned ii = 0; ii < edge_dual_indices.size(); ++ii )
					{
						if ( edge_no_dual[ ii ] )
							continue;
						center_pos = center_pos + this->dual_vts[ edge_dual_indices[ ii ] ];
						center_bt2 = center_bt2 + edge_dual_bt2[ ii ];
						center_bt3 += edge_dual_bt3[ ii ];
					}
					poly_dual_pos = vec( center_pos ) / (float)cnt_edge_duals;
					poly_dual_bt2 = center_bt2 / (float)cnt_edge_duals;
					poly_dual_bt3 = center_bt3 / (float)cnt_edge_duals;
				}
				else if ( _poly_dual_opt == PSEUDO_INV_DUAL || _poly_dual_opt == WEIGHT_CENTER_CLOSENESS_DUAL )
				{
					A.setZero();
					B.setZero();

					for ( auto q_it = poly_it->begin(); q_it != poly_it->end(); ++q_it )
					{
						TriPoint q = m_stSubdiv.getVert( *q_it );
						qi << q[ 0 ], q[ 1 ], q[ 2 ];
						const auto& dir_dist_pair = _vert_burnDirDistPair_map.find( *q_it )->second;
						const auto& burn_dir = dir_dist_pair.first;
						double di = dir_dist_pair.second;

						dir << burn_dir[ 0 ], burn_dir[ 1 ], burn_dir[ 2 ];
						Ai << P.transpose() * dir, -1.0;
						A += Ai * Ai.transpose();
						Ci = ( p3_cpy - qi ).dot( dir ) + di;
						B += -1 * Ci * Ai;
					}

					// the center of the cur poly
					TriPoint center( 0.0f, 0.0f, 0.0f );
					for ( auto idx_iter = edge_dual_indices.begin(); idx_iter != edge_dual_indices.end(); ++idx_iter )
					{
						if ( edge_no_dual[ idx_iter - edge_dual_indices.begin() ] )
							continue;
						const TriPoint& q = dual_vts[ *idx_iter ];
						center += q;
					}
					center /= edge_dual_indices.size();

					if ( _poly_dual_opt == PSEUDO_INV_DUAL )
					{
						util::find_barycentric( center, p1, p2, p3, weights[ 0 ], weights[ 1 ], weights[ 2 ] );
						snap_to_boundary( weights );

						// find the dist of the center using the distance functions
						// defined by each poly vert
						float center_d = std::numeric_limits<float>::max();
						for ( auto idx_iter = poly_it->begin(); idx_iter != poly_it->end(); ++idx_iter )
						{
							const auto& dir_dist_pair = _vert_burnDirDistPair_map.find( *idx_iter )->second;
							const auto& burn_dir = dir_dist_pair.first;
							float di = dir_dist_pair.second;

							center_d = std::min(
								( center - m_stSubdiv.getVert( *idx_iter ) ).dot( burn_dir ) + di,
								center_d );
						}
						/*center_d /= poly_it->size();*/

						Eigen::Vector3d c;
						c << weights[ 0 ], weights[ 1 ], center_d;
						linearSolveSVD<3>( A, B, c, sol );
					}
					else // poly_dual_opt == WEIGHT_CENTER_CLOSENESS_DUAL
					{
						A.topLeftCorner( 2, 2 ) += _w*poly_it->size()*P.transpose()*P;
						Eigen::Vector3d c_;
						c_ << center[ 0 ], center[ 1 ], center[ 2 ];
						B.head( 2 ) += -1 * _w*poly_it->size()*( ( p3_cpy - c_ ).transpose() )*P;

						sol = A.fullPivLu().solve( B );
					}
					// now linear solve is done. construct poly dual from solution
					weights[ 0 ] = (float)sol[ 0 ]; // s
					weights[ 1 ] = (float)sol[ 1 ]; // t
					weights[ 2 ] = 1 - weights[ 0 ] - weights[ 1 ]; // 1-s-t
					snap_to_boundary( weights );
					poly_dual_pos = ( weights[ 0 ] ) * p1 + ( weights[ 1 ] ) * p2 + ( weights[ 2 ] ) * p3;
					poly_dual_bt2 = sol[ 2 ];
					poly_dual_bt3 =
						( weights[ 0 ] ) * bt3_MC[ f[ 0 ] ] +
						( weights[ 1 ] ) * bt3_MC[ f[ 1 ] ] +
						( weights[ 2 ] ) * bt3_MC[ f[ 2 ] ];
				}
				else // new WEIGHT_DUAL
				{
					// find the angle based weight for each dual on the edges of cur poly
					angle_weights.clear();
					angle_weights.reserve( edge_dual_indices.size() );
					for ( unsigned ii = 0; ii < poly_it->size(); ++ii )
					{
						// skip burnt edge
						e = TriEdge( ( *poly_it )[ ii ], ( *poly_it )[ ( ii + 1 ) % poly_it->size() ] );
						if ( edge_no_dual[ ii ] )
							continue;

						const auto& dir_dist_1 = _vert_burnDirDistPair_map.find( e[ 0 ] )->second;
						const auto& dir_dist_2 = _vert_burnDirDistPair_map.find( e[ 1 ] )->second;
						angle_weights.push_back(
							std::max(
								std::abs( trimesh::angle( dir_dist_1.first, dir_dist_2.first ) ), 0.1f*3.14f / 360.0f
							)
						);
					}
					assert( angle_weights.size() == cnt_edge_duals );

					// normalize the weights
					float sum = std::accumulate( angle_weights.begin(), angle_weights.end(), 0.0f );

					// find the weighted sum of edge dual positions as the poly dual
					// and the dual bt2;
					// simply average to find bt3;
					poly_dual_pos[ 0 ] = poly_dual_pos[ 1 ] = poly_dual_pos[ 2 ] = 0.0f;
					poly_dual_bt2 = 0.0f;
					poly_dual_bt3 = 0.0f;
					for ( unsigned ii = 0; ii < edge_dual_indices.size(); ++ii )
					{
						if ( edge_no_dual[ ii ] )
							continue;
						poly_dual_pos += ( angle_weights[ ii ] / sum ) * this->dual_vts[ edge_dual_indices[ ii ] ];
						poly_dual_bt2 += ( angle_weights[ ii ] / sum ) * edge_dual_bt2[ ii ];
						poly_dual_bt3 += edge_dual_bt3[ ii ];
					}
					poly_dual_bt3 /= (float)cnt_edge_duals;
				}
			}
		} // when _f_part_of_pocket == false, dualize f in the usual way

		// now actually append dual pos/dist of the poly region to output lists
		// and connect the poly dual with edge duals (and, if any, dual at orig vts)

		dual_vts.push_back(poly_dual_pos);
		this->m_is_face_dual.push_back(_fi); // this is a face dual
		this->bt2_MC.push_back(poly_dual_bt2);
		this->bt3_MC.push_back(poly_dual_bt3);

		unsigned poly_dual_idx = dual_vts.size()-1;
		m_bt2_MC_perFace[trimesh::ivec2(poly_dual_idx, _fi)] = poly_dual_bt2;

		for ( auto ii = 0; ii < edge_dual_indices.size(); ++ii )
		{
			if ( edge_no_dual[ ii ] )
			{
				// STEINER_ONLY scheme uses all orig vts as dual. connect to them if any.
				if ( getBurnScheme() == STEINER_ONLY )
				{
					auto poly_vi = ( *poly_it )[ ii ];
					if ( !m_stSubdiv.isSteinerVert( poly_vi ) )
					{
						// NOTE: the index of the dual is just the index of the orig vert
						dual_edges.push_back( util::makeEdge( poly_dual_idx, poly_vi ) );
					}
				}
			}
			else
			{
				this->dual_edges.push_back(
					util::makeEdge( poly_dual_idx, edge_dual_indices[ ii ] )
					);
			}
		}

		// triangulate the poly region using the st vts and dual vts.
		// case 1: connect a st, edge-dual, and a face-dual into a triangle.
		// case 2: connect 2 st vts and a face-dual into a triangle.
		// (dual vts index are offset by # st vts. to differentiate from st vts).
		// also, record for each dual edge which fine tri. it's from
		int index_offset = m_stSubdiv.sizeOfVts(); // purpose: dual index encoding
		for (unsigned poly_ei = 0; poly_ei < poly_it->size(); ++poly_ei)
		{
			bool is_dual_0 = false, is_dual_1 = false;
			if ( edge_no_dual[ poly_ei ] )
			{// this edge doesn't have a dual. 
				// need to split cases according to burn scheme
				if ( getBurnScheme() == STEINER_ONLY )
				{
					is_dual_0 = !m_stSubdiv.isSteinerVert( ( *poly_it )[ poly_ei ] ) ? true : false;
					is_dual_1 = !m_stSubdiv.isSteinerVert(
						( *poly_it )[ ( poly_ei + 1 ) % poly_it->size() ]
					) ? true : false;
				}
				m_finerTri_byDual.push_back( util::makeFace(
					( *poly_it )[ poly_ei ] + (is_dual_0 ? index_offset : 0),
					( *poly_it )[ ( poly_ei + 1 ) % poly_it->size() ] + (is_dual_1 ? index_offset : 0),
					poly_dual_idx + index_offset // dual index encoding
				) );
				m_origTri_for_finerTri_byDual.push_back( _fi );
			}
			else // this edge has a dual. common to both burn schemes
			{
				// two tri.s are formed.
				m_finerTri_byDual.push_back( util::makeFace(
					( *poly_it )[ poly_ei ],
					edge_dual_indices[ poly_ei ] + index_offset, // dual index encoding
					poly_dual_idx + index_offset // dual index encoding
					) );
				m_finerTri_byDual.push_back( util::makeFace(
					( *poly_it )[ ( poly_ei + 1 ) % poly_it->size() ],
					edge_dual_indices[ poly_ei ] + index_offset, // dual index encoding
					poly_dual_idx + index_offset // dual index encoding
					) );
				m_fromFineTri_for_dualE.push_back( m_finerTri_byDual.size() - 1 );
				m_origTri_for_finerTri_byDual.push_back( _fi );
				m_origTri_for_finerTri_byDual.push_back( _fi );
			}
		}

		// TODO: this will not be needed when the H.S. is built on the fine tri.s directly.
		// collect each edge dual and poly dual into the H.S. structure 
		//unsigned k = 0;
		//for (unsigned poly_ei = 0; poly_ei < poly_it->size(); ++poly_ei)
		//{
		//	// skip if cur edge doesn't have a dual
		//	if ( edge_no_dual[poly_ei] )
		//		continue;

		//	// else add the edge dual & face dual to the H.S.
		//	qis[0] = *(poly_it->begin()+poly_ei);
		//	qis[1] = *(poly_it->begin()+(poly_ei+1)%poly_it->size());
		//	TriEdge e = util::makeEdge( qis[0], qis[1] );
		//	//this->recordDualwrtStEdge(_fi, e, edge_dual_indices[k], poly_dual_idx);
		//	k ++ ;
		//}

	} // for each poly face
} // SteinerGraph::dualize_poly_subdivision()

void SteinerGraph::process_orig_vts_as_duals()
{
	auto setup_for_adding_dual = [ & ]( int _size )
	{
		dual_vts.resize( _size );
		bt2_MC.resize( _size );
		bt3_MC.resize( _size );
		m_is_face_dual.resize( _size );
	};
	// add orig vert _old_vi to position _new_vi within dual vts list
	auto add_as_dual_vert = [ & ]( int _old_vi, int _new_vi )
	{
		dual_vts[ _new_vi ] = m_origG->vts[ _old_vi ];
		bt2_MC[_new_vi] = bt2MA_vert[ _old_vi ];
		bt3_MC[_new_vi] = bt3MA_vert[ _old_vi ];
		m_is_face_dual[ _new_vi ] = -2;
		// save per-face bt2 info for the dual
		const auto& nb_fs = m_origG->nbFacesOfVert[ _old_vi ];
		if ( nb_fs.empty() ) // isolate vts, or vts on isolate edge
			m_bt2_MC_perFace[ trimesh::ivec2( _new_vi, -1 ) ] = bt2MA_vert[ _old_vi ];
		else // with regular neighborhood
			for ( auto j = 0; j < nb_fs.size(); ++j )
			{
				auto fi = nb_fs[ j ];
				auto si = mapTopo( _old_vi, fi );
				float bt2 = bt2MA_vert_per_sheet[ _old_vi ][ si ];
				m_bt2_MC_perFace[ trimesh::ivec2( _new_vi, fi ) ] = bt2;
			}
	};
	auto isolate_edges = m_origG->getIsolateEdges();
	if ( getBurnScheme() == STEINER_ONLY )
	{
		// add all orig vts as dual vts
		setup_for_adding_dual( m_origG->vts.size() );
		for ( auto vi = 0; vi < m_origG->vts.size(); ++vi )
			add_as_dual_vert( vi, vi );
		// add isolated edges as dual
		auto pre_size = dual_edges.size();
		dual_edges.insert( dual_edges.end(), isolate_edges.begin(), isolate_edges.end() );
		for ( auto i = pre_size; i < dual_edges.size(); ++i )
			m_dualEdge_from_which_face.push_back( -1 );
	}
	else if ( getBurnScheme() == ORIGINAL_AND_STEINER )
	{
		// only add vts that are relevant to isolate edges as dual
		map<int, int> old_new_map;
		size_t presize;
		for ( auto i = 0; i < isolate_edges.size(); ++i )
		{
			auto e = isolate_edges[ i ];
			presize = old_new_map.size();
			auto& new_id0 = old_new_map[ e[ 0 ] ];
			if ( old_new_map.size() > presize )
				new_id0 = presize;
			presize = old_new_map.size();
			auto& new_id1 = old_new_map[ e[ 1 ] ];
			if ( old_new_map.size() > presize )
				new_id1 = presize;
		}
		setup_for_adding_dual( old_new_map.size() );
		for ( auto it = old_new_map.begin(); it != old_new_map.end(); ++it )
			add_as_dual_vert( it->first, it->second );
		// add isolated edges as dual
		auto pre_size = dual_edges.size();
		for ( auto i = 0; i < isolate_edges.size(); ++i )
		{
			auto e = isolate_edges[ i ];
			dual_edges.push_back( util::makeEdge( old_new_map[ e[ 0 ] ], old_new_map[ e[ 1 ] ] ) );
		}
		for ( auto i = pre_size; i < dual_edges.size(); ++i )
			m_dualEdge_from_which_face.push_back( -1 );
	}
}

void SteinerGraph::find_burn_dir_dist_of_face(
	const vector<unsigned>& _vts, unsigned _fi, 
	const vector<trimesh::vec>& _bndry_dirs, const map<unsigned, unsigned>& _vIdx_burnDirIdx_map, 
	map<unsigned, std::pair<trimesh::vec, float>>& _vert_burnDirDistPair_map,
	bool _verbose )
{
	auto f = this->m_origG->faces[_fi];
	TriEdge e;
	for (auto v_it = _vts.begin(); v_it != _vts.end(); ++v_it)
	{
		if ( _verbose )
			cout << "vid: " << *v_it <<". ";

		// skip vert that cannot be reached by fire (logic differs by burn-scheme)
		if ( getBurnScheme() == STEINER_ONLY && !m_stSubdiv.isSteinerVert( *v_it ) )
			continue;

		trimesh::vec burn_dir;
		// need to find out the burn dist corresponding to cur face _fi by examining topo info
		unsigned vi2;
		auto tt = getTopoType(*v_it);
		TopoGraph::EdgeIdx tei;
		/*if (isSteinerVert(*v_it))
		{
		e = m_stSubdiv.getResidingEdge(*v_it);
		vi2 = util::oppositeVert(f, e);

		// using e2 to find the topo edge of v_it that covers the cur face
		// so that we can use that information to get prev and corresponding dist
		TriEdge e2 = util::makeEdge(*v_it, vi2);
		tei = (TopoGraph::EdgeIdx)(tg.t_origEdge2TopoMap.find(e2)->second.second);
		}
		else // v is an orig. vert of face _fi
		{
		bool success = getAssocTEdgeOfVertForFace(*v_it, f, tei);
		assert(success);
		}*/
		int ret_code = mapTopo(*v_it, _fi);
		assert(ret_code >= 0);
		tei = (TopoGraph::EdgeIdx)ret_code;
		float burn_dist = this->bt2MA_vert_per_sheet[ *v_it ][ tei ];

		ret_code = getPrevVert(*v_it, tei);
		assert(ret_code >= 0);
		unsigned pre_v = (unsigned)ret_code;
		if ( _verbose )
		{
			cout << topo_type_str[ tt ] << ". "
				<< "tei " << (int)tei <<". "
				<< "pre_v " << pre_v << ". ";
		}
		// and the burn direction
		if ( TopoGraph::isOpenSheet( *v_it, tei ) )
		{
			auto find_it = _vIdx_burnDirIdx_map.find( *v_it );
			if ( find_it == _vIdx_burnDirIdx_map.end() )
				std::cout << "Error: couldn't find the burn-dir-id for vert-id " << *v_it << endl;
			unsigned dir_idx = find_it->second;
			burn_dir = _bndry_dirs.at(dir_idx);
		}
		else
		{
			burn_dir = m_stSubdiv.getVert(*v_it) - m_stSubdiv.getVert(pre_v);
			trimesh::normalize(burn_dir);
		}
		// TODO: this burn_dir needs to be rotated so that it lies in the plane of cur face f

		// add them to map
		_vert_burnDirDistPair_map[*v_it] = std::make_pair(burn_dir, burn_dist);
		if ( _verbose )
		{
			cout << "done." << endl;
		}
	}
}

void SteinerGraph::finalize_dual_vts(
	map<TriEdge, movingSum_edgeDual_struct>& _edge_dual_map)
{
	m_dualV_triEdgeIdx.assign(dual_vts.size(), -1);

	// iterate thru all edge duals
	for (auto iter = _edge_dual_map.begin(); iter != _edge_dual_map.end(); ++iter)
	{
		auto& tpl = iter->second;
		unsigned dual_idx = tpl.dual_v_id;
		dual_vts[dual_idx] = (1.0f / tpl.count) * tpl.sum;
		//dist2bndry_medialCurve[dual_idx] = (1.0f / tpl.count) * tpl.dist_sum;
		bt2_MC[dual_idx] = tpl.bt2;

		// complete the corresp: dual vert -> tri edge it's on
		m_dualV_triEdgeIdx[dual_idx] = 
			m_stSubdiv.isSteinerVert(iter->first[0]) ? m_stSubdiv.getResidingEdgeIdx(iter->first[0]) : 
			m_stSubdiv.isSteinerVert(iter->first[1]) ? m_stSubdiv.getResidingEdgeIdx(iter->first[1]) : 
			m_origG->getEdgeIdx(iter->first);
	}

	// ***IMPORTANT SUBTLE LOGIC:*** //
	// b.c. on singular curve, a st edge is only dualizable on a subset of faces,
	// there are invalid triangulation induced by duals 
	// that changes the topo around that st edge.
	// fix those invalid triangles.
	int temp0, temp1, from, to;
	vector<TriFace> new_fineTri;
	vector<int> new_origTri_for_finerTri_byDual;
	set<TriEdge> padded;
	for ( auto fi = 0; fi < m_finerTri_byDual.size(); ++fi )
	{
		const auto& f = m_finerTri_byDual[ fi ];
		auto from_tri = m_origTri_for_finerTri_byDual[ fi ];
		for ( int i = 0; i < 3; ++i )
		{
			auto e = util::makeEdge( f[ i ], f[ ( i + 1 ) % 3 ] );
			//auto oppo_v = f[ ( i + 2 ) % 3 ];
			if ( ( isDualVertInFineTri( e[ 0 ], temp0 ) || isDualVertInFineTri( e[ 1 ], temp1 ) ) /*||
				( getBurnScheme() == STEINER_ONLY && 
				( m_stSubdiv.isSteinerVert( temp0 ) || m_stSubdiv.isSteinerVert( temp1 ) ) )*/
				)
				continue;
			if ( isBurnPath( e, from_tri, from, to ) /*|| !burnt[ e[ 0 ] ] && !burnt[ e[ 1 ] ]*/ )
			{
				// there might be a dual created for this st edge on other faces;
				// for each such face, add a padding triangle at the st edge on that face.
				auto find_it = _edge_dual_map.find( e );
				if ( find_it != _edge_dual_map.end() /*&& padded.count( e ) == 0*/ )
				{
					// add a padding tri face using the dual
					const auto& tpl = find_it->second;
					auto dual_to_store = tpl.dual_v_id + m_stSubdiv.sizeOfVts();
					TriFace f1 = util::makeFace( e[ 0 ], e[ 1 ], dual_to_store );
					//TriFace f2 = util::makeFace( oppo_v, e[ 1 ], tpl.dual_v_id );
					new_fineTri.push_back( f1 );
					new_origTri_for_finerTri_byDual.push_back( from_tri );
					padded.insert( e );
				}
			}
		}
	}
	m_finerTri_byDual.insert( m_finerTri_byDual.end(), new_fineTri.begin(), new_fineTri.end() );
	m_origTri_for_finerTri_byDual.insert( m_origTri_for_finerTri_byDual.end(),
		new_origTri_for_finerTri_byDual.begin(), new_origTri_for_finerTri_byDual.end() );
}

bool SteinerGraph::is_edge_dualizable(
	unsigned _u, unsigned _v, unsigned _fi,
	const map<TriEdge, vector<int> >& _burntEdge_assocFaces_map,
	bool _verbose )
{
	TriEdge e = TriEdge( _u, _v );
	if ( _verbose )
		cout << "edge " << e << ": ";

	//// TODO: this check cannot differentiate between two possible situations.
	//// so get rid of this. 
	//// first, this edge could be part of pocket, then not dualizable
	//if ( !burnt[ e[ 0 ] ] || !burnt[ e[ 1 ] ] )
	//{
	//	if ( _verbose ) cout << "is part of pocket, thus not dualizable." << endl;
	//	return false;
	//}

	if (
		getBurnScheme() == STEINER_ONLY &&
		( !m_stSubdiv.isSteinerVert( e[ 0 ] ) || !m_stSubdiv.isSteinerVert( e[ 1 ] ) )
		)
	{
		if ( _verbose ) cout << "already contains at least one dual vert, not dualizable." << endl;
		return false;
	}

	list<int> assoc_faces;
	auto edge_face_pair = _burntEdge_assocFaces_map.find( e );
	if ( edge_face_pair != _burntEdge_assocFaces_map.end() )
		assoc_faces.insert(assoc_faces.end(), edge_face_pair->second.begin(), edge_face_pair->second.end());
	edge_face_pair = _burntEdge_assocFaces_map.find( TriEdge(e[1], e[0]) );
	if ( edge_face_pair != _burntEdge_assocFaces_map.end() )
		assoc_faces.insert(assoc_faces.end(), edge_face_pair->second.begin(), edge_face_pair->second.end());

	// no dual point is created
	// only if the edge is used by a burn path 
	// and the assoc face id matches cur face
	bool part_of_burntree = false;
	int assoc_f = -2; // by default means this edges is not a burnt edge
	for (auto a_f = assoc_faces.begin(); a_f != assoc_faces.end(); ++a_f)
	{
		assoc_f = *a_f;
		if ( *a_f == _fi || *a_f == -1 )
		{
			if (_verbose)
				cout << "assoc_f " << assoc_f <<". occupied & skip."<<endl;
			part_of_burntree = true;
			break;
		}
	}

	if (!part_of_burntree && _verbose)
		cout <<"assoc_f "<< assoc_f <<". is dualizable." << endl;

	return !part_of_burntree;
} // SteinerGraph::is_burnt_edge()

void SteinerGraph::computeIntersectionsForFaces(
	const std::shared_ptr<TriMesh> _mesh_orig,
	const vector<float>& _dist2surf
	)
{
	typedef double ldouble;
	typedef Vec<3,ldouble> DPoint;
	typedef Vec<3,ldouble> DVec;

	//bt3MA_face.clear();
	/*if (_dist2surf.empty())
		return;*/
	this->m_footPtsForMAFaces.resize(2*m_origG->faces.size());
	this->has_isects.resize(m_origG->faces.size(), true);

	TriFace f;
	DPoint a, b, c;
	ldouble ra, rb, rc;
	int non_isct_cnt = 0;
	int pt_isct_cnt = 0;
	int circ_isct_cnt = 0;
	int sph_isct_cnt = 0;
	float max_float = numeric_limits<float>::max();

	unsigned same_isect_cnt = 0;

	// build a kdtree of surface vts in case some faces don't have valid intersections
	KNNTree* tree = new KNNTree();
	TriPoint p;
	for (unsigned i = 0; i < _mesh_orig->vertices.size(); ++i)
	{
		p = _mesh_orig->vertices[i];
		tree->insert( Point_and_uint(P3(p[0], p[1], p[2]), i) );
	}

	TriPoint two_foot_pts[2]; // will contain 2 foot points returned by specialized routine
	Vec<3, double> foot_pts_double[2];
	array<DPoint, 2> isecs; // will contain 2 intersection points if successful
	for (unsigned fi = 0; fi < m_origG->faces.size(); ++fi)
	{
		f = m_origG->faces[fi];
		a = DPoint((ldouble)m_origG->vts[f[0]][0], (ldouble)m_origG->vts[f[0]][1], (ldouble)m_origG->vts[f[0]][2]);
		b = DPoint((ldouble)m_origG->vts[f[1]][0], (ldouble)m_origG->vts[f[1]][1], (ldouble)m_origG->vts[f[1]][2]);
		c = DPoint((ldouble)m_origG->vts[f[2]][0], (ldouble)m_origG->vts[f[2]][1], (ldouble)m_origG->vts[f[2]][2]);
		if (!_dist2surf.empty())
		{
			ra = _dist2surf[f[0]];
			rb = _dist2surf[f[1]];
			rc = _dist2surf[f[2]];
		}
		else
		{
			ra = 0.0f;
			rb = 0.0f;
			rc = 0.0f;
		}

		util::IsectType isec_type;
		/*isec_type = threeSphereIntersect(
		a, b, c, ra, rb, rc, isecs);*/

		// use cgal intersector
		//vector<CGAL::Object> isecs_cgal;
		//::IsectType isec_type = ::cgalSphereIntersect(
		//	a, b, c, ra, rb, rc, isecs_cgal
		//	);
		//if (isec_type == ::POINT_ISECT)
		//{
		//	if ( const ::PointRetType* point_ret = 
		//		CGAL::object_cast<::PointRetType>(&isecs_cgal[0]) )
		//	{
		//		isecs[0] = Vec<3, double>(
		//			CGAL::to_double(point_ret->first.x()), 
		//			CGAL::to_double(point_ret->first.y()), 
		//			CGAL::to_double(point_ret->first.z()) );
		//		if (point_ret->second == 2)
		//		{
		//			isecs[1] = isecs[0];
		//			same_isect_cnt ++;
		//		}
		//		else
		//		{
		//			point_ret = CGAL::object_cast<::PointRetType>(&isecs_cgal[1]);
		//			isecs[1] = Vec<3, double>(
		//				CGAL::to_double(point_ret->first.x()), 
		//				CGAL::to_double(point_ret->first.y()), 
		//				CGAL::to_double(point_ret->first.z()) );
		//		}
		//	}
		//}

		/*isec_type = util::computeFootPoints(
		a, b, c, ra, rb, rc, foot_pts_double
		);*/
		isec_type = find_two_foot_points_for_MA_face(tree, fi, two_foot_pts);

		switch(isec_type) // debug
		{
		case util::NONE_ISECT: 
			non_isct_cnt++; break;
		case util::POINT_ISECT:
			pt_isct_cnt++; break;
		case util::CIRCLE_ISECT:
			circ_isct_cnt++; break;
		case util::SPHERE_ISECT:
			sph_isct_cnt++; break;
		}
		//cout << "isect type: "<<isec_type<<endl;

		if ( isec_type != util::POINT_ISECT )
		{		
			m_footPtsForMAFaces[fi*2 + 0] = TriPoint(max_float, max_float, max_float);
			m_footPtsForMAFaces[fi*2 + 1] = TriPoint(max_float, max_float, max_float);

			// use a special routine to find two foot points
			/*find_two_foot_points_for_MA_face(tree, fi, two_foot_pts);
			isects[fi*2 + 0] = two_foot_pts[0];
			isects[fi*2 + 1] = two_foot_pts[1];*/

			has_isects[fi] = false;
		}
		else
		{
			/*isects[fi*2 + 0] = TriPoint(foot_pts_double[0]);
			isects[fi*2 + 1] = TriPoint(foot_pts_double[1]);*/
			m_footPtsForMAFaces[fi*2 + 0] = TriPoint(two_foot_pts[0]);
			m_footPtsForMAFaces[fi*2 + 1] = TriPoint(two_foot_pts[1]);

			has_isects[fi] = true;

		}

		// use a special routine to find two foot points
		/*bool valid_foot_pts = find_two_foot_points_for_MA_face(tree, fi, two_foot_pts);
		if ( valid_foot_pts )
		{
		(*_isects_per_face)[fi*2 + 0] = two_foot_pts[0];
		(*_isects_per_face)[fi*2 + 1] = two_foot_pts[1];
		}
		else
		{
		(*_isects_per_face)[fi*2 + 0] = TriPoint(max_float, max_float, max_float);
		(*_isects_per_face)[fi*2 + 1] = TriPoint(max_float, max_float, max_float);
		}*/
	}

	delete tree;
}

util::IsectType SteinerGraph::find_two_foot_points_for_MA_face( const KdTree* _tree, unsigned _fi, TriPoint* _pts )
{
	/* compute centroid of the face */
	const auto& f = m_origG->faces[_fi];
	TriPoint cent(0.0f,0.0f,0.0f);
	cent += m_origG->vts[f[0]];
	cent += m_origG->vts[f[1]];
	cent += m_origG->vts[f[2]];
	cent /= 3.0f;
	/*cout << "face "<<_fi<<": "
	<<m_origG->vts[f[0]]<<", "
	<<m_origG->vts[f[1]]<<", "
	<<m_origG->vts[f[2]]<<endl;*/

	/* use the centroid to find the closest point in kdtree, i.e. our first foot point */  
	P3 query_p(cent[0], cent[1], cent[2]);
	K_neighbor_search search1( *_tree, query_p, 1 );
	auto pr = search1.begin()->first;
	auto v_on_surf = boost::get<0>(pr);
	auto p = TriPoint(v_on_surf[0], v_on_surf[1], v_on_surf[2]);
	//cout << "1st closest point: "<<p<<endl;
	auto vec_to_p = p - cent;
	auto dir_to_p = vec_to_p;
	trimesh::normalize(dir_to_p);
	p = cent + bt3MA_face[_fi]*dir_to_p;
	vec_to_p = p - cent;
	_pts[0] = p;

	/* reflect the closest point to get the second foot point */
	// is the face quality good enough to be used as reflectance plane? 
	bool bad_face = util::is_degenerate(
		m_origG->vts[f[0]],
		m_origG->vts[f[1]],
		m_origG->vts[f[2]],
		0.0f, 0.001f
		);
	if ( bad_face ) // normal unusable.
	{
		_pts[1] = _pts[0];
	}
	else // good face, normal usable.
	{
		auto nml = trimesh::normalize( trimesh::trinorm(
			m_origG->vts[f[0]],
			m_origG->vts[f[1]],
			m_origG->vts[f[2]]
		) );
		auto mirror_len = 2*vec_to_p.dot(nml);
		// locate the second closest
		auto q = p - mirror_len * nml; 

		/*if (_fi == 45)
		{
		cout << "foot info of face "<<_fi<<":"<<endl;
		cout << "cent: "<<cent<<endl;
		cout << "normal: "<<nml<<endl;
		cout << "vec_to_p: " <<vec_to_p<<endl;
		cout << "mirror_len: "<<mirror_len<<endl;
		}*/
		_pts[1] = q;
	}

	//cout << "foot points: "<<_pts[0]<<", "<<_pts[1]<<endl;
	return bad_face ? util::NONE_ISECT : util::POINT_ISECT;
}

void SteinerGraph::output_footPts_info( 
	const vector<TriPoint>& _foot_pts,
	const vector<bool>& _has_valid_isects
	)
{
	std::stringstream out_string;
	out_string << "{"<< endl; // Big {

	// Begin: output mesh vts
	out_string << "{ (*Begin: mesh vts*)"<< endl; 
	for (unsigned i = 0; i < m_origG->vts.size(); ++i)
	{
		const auto & p = m_origG->vts[i];
		out_string << p[0]<<","<<p[1]<<","<<p[2];
		if (i != m_origG->vts.size()-1)
		{
			out_string <<",";
		}
	}
	out_string << endl;
	out_string << "}, (*End: mesh vts*)"<< endl; 
	// End: output mesh vts

	// Begin: output mesh faces
	out_string << "{ (*Begin: mesh faces*)"<< endl; 
	for (unsigned i = 0; i < m_origG->faces.size(); ++i)
	{
		const auto & f = m_origG->faces[i];
		out_string << f[0]+1<<","<<f[1]+1<<","<<f[2]+1;
		if (i != m_origG->faces.size()-1)
		{
			out_string <<",";
		}
	}
	out_string << endl;
	out_string << "}, (*End: mesh faces*)"<< endl; 
	// End: output mesh faces

	// Begin: output MA radii
	out_string << "{ (*Begin: MA radii*)"<< endl; 
	unsigned num_origVts = m_origG->vts.size();
	for (unsigned i = 0; i < num_origVts; ++i)
	{
		out_string << bt3MA_vert[i] << ",";
	}
	out_string.seekp(-1, out_string.cur);
	out_string << "}, (*End: MA radii*)"<< endl; 
	// End: output MA radii

	// Begin: output whether face has valid foot points or not
	out_string << "{ (*Begin: has valid intersections or not*)"<< endl; 
	unsigned num_faces = m_origG->faces.size();
	for (unsigned i = 0; i < num_faces; ++i)
	{
		out_string << (_has_valid_isects[i] ? 1 : 0) << ",";
	}
	out_string.seekp(-1, out_string.cur);
	out_string << "}, (*End: has valid intersections or not*)"<< endl; 
	// End: output whether face has valid foot points or not

	// Begin: output foot points for each face
	out_string << "{ (*Begin: foot points (2 for each face) *)"<< endl; 
	for ( int i=0; i < num_faces; ++i )
	{
		for (int j = 0; j < 2; ++j)
		{
			auto p = _foot_pts[2*i+j];
			out_string << "{"<<p[0]<<","<<p[1]<<","<<p[2]<<"}"<<",";
		}
	}
	out_string.seekp(-1, out_string.cur);
	out_string << "} (*End: foot points (2 for each face) *)"<< endl; 
	// End: output foot points for each face

	out_string << "}" << endl;// Big }

	ofstream o_file("footPtsInfo.txt"); 
	o_file << out_string.rdbuf();
	o_file.close();	
}

void SteinerGraph::setupDiffuserForMA( GraphDiffuser& _diffuser )
{
	/* faces as nodes of the graph */
	vector<TriPoint> g_vts(m_origG->faces.size());
	for ( unsigned fi = 0; fi < m_origG->faces.size(); ++fi )
	{
		const auto& f = m_origG->faces[fi];
		g_vts[fi] =  
			(m_origG->vts[f[0]] + m_origG->vts[f[1]] + m_origG->vts[f[2]]) / 3.0f;
	}

	/* a graph edge is connected between two neighbor faces sharing an edge e
	whenever they are finalized together w.r.t. at least one steiner vert on e*/
	vector<vector<int>> g_adj(m_origG->faces.size());
	vector<unsigned> st_vts_e;
	vector<int> assoc_faces_i, assoc_faces_j;
	set<TriEdge> g_edge_union;
	for ( unsigned ei = 0; ei < m_origG->edges.size(); ++ei )
	{
		/*if (ei % 5000 == 0)
			cout << "# MA edges processed: "<<ei << endl;*/

		st_vts_e.clear();
		g_edge_union.clear();
		this->m_stSubdiv.getStVertIndicesOnTriEdge(ei, st_vts_e);

		for ( unsigned j = 0; j < st_vts_e.size(); ++j )
		{
			auto st_v = st_vts_e[j];
			auto tt = getTopoType(st_v);

			if ( tt == TopoType::MANIFOLD_2D )
			{
				// in this case, only two faces assoce.ed with the sole sheet
				// simply save them
				assoc_faces_i.clear();
				getAssocFaces(st_v, 0, assoc_faces_i);
				assert(assoc_faces_i.size() == 2);
				g_edge_union.insert( util::makeEdge(assoc_faces_i[0], assoc_faces_i[1]) );
			}
			else if ( tt == TopoType::NONMANIFOLD_2D )
			{
				for ( TopoGraph::EdgeIdx si = 0; si < t_n_sector[st_v] - 1; ++si )
				{
					if ( !is_leading_topo_edge(st_v, si) )
						continue;
					assoc_faces_i.clear();
					getAssocFaces(st_v, si, assoc_faces_i);
					// should have only one assoc. face per sheet for a steiner vert
					auto fi = assoc_faces_i[0]; 
					for ( TopoGraph::EdgeIdx sj = si+1; sj < t_n_sector[st_v]; ++sj )
					{
						if ( !is_leading_topo_edge(st_v, sj) && 
							prev_tedge[st_v][sj] == si
							)
						{
							assoc_faces_j.clear();
							getAssocFaces(st_v, sj, assoc_faces_j);
							// again should have only one assoc. face 
							auto fj = assoc_faces_j[0];
							// save to the union
							g_edge_union.insert( util::makeEdge(fi, fj) );
						}
					}
				}
			}
		}

		// save edges to graph edge list
		for ( auto it = g_edge_union.begin(); it != g_edge_union.end(); ++it )
		{
			const auto& g_e = *it;
			g_adj[g_e[0]].push_back( g_e[1] );
			g_adj[g_e[1]].push_back( g_e[0] );
		}
	}

	vector<vector<TriEdge>> cc_elist;
	vector<vector<unsigned>> cc_vlist;
	this->find_conn_cmpnts(g_vts.size(), g_adj, cc_elist, cc_vlist);
	for (auto cc_it = cc_vlist.begin(); cc_it != cc_vlist.end(); ++cc_it)
	{
		// if this cc has no known value, enforce one.
		bool has_seed = false;
		for (unsigned i = 0; i < cc_it->size(); ++i)
		{
			auto vi = *(cc_it->begin()+i);
			if ( has_isects[ vi ] )
			{
				has_seed = true;
				break;
			}
		}
		if ( !has_seed ) 
		{
			// enforce seed...
			has_isects[*cc_it->begin()] = true;
		}
	}
	cc_elist.clear();
	cc_vlist.clear();

	/* setup the diffuser */
	_diffuser.setupDiffuseSystem(g_vts, g_adj, UniformWeight(), has_isects);
}

void SteinerGraph::computeDist2SurfForFaces(const vector<float>& _dist2surf)
{
	bt3MA_face.resize(m_origG->faces.size());
	TriFace f;
	TriPoint a, b, c;
	float ra, rb, rc;
	TriPoint isec1, isec2;
	for (unsigned fi = 0; fi < m_origG->faces.size(); ++fi)
	{
		f = m_origG->faces[fi];
		a = m_origG->vts[f[0]];
		b = m_origG->vts[f[1]];
		c = m_origG->vts[f[2]];
		ra = _dist2surf.empty() ? 0.0f : _dist2surf[f[0]];
		rb = _dist2surf.empty() ? 0.0f : _dist2surf[f[1]];
		rc = _dist2surf.empty() ? 0.0f : _dist2surf[f[2]];

		//if ( _isects_per_face[fi*2 + 0][0] == numeric_limits<float>::max() )
		//{
		//	dist2surf_face[fi] = (ra + rb + rc) / 3.0f;
		//}
		//else
		//{
		//	// distance from intersection (point on surface) to face center
		//	dist2surf_face[fi] = trimesh::dist(
		//		(a + b + c ) / 3.0f, _isects_per_face[fi*2 + 0]
		//	);
		//}

		bt3MA_face[fi] = (ra + rb + rc) / 3.0f;
	}
}

void SteinerGraph::computeDist2SurfForFaces(const vector<TriPoint>& _closest_2pts)
{
	TriFace f;
	TriPoint center;
	bt3MA_face.resize(m_origG->faces.size());
	for (unsigned fi = 0; fi < this->m_origG->faces.size(); ++fi)
	{
		f = m_origG->faces[fi];
		center = 1.0f / 3.0f * (m_origG->vts[f[0]] + m_origG->vts[f[1]] + m_origG->vts[f[2]]);
		//cout << "2*fi: " << 2*fi << ". _closest size: " << _closest_2pts.size() << endl; //debug
		bt3MA_face[fi] = trimesh::dist(center, _closest_2pts[2*fi]);
		//cout << "_closest[2*fi][0]="<<_closest_2pts[2*fi][0]<< endl; // debug
	}
}

void SteinerGraph::computeAngleMetricForFaces(
	const vector<TriPoint>& _isects_per_face, 
	vector<float>& _angle_per_face
	)
{
	_angle_per_face.assign(m_origG->faces.size(), -1.0f);
	TriFace f;
	TriPoint center;
	vec u, v;
	float cos_uv;
	for (unsigned fi = 0; fi < m_origG->faces.size(); ++fi)
	{
		f = m_origG->faces[fi];
		if ( _isects_per_face[2*fi + 0][0] == numeric_limits<float>::max() )
		{
			_angle_per_face[fi] = -1.0f; // invalid value
			//_angle_per_face[fi] = 0.0f;
		}
		else
		{
			center = (
				m_origG->vts[f[0]] + 
				m_origG->vts[f[1]] + 
				m_origG->vts[f[2]]
			) / 3.0f;

			u = trimesh::normalize( vec(_isects_per_face[2*fi + 0] - center) );
			v = trimesh::normalize( vec(_isects_per_face[2*fi + 1] - center) );
			cos_uv = u DOT v;
			cos_uv = std::min(1.0f, std::max(-1.0f, cos_uv));
			_angle_per_face[fi] = acos(cos_uv) * 180.0f / M_PI; // convert to degrees
		}
	}

	/// TODO: factor this out & finish this
	// guess a valid value for each face that gets a invalid value above
	// via interpolation
	/*vector<float> new_values;
	vector<unsigned> face_idx;
	struct FaceCompare
	{
	bool operator () (const TriFace& a, const TriFace& b)
	{
	if (a[0] == b[0])
	if (a[1] == b[1])
	return a[2] < b[2];
	else
	return a[1] < b[1];
	else
	return a[0] < b[0];
	}
	};
	set<TriFace, FaceCompare> nbFaces_set;
	for (unsigned i = 0; i < _angle_per_face.size(); ++i)
	{
	if ( !(_angle_per_face[i] < 0.0f) )
	continue;

	face_idx.push_back(i);
	nbFaces_set.clear();
	f = orig_g.faces[i];
	for (unsigned ii = 0; ii < 3; ++ii)
	{
	auto* nbFaces = &(orig_g.nbFacesOfVert[f[ii]]);
	for (unsigned ni = 0; ni < nbFaces->size(); ++ni)
	{
	if ( _angle_per_face[(*nbFaces)[ni]] < 0.0f )
	continue;
	TriFace nf = orig_g.faces[(*nbFaces)[ni]];
	nbFaces_set.insert(nf);
	}
	}
	}*/
}

void SteinerGraph::computeLambdaMetricForFaces(
	const vector<TriPoint>& _isects_per_face, 
	vector<float>& _lambda_per_face)
{
	_lambda_per_face.assign(m_origG->faces.size(), -1.0f);
	TriFace f;
	TriPoint center;
	vec u, v;
	float cos_uv;
	for (unsigned fi = 0; fi < m_origG->faces.size(); ++fi)
	{
		f = m_origG->faces[fi];
		if ( _isects_per_face[2*fi + 0][0] == numeric_limits<float>::max() )
		{
			_lambda_per_face[fi] = -1.0f; // invalid value
			//_lambda_per_face[fi] = 2.0f*this->bt3MA_face[fi];
		}
		else
		{
			_lambda_per_face[fi] = 
				trimesh::dist(_isects_per_face[2*fi + 0], _isects_per_face[2*fi + 1]);
		}
	}
}

void SteinerGraph::computeGeodMetricForFaces(
	shared_ptr<TriMesh> _m3d, // the orig 3d mesh
	const vector<TriPoint>& _isects_per_face, 
	vector<float>& _geoD_per_face, 
	string& _cache_file)
{
	// try opening the cache file and read geodesic distances
	// if failed to open, recompute distances
	cout << "trying to find cache file "<<_cache_file<<"...\n";
	ifstream in_file;
	in_file.open(_cache_file);
	if (in_file.is_open())
	{
		cout << "found a cache file "<<_cache_file<<endl;
		cout << "reading in geo.d. ..." << endl;
		// read in distances
		_geoD_per_face.resize(m_origG->faces.size());
		in_file.exceptions(std::ifstream::failbit | std::ifstream::badbit);
		try
		{
			for (unsigned fi = 0; fi < m_origG->faces.size(); ++fi)
			{
				if ( _isects_per_face[fi*2 + 0][0] == numeric_limits<float>::max() )
				{
					float dummy;
					in_file >> dummy;
					_geoD_per_face[fi] = -1;
				}
				else
				{
					in_file >> _geoD_per_face[fi];
				}
			}
		}
		catch (std::ifstream::failure e)
		{
			std::cerr << "Error occurred while reading file "<<_cache_file<<endl;
		}

		cout << "geo.d. reading done!" <<endl;
		in_file.close();
	}
	else
	{
		cout << "couldn't find cache file. recomputing... " << endl;

		// begin timing:
		auto t_start = clock();

		// recompute distances and then write them to cache file.
		// kd tree built upon m3d's vertices
		KNNTree* tree = new KNNTree();
		TriPoint p;
		for (unsigned i = 0; i < _m3d->vertices.size(); ++i)
		{
			p = _m3d->vertices[i];
			tree->insert( Point_and_uint(P3(p[0], p[1], p[2]), i) );
		}

		// closest points for each face of MA mesh on orig 3d mesh
		vector<int> closest(m_origG->faces.size()*2); 
		for (unsigned fi = 0; fi < m_origG->faces.size(); ++fi)
		{
			if ( _isects_per_face[fi*2 + 0][0] == numeric_limits<float>::max() )
			{
				// invalid closest
				closest[fi*2 + 0] = -1;
			}
			else
			{
				p = _isects_per_face[fi*2 + 0];
				K_neighbor_search search1(*tree, P3(p[0], p[1], p[2]), 1);
				K_neighbor_search::iterator nb_it = search1.begin();
				unsigned vi = boost::get<1>(nb_it->first);
				closest[fi*2 + 0] = (int)vi;

				p = _isects_per_face[fi*2 + 1];
				K_neighbor_search search2(*tree, P3(p[0], p[1], p[2]), 1);
				nb_it = search2.begin();
				vi = boost::get<1>(nb_it->first);
				closest[fi*2 + 1] = (int)vi;
			}
		}
		delete tree;

		// remove duplicate closest pair
		typedef TriEdge index_pair;
		struct index_pair_comp
		{
			bool operator () (const index_pair& _p1, const index_pair& _p2) const
			{
				return _p1[0] == _p2[0] ? 
					_p1[1] < _p2[1]
				: _p1[0] < _p2[0];
			}
		};
		set <index_pair, index_pair_comp> index_pair_set;
		for (unsigned fi = 0; fi < m_origG->faces.size(); ++fi)
		{
			int vi = closest[fi*2 + 0];
			int vj = closest[fi*2 + 1];
			if (vi == -1)
				continue;

			index_pair_set.insert(util::makeEdge(vi, vj));
		}
		cout << "unique pairs need to compute geodesic for: " << index_pair_set.size() << endl; // debug

		// build a graph (neighbor relationship) using closest pairs
		vector< vector<unsigned> > nbrs(_m3d->vertices.size());
		for (auto pr_it = index_pair_set.begin(); pr_it != index_pair_set.end(); ++pr_it)
		{
			nbrs[(*pr_it)[0]].push_back((*pr_it)[1]);
			//nbrs[(*pr_it)[1]].push_back((*pr_it)[0]);
		}
		index_pair_set.clear();

		// randomize the vertices sequence
		vector<unsigned> vts_seq;
		for (unsigned i = 0; i < nbrs.size(); ++i)
		{
			if (nbrs[i].empty())
				continue;
			vts_seq.push_back(i);
		}
		std::random_shuffle(vts_seq.begin(), vts_seq.end());
		cout << "# propagations to do: " << vts_seq.size() << endl; // debug

		// compute geodesic for each vert in the vts seq.
		cout << "computing geodesics..." <<endl;
		geodesic::Mesh geo_m;
		vector<double> pts;
		vector<unsigned> faces;
		for (unsigned i = 0; i < _m3d->vertices.size(); ++i)
		{
			pts.push_back(_m3d->vertices[i][0]);
			pts.push_back(_m3d->vertices[i][1]);
			pts.push_back(_m3d->vertices[i][2]);
		}
		for (unsigned i = 0; i < _m3d->faces.size(); ++i)
		{
			faces.push_back(_m3d->faces[i][0]);
			faces.push_back(_m3d->faces[i][1]);
			faces.push_back(_m3d->faces[i][2]);
		}

		geo_m.initialize_mesh_data(pts, faces);
		geodesic::GeodesicAlgorithmExact geo_algo(&geo_m);

		vector<geodesic::SurfacePoint> stop_pts;
		vector<geodesic::SurfacePoint> sources(1);
		map<index_pair, double> index_geod_map;

		// estimation of remaining time
		long t_interval_start, t_interval;
		bool do_timing = true;
		int interval_size = 100;
#pragma omp parallel default(shared) private(sources, stop_pts, geo_m)
		{
#pragma omp parallel for 
			for (int i = 0; i < vts_seq.size(); ++i)
			{
				if (do_timing)
				{
					t_interval_start = clock();
					do_timing = false;
				}

				int si = vts_seq[i];
				sources[0] = &geo_m.vertices()[si];
				// prepare target points 
				stop_pts.clear();
				for (int ni = 0; ni < nbrs[si].size(); ++ni)
				{
					unsigned ti = nbrs[si][ni];
					stop_pts.push_back(&geo_m.vertices()[ti]);
				}

				geo_algo.propagate(sources, geodesic::GEODESIC_INF, &stop_pts);	

				// get geodesic distance for each target
				double geo_dist;
				for (int ii = 0; ii < nbrs[si].size(); ++ii)
				{
					unsigned ti = nbrs[si][ii];
					geo_algo.best_source(geodesic::SurfacePoint(&geo_m.vertices()[ti]), geo_dist);
					index_geod_map[util::makeEdge(si, ti)] = geo_dist;
				}

				if (i % interval_size == 0 || i == vts_seq.size() - 1)
				{
					// reached the end of an interval. stop timer. reset flags.
					do_timing = true;
					// calculate time used for this interval.
					cout << "# of propagations done: " << i;
					t_interval = clock() - t_interval_start;
					auto t_used = t_interval * 1000.0f / CLOCKS_PER_SEC;
					if (i < vts_seq.size() -1)
					{
						// estimate how much time in total will be used
						auto t_total = t_used / (float)interval_size * vts_seq.size();
						t_total = t_total * 1000.0f / CLOCKS_PER_SEC;
						cout << " estimated total time: " << t_total <<"ms"<< endl;
					}
					cout << endl;
				}
			}
		}

		cout << "geodesics done!" << endl;

		// assign geodesic of each pair to corresponding face
		_geoD_per_face.resize( m_origG->faces.size(), -1.0f );//by default invalid value
		for (unsigned fi = 0; fi < m_origG->faces.size(); ++fi)
		{
			if (closest[fi*2+0] == -1)
				continue;
			index_pair pr = util::makeEdge(closest[fi*2+0], closest[fi*2+1]);
			_geoD_per_face[fi] = (float)(index_geod_map[pr]);
			//cout << "geodesics for " <<fi<<": "<<_geoD_per_face[fi]<< endl;
		}

		// end timing:
		auto t_duration = clock() - t_start;
		cout << "MGF computation finished within " << t_duration * 1000.0f / CLOCKS_PER_SEC << "ms."<<endl;

		// after recomputation, save results to cache file
		ofstream out_file;
		out_file.open(_cache_file);
		if (!out_file.is_open())
		{
			cout << "Error: couldn't open " << _cache_file << endl;
			return;
		}
		cout << "writing geo.d. to cache file "<<_cache_file<<endl;
		out_file.exceptions(std::ofstream::failbit | std::ofstream::badbit);
		try
		{
			for (unsigned fi = 0; fi < _geoD_per_face.size(); ++fi)
			{
				out_file << _geoD_per_face[fi] << endl;
			}
		}
		catch (std::ofstream::failure e)
		{
			std::cerr << "Error occurred while reading file "<<_cache_file<<endl;
		}

		cout << "writing geo.d. done!"<<endl;
		out_file.close();
	} // else
}

void SteinerGraph::computeDiffDist()
{
	if (bt3MA_face.empty())
		return;

	bt2MA_face.clear();
	bt2MA_face.resize(bt3MA_face.size());
	bt2bt3diffMA_face.clear();
	bt2bt3diffMA_face.resize(bt3MA_face.size());

	// compute burn dist for face by interpolating from its vertices' burn dist
	min_diffDistOrigFaces = numeric_limits<float>::max();
	max_diffDistOrigFaces = numeric_limits<float>::min();
	TriFace f;
	float burnD_f;
	float burn_v[3];
	// intermediate variables
	float s; 
	for (unsigned fi = 0; fi < m_origG->faces.size(); ++fi)
	{
		f = m_origG->faces[fi];
		// skip this face if at least one vert is not burnt
		if (!burnt[f[0]] && !burnt[f[1]] && !burnt[f[2]])
			continue;

		// obtain burn time of this sheet for each vert. 
		for (unsigned i = 0; i < 3; ++i)
		{
			unsigned vi = f[i];
			TriEdge e = util::oppositeEdge(f, vi);
			TopoGraph::EdgeIdx tei;
			int ret_code = mapTopo(vi, fi);
			assert(ret_code >= 0);
			tei = (TopoGraph::EdgeIdx)ret_code;
			burn_v[i] = bt2MA_vert_per_sheet[vi][tei];
		}

		burnD_f = (burn_v[0] + burn_v[1] + burn_v[2]) / 3.0f;
		bt2MA_face[fi] = burnD_f;

		// compute the difference between burn dist and dist2surf for this face
		bt2bt3diffMA_face[fi] = abs(bt3MA_face[fi] - burnD_f);

		// update stats
		min_diffDistOrigFaces = std::min(min_diffDistOrigFaces, bt2bt3diffMA_face[fi]);
		max_diffDistOrigFaces = std::max(max_diffDistOrigFaces, bt2bt3diffMA_face[fi]);
	}

	//cout << "min diff dist for orig faces: " << min_diffDistOrigFaces << endl; //debug
	//cout << "max diff dist for orig faces: " << max_diffDistOrigFaces << endl; //debug
}

void SteinerGraph::filterFaceByDiffDist(float _t, vector<unsigned>& _faces)
{
	_faces.clear();
	TriFace f;
	for (unsigned fi = 0; fi < m_origG->faces.size(); ++fi)
	{
		f = m_origG->faces[fi];
		if ( !burnt[f[0]] && !burnt[f[1]] && !burnt[f[2]] )
			continue;

		if (bt2bt3diffMA_face[fi] >= _t)
			_faces.push_back(fi);
	}
}

void SteinerGraph::dualizeFacesWithBurnPaths( 
	vector<std::pair<TriEdge, int> >& _burntE_assocF_list,
	DualOpt _edge_dual_opt,
	DualOpt _poly_dual_opt,
	float _eps)
{
	this->dual_vts.clear();
	this->dual_edges.clear();
	this->m_dualEdge_from_which_face.clear();
	this->bt2_MC.clear();
	this->m_bt2_MC_perFace.clear();
	this->m_is_face_dual.clear();
	this->m_finerTri_byDual.clear();
	this->m_origTri_for_finerTri_byDual.clear();
	//this->m_eDualFaceIdPr_polyDual_map.clear();
	this->m_fromFineTri_for_dualE.clear();

	struct TriEdgeLessThan
	{
		bool operator () (const TriEdge& _e1, const TriEdge& _e2) const
		{
			return _e1[0] == _e2[0] ? 
				_e1[1] < _e2[1] : _e1[0] < _e2[0];
		}
	};

	map <TriEdge, vector<int> > burntEdge_assocFace_map;
	for (auto pair_it = _burntE_assocF_list.begin(); 
		pair_it != _burntE_assocF_list.end(); ++pair_it)
	{
		TriEdge e = pair_it->first;
		burntEdge_assocFace_map[e].push_back(pair_it->second);
	}

	// for each face build a directed graph from steiner edges 
	// & the burnt edges within that face. 
	// each steiner edges are traversed as one half edge in a consistent winding direction;
	// while each burnt edge is traversed as two half edges
	bool reversed[] = {false, false, false};
	TriEdge edges[3];
	TriEdge e;
	TriFace f;
	// a list of such order info is stored for each vertex of the directed graph
	// it corresponds to the vertex pointed to by a directed edge:
	// <actual index of the vertex on the vts2traverse list, the order used for sorting>
	typedef pair<unsigned, unsigned> order_sortOrder_pair;
	typedef vector< list<order_sortOrder_pair> > vOrder_outVOrder_map; // vert order -> a list of out going vts orders
	vOrder_outVOrder_map directed_neighbors;
	// less-than comparator for a pair of <vid, its relative order on face boundary>
	struct v_order_pair_lessthan  
	{
		bool operator () (const order_sortOrder_pair& _x, const order_sortOrder_pair& _y)
		{
			return _x.second > _y.second;
		}
	};

	vector<unsigned> stVts;//vts on cur edge
	vector<int> stVts_other; // vts on other edge(s)
	vector<unsigned> vts2traverse; // traverse vts in this order later
	int nVts_eachEdge[3]; // # of st vts on each edge
	list<TriEdge> st_edges_in_face; // path edges within cur face
	vector<int> vOrder_eID_map; // st vert -> edge id
	//bool is_nonManifold_edge[3];

	vector< vector<unsigned> > poly_faces; // poly faces from planar subdivision of the cur. triangle face
	map<TriEdge, movingSum_edgeDual_struct> edge_dualV_map; // triedge -> the index of its dual vert
	vector< unsigned > dualVts_curPolyFace; 
	vector<int> assoc_faces; // one un-directed edge could associate with multi faces

	auto cout_uint_vector = [] (const char* _msg, const vector<unsigned>& _v, bool _condition) 
	{
		if (_condition)
		{
			cout << _msg;
			std::for_each( _v.begin(), _v.end(), [](unsigned _x){cout << _x<<" ";} );
			cout <<endl;
		}
	};
	typedef vOrder_outVOrder_map::const_iterator const_iter;
	auto cout_gNode = [&directed_neighbors] (const char* _msg, const_iter _it, bool _condition)
	{
		if (_condition)
		{
			cout << _msg;
			cout << _it - directed_neighbors.begin() <<"->";
			for (auto iit = _it->begin(); iit != _it->end(); ++iit)
			{
				cout << iit->first<<","<<iit->second<<" ";
			}
			cout << endl;
		}
	};

	set<unsigned> f_debug;
	//f_debug.insert(7276);
	unsigned bad_faces_cnt = 0; //debug
	unsigned cnt_fixed_faces = 0;

	vector<trimesh::vec> bndry_dirs;
	map<unsigned, unsigned> vIdx_burnDirIdx_map;
	findBoundaryBurnDir(bndry_dirs, vIdx_burnDirIdx_map);

	map<unsigned, std::pair<trimesh::vec, float>> vert_burnDirDistPair_map;

	// debug
	//auto os = ofstream("dual_struct.txt");
	//os << "# there are group(s) in this file, each group contains coordinates of a tri face f\n"
	//	<<"# and the polygon subdivision of this face.\n"
	//	<<"#\n"
	//	<<"# group :=\n"
	//	<<"# 	coordinates of tri face (always 9 numbers)\n"
	//	<<"# 	a number n, i.e. the number of polygons in the following poly subdivision of the face\n"
	//	<<"# 	n polygons' info\n"
	//	<<"# \n"
	//	<<"# polygon :=\n"
	//	<<"# 	coordinates of polygon's vertices (always integer multiples of 3)\n"
	//	<<"# 	dir (3 numbers) and distance of each vertex\n"
	//	<<"#	dual point of this polygon (3 numbers) and the distance of the dual point"
	//	<<endl;
	//os << "{"; // the start of the whole big list

	// if "steiner_only" burn scheme is used, put all orig. vts into dual list
	process_orig_vts_as_duals();

	// where do st vts on each edge start & end?
	pair<unsigned, unsigned> range_eachEdge[3];
	for (unsigned fi = 0; fi < m_origG->faces.size(); ++fi)
	{
		if (fi % 100000 == 0 /*1*/) // debug
			cout << "# of faces dualized: " << fi << endl;

		f = m_origG->faces[fi];
		directed_neighbors.clear();
		st_edges_in_face.clear();
		vts2traverse.clear();

		//if (st_vts[f[0]]==st_vts[f[1]] && st_vts[f[0]]==st_vts[f[2]]) // debug
		//{
		//	cout << "a point triangle: " << fi << endl;
		//	exit(1);
		//}

		// skip degenerate or not-to-dualize triangles
		/*if ( !util::valid_triangle(st_vts[f[0]], st_vts[f[1]], st_vts[f[2]]) ||
		( dualize && dualize->size() == orig_g.faces.size() && (*dualize)[fi] == false ) )
		continue;*/

		// whether we need to reverse the steiner points order
		// for each edge
		for ( unsigned ei = 0; ei < 3; ei++ )
		{
			e = util::makeEdge( f[ ei ], f[ ( ei + 1 ) % 3 ] );
			edges[ ei ] = e;
			if ( e[ 0 ] != f[ ei ] )
				reversed[ ei ] = true;
		}

		/*
		// check non-manifold-ness of each edge
		for (unsigned ei = 0; ei < 3; ei ++)
		{
		e = edges[ei];
		auto* nb_faces = &( this->orig_g.nbFacesOfEdge[e] );
		is_nonManifold_edge[ei] = nb_faces->size() > 2;
		}
		*/

		// save the traversal order as we traverse each st&orig vert
		for ( unsigned ei = 0; ei < 3; ei++ )
		{
			e = edges[ ei ];
			stVts.clear();
			m_stSubdiv.getStVertIndicesOnTriEdge( e, stVts );
			nVts_eachEdge[ ei ] = stVts.size();

			if ( reversed[ ei ] )
			{
				std::reverse( stVts.begin(), stVts.end() );
				stVts.insert( stVts.begin(), e[ 1 ] );
			}
			else
			{
				stVts.insert( stVts.begin(), e[ 0 ] );
			}
			std::copy( stVts.begin(), stVts.end(), back_inserter( vts2traverse ) );

		}

		// check if this face is part of "pocket" or not
		bool f_part_of_pocket = false;
		for ( auto v_it = vts2traverse.begin(); v_it != vts2traverse.end(); ++v_it )
		{
			// only care about whether burnable verts are burned or not on this face
			if ( getBurnScheme() == STEINER_ONLY && !m_stSubdiv.isSteinerVert( *v_it ) )
				continue;
			int ret_code = mapTopo( *v_it, fi );
			assert( ret_code >= 0 );
			auto tei = ( TopoGraph::EdgeIdx )ret_code;
			float burn_dist = this->bt2MA_vert_per_sheet[ *v_it ][ tei ];
			f_part_of_pocket = burn_dist == infiniteBurnDist();
			if ( f_part_of_pocket )
				break;
		}

		directed_neighbors.resize( vts2traverse.size() );

		// traverse each edge along the "face orientation": f[0->1->2]
		// [start_st_v, ..., end_st_v)
		// save this traversal order in vts2traverse list.
		// 
		// first, add to the directed graph using steiner edges
		// on each face's orig edges
		for ( unsigned i = 0; i < vts2traverse.size() - 1; ++i )
		{
			directed_neighbors[ i ].push_back( make_pair( i + 1, i + 1 ) );
		}
		directed_neighbors[ vts2traverse.size() - 1 ].push_back( make_pair( 0, vts2traverse.size() ) );

		// second, add to the directed graph using steiner edges
		// lying within each face, only if cur face is not part of pocket
		if ( !f_part_of_pocket )
		{
			unsigned start_v = 1;
			unsigned end_v = 1 + nVts_eachEdge[ 0 ];
			unsigned next_start_v, next_end_v;
			// where do st vts on each edge start & end?
			range_eachEdge[ 0 ].first = start_v;
			range_eachEdge[ 0 ].second = end_v;
			range_eachEdge[ 1 ].first = range_eachEdge[ 0 ].second + 1;
			range_eachEdge[ 1 ].second = range_eachEdge[ 0 ].second + 1 + nVts_eachEdge[ 1 ];
			range_eachEdge[ 2 ].first = ( range_eachEdge[ 1 ].second + 1 ) % vts2traverse.size();
			range_eachEdge[ 2 ].second = ( range_eachEdge[ 1 ].second + 1 + nVts_eachEdge[ 2 ] ) % vts2traverse.size();
			// now form st edges in face
			// 1. vts on e1 connected to those on e2 and e3 and oppo_v
			for ( unsigned ei = 0; ei < 3; ei++ )
			{
				// update the start / end index into vts2traverse for next edge
				const auto& s_e_pr = range_eachEdge[ ei ];
				const auto& s_e_pr_next = range_eachEdge[ ( ei + 1 ) % 3 ];
				start_v = s_e_pr.first;
				end_v = s_e_pr.second;
				next_start_v = s_e_pr_next.first;
				next_end_v = s_e_pr_next.second;
				unsigned oppo_v_order = next_end_v;

				for ( unsigned order = start_v; order != end_v; order = ( order + 1 ) % vts2traverse.size() )
				{
					for ( unsigned other_order = next_start_v;
						other_order != next_end_v; other_order = ( other_order + 1 ) % vts2traverse.size() )
					{
						e = TriEdge( vts2traverse[ order ], vts2traverse[ other_order ] );
						if ( getBurnScheme() == STEINER_ONLY &&
							( !m_stSubdiv.isSteinerVert( e[ 0 ] ) || !m_stSubdiv.isSteinerVert( e[ 1 ] ) )
							)
							continue;
						if ( burntEdge_assocFace_map.count( e ) > 0 ||
							burntEdge_assocFace_map.count( TriEdge( e[ 1 ], e[ 0 ] ) ) > 0 )
						{
							st_edges_in_face.push_back( util::makeEdge( order, other_order ) );
						}
					}
					if ( getBurnScheme() == ORIGINAL_AND_STEINER )
					{
						e = TriEdge( vts2traverse[ order ], vts2traverse[ oppo_v_order ] );
						if ( burntEdge_assocFace_map.count( e ) > 0 ||
							burntEdge_assocFace_map.count( TriEdge( e[ 1 ], e[ 0 ] ) ) > 0 )
						{
							st_edges_in_face.push_back( util::makeEdge( order, oppo_v_order ) );
						}
					}
				}
			}

			/*if (!st_edges_in_face.empty())
			cout << "in-face burnt st edges found in face " << fi << endl;*/

			// add these new paths to directed graph
			for ( auto path_it = st_edges_in_face.begin();
				path_it != st_edges_in_face.end(); ++path_it )
			{
				int order = ( *path_it )[ 0 ];
				int order_other = ( *path_it )[ 1 ];
				directed_neighbors[ order ].push_back( make_pair( order_other,
					order_other < order ? order_other + vts2traverse.size() : order_other ) );
				////printf("%d->%d, cosine: %f\n", stVts[j], stVts_other[i], cosine); //debug
				directed_neighbors[ order_other ].push_back( make_pair( order,
					order < order_other ? order + vts2traverse.size() : order ) );
				////printf("%d->%d, cosine: %f\n", stVts_other[i], stVts[j], cosine); //debug
			}
		}

		// now directed graph is ready. sort outgoing edges.
		for (auto outVOrder_list_it = directed_neighbors.begin(); 
			outVOrder_list_it != directed_neighbors.end(); ++outVOrder_list_it)
		{
			if (outVOrder_list_it->size() < 2)
				continue;
			outVOrder_list_it->sort(v_order_pair_lessthan());
		}

		// debug
		if (f_debug.count(fi))
		{
			std::cout << "----dualizing face " <<fi <<"----: "<<endl
				<</*m_stSubdiv.getVert*/(f[0])<<","
				<</*m_stSubdiv.getVert*/(f[1])<<","
				<</*m_stSubdiv.getVert*/(f[2])<< std::endl;

			std::cout << "traverse order: ";
			std::for_each( vts2traverse.begin(), vts2traverse.end(), 
				[this] (unsigned _vi) { 
					cout << _vi << "("<<topo_type_str[this->getTopoType(_vi)]<<") "; 
			} );
			std::cout << std::endl;

			std::cout << "in face graph edges: ";
			std::for_each( st_edges_in_face.begin(), st_edges_in_face.end(),
				[] (const TriEdge& _e) { cout << _e<<" "; });
			std::cout << std::endl;
		}

		// trace polygonal faces by following the directed graph
		trimesh::vec v;
		vector<unsigned> cur_polyFace;
		poly_faces.clear();
		bool polyFaceTrace_failed = false;
		for (unsigned order = 0; order < vts2traverse.size(); ++order)
		{
			cur_polyFace.clear();
			// init the first vert to trace from
			auto outgoing_list = directed_neighbors.begin()+order;
			cur_polyFace.push_back( vts2traverse[order] );

			// goto next vert. delete the neighbor once it is traced
			if (!outgoing_list->empty())
			{
				auto next_it = outgoing_list->begin();
				unsigned next_order = next_it->first;
				unsigned cur_order = order;
				auto outgoing_list_next = directed_neighbors.begin()+next_order;
				outgoing_list->erase(next_it);

				// loop until terminal node is reached or looping back
				unsigned iter_limit = 200, cnt = 0;//debug
				while( !outgoing_list_next->empty() && vts2traverse[next_order] != cur_polyFace[0] )
				{
					//if (cnt == iter_limit) //debug
					//{
					//	cout << "exiting on face " << fi << endl;
					//	exit(1);

					//	polyFaceTrace_failed = true;
					//	bad_faces_cnt++;
					//	goto NEXT_FACE;
					//}
					//cnt++;//debug

					// the relative order used when comparing
					unsigned relative_order = cur_order < next_order ? cur_order + vts2traverse.size() : cur_order;

					// find the neighbor (of next_order) to trace to
					for (auto nb_it = outgoing_list_next->begin(); nb_it != outgoing_list_next->end(); ++nb_it)
					{
						if (nb_it->first == cur_order || nb_it->second > relative_order)
							continue;
						outgoing_list = outgoing_list_next;
						cur_polyFace.push_back( vts2traverse[next_order] );
						outgoing_list_next = directed_neighbors.begin()+nb_it->first;
						cur_order = next_order;
						next_order = nb_it->first;
						outgoing_list->erase(nb_it);

						if ( f_debug.count( fi ) )
						{
							cout_uint_vector( "cur_polyFace:", cur_polyFace, true );
							cout_gNode( "cur outgoing list: ", outgoing_list, true );
							cout_gNode( "next outgoing list: ", outgoing_list_next, true );
							std::cout << endl;
						}

						break;
					}
				}

				if (vts2traverse[next_order] == cur_polyFace[0])
				{
					poly_faces.push_back(cur_polyFace);

					// debug
					/*cout << "cur poly face: ";
					std::for_each(cur_polyFace.begin(), cur_polyFace.end(), [] (unsigned _vi) { cout << _vi << " "; });
					cout << endl;*/
				}

			}
		}// tracing poly faces done.

		if ( f_debug.count( fi ) )
		{
			std::cout << "finished tracing poly faces." << endl;
		}

		/* dualize poly face */

		unsigned start_newDualEdges = dual_edges.size();
		vert_burnDirDistPair_map.clear();
		if (!f_part_of_pocket )
		{
			find_burn_dir_dist_of_face( vts2traverse, fi,
				bndry_dirs, vIdx_burnDirIdx_map, vert_burnDirDistPair_map,
				f_debug.count( fi ) );
		}

		if ( f_debug.count( fi ) )
		{
			std::cout << "finished find_burn_dir_dist_of_face()." << endl;
		}

		dualize_poly_subdivision(vts2traverse, fi, f_part_of_pocket,
			poly_faces,
			vert_burnDirDistPair_map, 
			burntEdge_assocFace_map, edge_dualV_map,
			_edge_dual_opt, _poly_dual_opt, 0.080001f, _eps, 
			/*false*/f_debug.count(fi) > 0, /*&os*/nullptr);

		if ( f_debug.count( fi ) )
		{
			std::cout << "finished dualize_poly_subdivision()." << endl;
		}

		/*if (fi != m_origG->faces.size()-1)
		os << ",";
		os << endl;*/

		//if (this->dual_vts.size()-1 == 69716) // debug
		//	cout << "dual v "<< 69716 << " is from face "<<fi<<endl;

		//if (_dual_method == SIMPLE_DUAL || _dual_method == WEIGHT_DUAL)
		//{
		//	dualize_poly_subdivision_simple(fi, poly_faces, 
		//		vert_burnDirDistPair_map, burntEdge_assocFace_map, edge_dualV_map, 
		//		_dual_method,
		//		f_debug.count(fi) > 0);
		//}
		//else 
		//{
		//	dualize_poly_subdivision(vts2traverse, fi, poly_faces,
		//		vert_burnDirDistPair_map, 
		//		burntEdge_assocFace_map, edge_dualV_map, stV_triEdge_map,
		//		_dual_method, 0.080001f,
		//		false/*f_debug.count(fi) > 0*/, &os);
		//	if (fi != orig_g.faces.size()-1)
		//		os << ",";
		//	os << endl;
		//}
		// record cur face as "from-face" for the just added dual edges
		for (unsigned j = 0; j < dual_edges.size() - start_newDualEdges; ++j)
			m_dualEdge_from_which_face.push_back(fi);

NEXT_FACE:
		;
	}
	//os << "}"; // the end of the whole big list
	//os.close();

	// don't forget to finalize dual vertices' coordinates values and distance
	finalize_dual_vts(edge_dualV_map);
	cout << "after dual vts finalization, # dual vts: "<<dual_vts.size() << endl;

	// clear up 
	edge_dualV_map.clear();
	burntEdge_assocFace_map.clear();

	cout << "# faces dualization failed on: " << bad_faces_cnt << endl; // debug

	/// compute topology numbers for dual structure
	cout << "checking topological properties of dual structure ..." << endl;

	/// check if there are st vts that cause any disconnection in the dual structure
	// they are characterized by the following properties:
	// - topo graphs have more than 1 conn. component
	unsigned cnt_isolateVts = 0;
	for (unsigned vi = 0; vi < this->m_stSubdiv.sizeOfVts(); ++vi)
	{
		if (simpleTopoAt(vi) || !TopoGraph::has2dNbhood(vi))
			continue;

		const auto &tg = this->t_graphs[vi];
		if (tg->nConnCmpnt() > 1)
			cnt_isolateVts++;
	}

	cout << "# isolation st vts (causing disconnection in dual structure): " << cnt_isolateVts << endl; 
	cout << endl;
}

bool SteinerGraph::isDualValid() const
{
	return m_dual_valid;
}

bool SteinerGraph::isFaceDual(int _dual_id) const
{
	return m_is_face_dual[_dual_id] >= 0;
}

int SteinerGraph::isEdgeDual( int _dual_id ) const
{
	int ret_id = m_dualV_triEdgeIdx[ _dual_id ];
	return ret_id >= 0 ? ret_id : -1;
}

const vector<TriPoint>& SteinerGraph::getDualVts() const
{
	return dual_vts;
}

vector<int> SteinerGraph::checkUnburnt()
{
	map<int, bool> vert_visited;
	map<int, bool> face_visited;
	for (unsigned i = 0; i < m_origG->vts.size(); ++i)
		if (!burnt[i])
			vert_visited[i] = false;
	//set<TriEdge, EdgeCompare> component;
	vector<TriEdge> component_edges;
	set<int> component_faces;
	map<TriEdge, vector<int>, EdgeCompare> nbFacesOfEdge;

	// ideally removing one bad face for each bad component will make it burnable.
	vector<int> bad_faces; 
	for (unsigned vi = 0; vi < this->m_origG->vts.size(); ++vi)
	{
		if (burnt[vi] || vert_visited[vi])
			continue;

		queue<int> q;
		q.push(vi);
		vert_visited[vi] = true;
		component_edges.clear();
		component_faces.clear();

		// use flooding to find an unburnt connected component
		while(!q.empty())
		{
			int cur_v = q.front();
			q.pop();
			vector<int>& nbs = m_origG->nbVtsOfVert[cur_v];
			vector<int>& nbFaces = m_origG->nbFacesOfVert[cur_v];

			if (cur_v == 82)
			{
				cout << cur_v <<"'s nbs: ";
				for (unsigned ni = 0; ni < nbs.size(); ++ni)
				{
					int nb = nbs[ni];
					cout << nb << ", ";
				}
				cout << endl;
			}

			for (unsigned ni = 0; ni < nbs.size(); ++ni)
			{
				int nb = nbs[ni];
				if (!burnt[nb])
					component_edges.push_back(TriEdge(cur_v, nb));
				if (!burnt[nb])
				{
					if (!vert_visited[nb])
					{
						q.push(nb);
						vert_visited[nb] = true;
					}
				}
			}

			for (auto fi = nbFaces.begin(); fi != nbFaces.end(); ++fi)
			{
				TriFace f = m_origG->faces[*fi];
				if (!burnt[f[0]] + !burnt[f[1]] + !burnt[f[2]] == 3)
				{
					component_faces.insert(*fi);
				}
			}
		}

		// now we have a component. print it?
		// TODO: do something smart to see if it's a "pocket"
		if (!component_edges.empty())
		{
#ifdef PRINT_UNBURNT
			cout << "-------------An unburnt component-------------: " << endl;
			for (auto f_it = component_faces.begin(); f_it != component_faces.end(); ++f_it)
			{
				TriFace f = m_origG->faces[*f_it];
				cout << "<" <<f[0]<<", "<<f[1]<<", "<<f[2] <<">\t";
			}
			cout << endl;
#endif // PRINT_UNBURNT

			// group faces according to shared edge
			// and set each face unvisited
			nbFacesOfEdge.clear();
			for (auto f_it = component_faces.begin(); f_it != component_faces.end(); ++f_it)
			{
				vector<TriEdge> edges = util::edgesFromFace(m_origG->faces[*f_it]);
				for (auto e_it = edges.begin(); e_it != edges.end(); ++e_it)
				{
					nbFacesOfEdge[*e_it].push_back(*f_it);
				}

				face_visited[*f_it] = false;
			}

			// now perform open operation on current component
			for (auto f_it = component_faces.begin(); f_it != component_faces.end(); ++f_it)
			{
				if (face_visited[*f_it])
					continue;

				int seed_fi = *f_it;// start from this face opening
				bad_faces.push_back(seed_fi);

				queue<int> face_q;
				face_q.push(seed_fi);
				vector<TriEdge> edges;
				vector<int> unvisited_faces;
				while (!face_q.empty())
				{
					int cur_fi = face_q.front();
					face_q.pop();
					if (face_visited[cur_fi])
						continue;

					TriFace f = m_origG->faces[cur_fi];
					face_visited[cur_fi] = true;

					// add neighbor face to q iff: 
					// the face is the only unvisited face sharing the edge
					edges = util::edgesFromFace(f);
					for (auto eit = edges.begin(); eit != edges.end(); ++eit)
					{
						vector<int>& nbFaces = nbFacesOfEdge[*eit];
						unvisited_faces.clear();

						// count unvisited nb faces
						for (auto nbit = nbFaces.begin(); nbit != nbFaces.end(); ++nbit)
						{
							auto visited = face_visited[*nbit];
							if (!visited)
								unvisited_faces.push_back(*nbit);
						}

						if (unvisited_faces.size() == 1)
							face_q.push(unvisited_faces.front());
					}
				}
			}
		}
	}

	return bad_faces;
}

bool SteinerGraph::isBurnPath(const TriEdge& e, int _fi, int& _from, int& _to)
{
	unsigned u = e[0];
	unsigned v = e[1];

	TopoGraph::EdgeIdx te_u = mapTopo(u, _fi);
	TopoGraph::EdgeIdx te_v = mapTopo(v, _fi);
	if (is_leading_topo_edge(u, te_u) && 
		prev_vert[u][te_u] == v && prev_tedge[u][te_u] == te_v) // v is pre to u
	{
		_from = v; _to = u;
		return true;
	}
	else if (is_leading_topo_edge(v, te_v) && 
		prev_vert[v][te_v] == u && prev_tedge[v][te_v] == te_u) // u is pre to v
	{
		_from = u; _to = v;
		return true;
	}
	else
		return false;
}

//////////////////////////////////////////////////////////////////////////
/// -----------definition of TopoGraph-----------
//////////////////////////////////////////////////////////////////////////
//TopoGraph::TopoGraph()
//{
//
//}

namespace util
{
	float rescale(float v, float vmin, float vmax, float _exp, float _new_min, float _new_max)
	{
		float dv;

		if (v < vmin)
			v = vmin;
		if (v > vmax)
			v = vmax;
		dv = vmax - vmin;

		float w = pow( (v - vmin) / dv, _exp );
		return (1.0f - w) * _new_min + w * _new_max;
	}

#define N_KEY_COLORS 9
	const TriColor jet_map[N_KEY_COLORS] = {
		TriColor ( 0.0f, 0.0f, 0.6f ), // dark blue
		TriColor ( 0.0f, 0.0f, 1.0f ), // blue
		TriColor ( 0.0f, 0.5f, 1.0f ), // azure
		TriColor ( 0.0f, 1.0f, 1.0f ), // cyan
		TriColor ( 0.5f, 1.0f, 0.5f ), // light green
		TriColor ( 1.0f, 1.0f, 0.0f ), // yellow
		TriColor ( 1.0f, 0.5f, 0.0f ), // orange
		TriColor ( 1.0f, 0.0f, 0.0f ), // red
		TriColor ( 0.5f, 0.0f, 0.0f ) // dark red
	};

	TriColor GetColour(float v,float vmin,float vmax)
	{
		TriColor c(1.0f,1.0f,1.0f); // white

		// get subrange v falls into
		float range, dv;

		if (v < vmin)
			v = vmin;
		if (v > vmax)
			v = vmax;
		range = std::max( vmax - vmin, 0.0000001f );
		dv = (1.0f / (N_KEY_COLORS - 1) );

		float multiples = ( (v - vmin) / range ) / dv;
		int min_idx = multiples;
		const auto& min_c = jet_map[min_idx];
		const auto& max_c = jet_map[min_idx + 1];
		c = trimesh::mix( min_c, max_c, (multiples - min_idx) );

		/*float dv;

		if (v < vmin)
		v = vmin;
		if (v > vmax)
		v = vmax;
		dv = vmax - vmin;

		if (v < (vmin + 0.25f * dv)) {
		c[0] = 0;
		c[1] = 4 * (v - vmin) / dv;
		} else if (v < (vmin + 0.5f * dv)) {
		c[0] = 0;
		c[2] = 1 + 4 * (vmin + 0.25f * dv - v) / dv;
		} else if (v < (vmin + 0.75f * dv)) {
		c[0] = 4 * (v - vmin - 0.5f * dv) / dv;
		c[2] = 0;
		} else {
		c[1] = 1 + 4 * (vmin + 0.75f * dv - v) / dv;
		c[2] = 0;
		}*/

		return(c);
	}

	void find_barycentric(const TriPoint& _p, 
		const TriPoint& _a, const TriPoint& _b, const TriPoint& _c,
		float& _s, float& _t, float& _w)
	{
		float area_abc = trimesh::len((_b - _a).cross(_c - _a));
		float area_pbc = trimesh::len((_b - _p).cross(_c - _p));
		float area_pab = trimesh::len((_a - _p).cross(_b - _p));

		_s = area_pbc / area_abc;
		_w = area_pab / area_abc;
		_t = 1.0f - _s - _w;
	}
}