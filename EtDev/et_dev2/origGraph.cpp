// A simple 2-simplicial complex graph
// 
// Copyright (C) 2018 Yajie Yan <danielyan86129@hotmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.

//////////////////////////////////////////////////////////////////////////
/// -----------definition of MyGraph-----------
//////////////////////////////////////////////////////////////////////////

#include <queue>
#include "origGraph.h"
#include "commonUtils.h"

/// ------------------ public non-static of MyGraph-----------
MyGraph::MyGraph( const vector<TriPoint>& _vts, const vector<TriEdge>& _edges, const vector<TriFace>& _faces )
{
	// store vertices
	this->vts = _vts;
	this->nbVtsOfVert.resize(vts.size());

	// store faces	
	this->faces = _faces;
	this->nbFacesOfVert.resize(vts.size());

	// store any pre-existing edges
	this->edges = _edges;
	for ( auto i = 0; i < edges.size(); ++i )
	{
		auto e = edges[ i ];
		m_triEdge_idx_map[ util::makeEdge( e[ 0 ], e[ 1 ] ) ] = i;
	}

	// extract edges from faces
	// add face to edge's nbface list
	// add face to vertex's nbface list
	vector<TriEdge> es;
	for (unsigned fi = 0; fi < this->faces.size(); ++fi)
	{
		TriFace& f = this->faces[fi];
		es = util::edgesFromFace(f);
		for (unsigned j = 0; j < es.size(); ++j)
		{
			int presize = m_triEdge_idx_map.size();
			auto& e_idx = this->m_triEdge_idx_map[ es[ j ] ];
			if ( m_triEdge_idx_map.size() != presize)
			{
				this->edges.push_back(es[j]);
				e_idx = this->edges.size() - 1;
			}
		}

		this->nbFacesOfVert[f[0]].push_back(fi);
		this->nbFacesOfVert[f[1]].push_back(fi);
		this->nbFacesOfVert[f[2]].push_back(fi);

		this->triFace_faceIdx_map[util::makeFace(f[0], f[1], f[2])] = fi;
	}

	// build edge-face-adjacency
	this->nbFacesOfEdge.resize( edges.size() );
	for ( unsigned fi = 0; fi < this->faces.size(); ++fi )
	{
		es = util::edgesFromFace( this->faces[ fi ] );
		for (auto j = 0; j < es.size(); ++j )
		{
			auto ei = m_triEdge_idx_map.find( es[ j ] )->second;
			nbFacesOfEdge[ ei ].push_back( fi );
		}
	}

	// build vert-vert-adjacency
	for (unsigned i = 0; i < this->edges.size(); ++i)
	{
		TriEdge& e = this->edges[i];
		this->nbVtsOfVert[e[0]].push_back(e[1]);
		this->nbVtsOfVert[e[1]].push_back(e[0]);
	}

	int n_conn_compnts = 0;
	findConnCompntsInGraph(n_conn_compnts);
	cout << "# connected components in orig. graph: " << n_conn_compnts << endl;
}

std::shared_ptr<MyGraph> MyGraph::subdivideOnce()
{
	return subdivide(1);
}

std::shared_ptr<MyGraph> MyGraph::subdivide(int _n)
{
	// vts, edges, faces of previous subdivision and curr subdivision
	vector<TriPoint>* pre_vts = new vector<TriPoint>(vts);
	vector<TriPoint>* new_vts = new vector<TriPoint>();
	vector<TriEdge>* pre_edges = new vector<TriEdge>(edges);
	vector<TriEdge>* new_edges = new vector<TriEdge>();
	vector<TriFace>* pre_faces = new vector<TriFace>(faces);
	vector<TriFace>* new_faces = new vector<TriFace>();
	// edge -> new vert id on it
	map<TriEdge, int> newVertOnEdge;

	for (int c = 0; c < _n; ++c)
	{
		newVertOnEdge.clear();
		new_vts->reserve(pre_vts->size() + pre_edges->size());
		new_faces->reserve(pre_faces->size() * 4);

		for (unsigned ei = 0; ei < pre_edges->size(); ++ei)
		{
			TriPoint v1 = vts[(*pre_edges)[ei][0]];
			TriPoint v2 = vts[(*pre_edges)[ei][1]];
			new_vts->push_back(trimesh::mix(v1, v2, 0.5f));
			newVertOnEdge[(*pre_edges)[ei]] = new_vts->size()-1;
		}

		// within each old face, connect new vts to form new edges and new faces
		TriFace new_f, f;
		vector<TriEdge> es;
		for (unsigned fi = 0; fi < pre_faces->size(); ++fi)
		{
			f = faces[fi];
			es = util::edgesFromFace(f);
			int v1 = f[0], v2 = f[1], v3 = f[2];
			int v12 = newVertOnEdge[es[0]];
			int v23 = newVertOnEdge[es[1]];
			int v13 = newVertOnEdge[es[2]];

			new_f = TriFace(v1, v12, v13);
			new_faces->push_back(new_f);
			new_f = TriFace(v12, v2, v23);
			new_faces->push_back(new_f);
			new_f = TriFace(v13, v23, v3);
			new_faces->push_back(new_f);
			new_f = TriFace(v23, v13, v12);
			new_faces->push_back(new_f);

			new_edges->insert(new_edges->end(), es.begin(), es.end());
		}

		// lastly, switch pre_ new_ pointer, so that curr elements
		// become pre_ in next subdivision 
		pre_vts->insert(pre_vts->end(), new_vts->begin(), new_vts->end());
		std::swap(pre_edges, new_edges);
		std::swap(pre_faces, new_faces);
		new_vts->clear();
		new_edges->clear();
		new_faces->clear();
	}// end of iteration of subdivision

	std::shared_ptr<MyGraph> g( new MyGraph(*pre_vts, *pre_edges, *pre_faces) );

	delete pre_vts, pre_edges, pre_faces;
	delete new_vts, new_edges, new_faces;

	return g;
}

void MyGraph::findConnCompntsInGraph(int& _n)
{
	_n = 0;
	vector<bool> visited_face(faces.size(), false);
	set<TriEdge> visited_edges;
	std::queue<unsigned> q;
	for (unsigned fi = 0; fi < faces.size(); ++fi)
	{
		if (visited_face[fi])
			continue;

		_n ++; // starting a new conn. compnt

		visited_face[fi] = true;
		q.push(fi);
		while (!q.empty())
		{
			unsigned cur_fi = q.front();
			q.pop();
			auto cur_f = faces[cur_fi];
			auto edges = util::edgesFromFace(cur_f);
			for (auto e_it = edges.begin(); e_it != edges.end(); ++e_it)
			{
				if (visited_edges.count(*e_it) > 0)
					continue;
				visited_edges.insert(*e_it);
				const auto& nb_faces = getNbFaces(*e_it);
				for (auto f_it = nb_faces.begin(); f_it != nb_faces.end(); ++f_it)
				{
					visited_face[*f_it] = true;
					if (*f_it != cur_fi)
						q.push(*f_it);
				}
			}
		}
	}
}
