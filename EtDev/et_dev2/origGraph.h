#ifndef ORIG_GRAPH_H
#define ORIG_GRAPH_H

#include <map>
#include <vector>

#include "allTypes.h"

using namespace std;

class MyGraph
{
	/// types
public:
	typedef vector<vector<int>> NbFaceForEdge;
	typedef map<TriEdge, float, EdgeCompare> WeightMap;
	typedef WeightMap::iterator WeightMapItor;

public:
	// construct graph from given vertices, edges and faces (connectivity information)
	MyGraph() {}
	MyGraph(const vector<TriPoint>& _vts, const vector<TriEdge>& _edges, const vector<TriFace>& _faces);
	~MyGraph() {}

	// subdivide the graph by subdividing each face n times 
	std::shared_ptr<MyGraph> subdivide(int _n);
	std::shared_ptr<MyGraph> subdivideOnce();

	bool isManifoldEdge(const TriEdge& _e) const;
	/*bool isManifoldEdge(const unsigned _ei) const;*/

	// return: 
	// - -1 if _f is not a valid face in graph;
	// - positive index number if _f is found 
	int getFaceIdx(const TriFace _f) const;

	// -1 if _e is not a valid edge in graph;
	int getEdgeIdx(const TriEdge& _e) const;

	// return list of neighboring faces' idx
	const vector<int>& getNbFaces( const TriEdge& _e ) const;
	const vector<int>& getNbFaces( int _ei ) const;

	// return the isolate edges in the graph (i.e. not used by any faces)
	vector<TriEdge> getIsolateEdges() const;

public:
	string m_name;
	vector<TriPoint> vts;
	vector<TriEdge> edges;
	vector<TriFace> faces;
	vector<vector<int>> nbVtsOfVert;
	vector<vector<int>> nbFacesOfVert;
	// map a given edge to its neighbor faces
	vector<vector<int>> nbFacesOfEdge;
	// map a given edge to its weight
	WeightMap weights;
	// triface -> its idx in faces list
	map<TriFace, unsigned, FaceCompare> triFace_faceIdx_map;
	// edge -> its idx in edge list
	map<TriEdge, unsigned> m_triEdge_idx_map;

private:
	// compute # of connected components of orig. graph
	// by flooding along neighbor faces
	void findConnCompntsInGraph(int& _n);
};

inline bool MyGraph::isManifoldEdge(const TriEdge& _e) const
{
	return getNbFaces(_e).size() == 2;
}

//inline bool MyGraph::isManifoldEdge(const unsigned _ei) const
//{
//	auto iter = nbFacesOfEdge.find(_e);
//	if (iter == nbFacesOfEdge.end())
//		return false;
//	return iter->second.size() == 2;
//}

inline int MyGraph::getFaceIdx(const TriFace _f) const
{
	auto find_it = this->triFace_faceIdx_map.find(_f);
	
	return find_it == triFace_faceIdx_map.end() ? -1 : find_it->second;
}

inline int MyGraph::getEdgeIdx(const TriEdge& _e) const
{
	auto find_it = this->m_triEdge_idx_map.find(_e);

	return find_it == m_triEdge_idx_map.end() ? -1 : find_it->second;
}

inline const vector<int>& MyGraph::getNbFaces( const TriEdge & _e ) const
{
	return nbFacesOfEdge[ getEdgeIdx( _e ) ];
}

inline const vector<int>& MyGraph::getNbFaces( int _ei ) const
{
	return nbFacesOfEdge[ _ei ];
}

inline vector<TriEdge> MyGraph::getIsolateEdges() const
{
	vector<TriEdge> res; 
	res.reserve( edges.size() );
	for ( auto i = 0; i < edges.size(); ++i )
		if ( getNbFaces( i ).size() == 0 )
			res.push_back( edges[ i ] );
	res.shrink_to_fit();
	return res;
}

#endif