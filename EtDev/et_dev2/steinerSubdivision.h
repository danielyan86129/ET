// Steiner subdivision on graph
// 
// Copyright (C) 2018 Yajie Yan <danielyan86129@hotmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef STEINER_SUBDIV_H
#define STEINER_SUBDIV_H

#include <vector>
#include <map>

#include "allTypes.h"

using namespace std;

class SteinerSubdivision
{
public:
	enum SubdivideScheme {
		FIXED,		// put m_nFixed steiner points on each edge
		ADAPTIVE,	// distance between 2 neighboring st vts not exceed d
		ADAPTIVE_MIDPOINT
	};

	enum ConnectScheme 
	{
		// st vts will connect to opposite orig vert within a triangle
		CONNECT_TO_ORIG_VTS, 
		// st vts will not connect to opposite orig vert within a triangle
		NOT_CONNECT_TO_ORIG_VTS 
	};

	/***** cons and de-cons. *****/
	SteinerSubdivision();
	~SteinerSubdivision();

	/***** public interfaces *****/
	
	/*return currently set connect scheme*/
	ConnectScheme CurrentConnectScheme() const;

	/* fixed subdivision */
	bool fixedSubdivide(const shared_ptr<MyGraph>& _g, int _nSamples);

	/* adaptive subdivision */
	bool adaptiveSubdivide(const shared_ptr<MyGraph>& _g_ptr, const double _param);

	/* adaptive (mid-point) subdivision */
	bool adaptiveMidPointSubdivide(const shared_ptr<MyGraph>& _g_ptr, const double _param);

	//
	// managing steiner/all vts
	//
	/* get/set the position of a vert with the given index (into the all-vts list) */
	TriPoint getVert(const unsigned _id) const;
	void setVert(const unsigned _id, TriPoint _pos);

	/* get/set specifically to a steiner vertex (id in [0, size-of-st-vts) */
	TriPoint getStVert(const unsigned _i) const;
	void setStVert(const unsigned _i, TriPoint _pos);

	/* return size of st vts */
	unsigned sizeOfStVts() const;

	/* Given an index into the all-vts list, does it correspond to a steiner vert (or orig vert) */
	bool isSteinerVert(const unsigned _i) const;

	/* return the tri edge that the given st vert resides on */
	TriEdge getResidingEdge(const unsigned _vi) const;
	unsigned getResidingEdgeIdx(const unsigned _vi) const;
	void getResidingEdgeIdx(const unsigned _vi, unsigned& _ei, int& _idx) const;

	/* return the residing face if given a st edge, i.e. both ends are st vert
	   return -1 if _e is not a valid st edge
	   */
	int getResidingFaceIdx(const TriEdge _e) const;

	/* return the range of steiner vts within within all vts stored */
	unsigned getStartOfStVts() const;
	unsigned getEndOfStVts() const;

	/* get/set st edge by index */
	TriEdge getStEdge(const unsigned _ei) const;
	TriEdge getOrigTriEdge(const unsigned _ei) const;

	/* get the index of the given steiner edge.
	   -1 returned if not found 
	   */
	int getStEdgeIndex(const TriEdge _st_e) const;

	/* is this an orig tri edge? */
	bool isTriEdge(const unsigned _ei) const;
	/* whether the query tri edge has any st vts on it */
	bool triEdgeSubdivided(const unsigned _i) const;

	/*
	append the steiner vts for the query edge to _stVts
	*/
	void getStVertIndicesOnTriEdge(
		const unsigned _ei, vector<unsigned>& _stVts) const;
	// TODO: OPTIMIZE_AWAY_THIS
	void getStVertIndicesOnTriEdge(
		const TriEdge& _e, vector<unsigned>& _stVts) const;

	/* get all vts (orig & steiner) on the given edge */
	void getVtsOnTriEdge(const unsigned _ei, vector<unsigned>& _vts) const;
	
	/*
	the position of the st vts on the query edge
	*/
	void getStVtsPosOnTriEdge(
		const unsigned _ei, vector<TriPoint>& _stVts) const;
	void getStVtsPosOnTriEdge(
		const TriEdge& _e, vector<TriPoint>& _stVts) const;

	/*
	given edge <u, v>, and one of its end, return the closest vert index
	Precondition: _which must be 0 or 1. 
	0 = u, 1 = v
	*/
	int getClosestVertOnEdge(unsigned _u, unsigned _v, unsigned _which) const;

	/*
	append to the output list the neighbor vts of _vi that are restricted within the given face
	@param _exclude_orig: when _vi is steiner, do not connect to orig vts.
	return code:
	- -1 if failed.
	*/
	int getNbVtsInFace(const unsigned _vi, const unsigned _fi, bool _exclude_orig, vector<int>& _nb_vts) const;

	//
	// insert st edges on an original triangle edge to the list
	//
	void getStEdgesOfEdge(const unsigned _ei, vector<TriEdge>& _st_edges_on_e) const;
	//
	// insert st edges within an original triangle face to the list
	//
	void getStEdgesWithinFace(const unsigned _fi, vector<TriEdge>& _st_edges_within_f) const;
	//
	// get the st edges on a face's interior and boundary 
	//
	void getStEdgesOfFace(const unsigned _fi, vector<TriEdge>& _st_edges_of_f) const;

	/*
	get the neighbor vertices 
	return code:
	- -1 if failed.
	*/
	void getNbVerts(const unsigned _vi, vector<unsigned>& _nbvts);

	/* return size of important members */
	unsigned sizeOfVts() const;
	unsigned sizeOfStEdges() const;

	/* getters to data member themselves */
	const vector<TriPoint>& getAllVts() const;
	// DEPRECATED
	const vector<TriEdge>& getStEdges() const;
	float getIntervalSize() const;

	/***** private helpers *****/	
private:
	void clear();

	/* create all st edges of a face */
	void connect_st_edges_of_face(const TriFace _f, set<TriEdge>& _visited_origedges, vector<TriEdge>& _edges);

	/* build mapping from a steiner edge to its index in the st edge list */
	void build_st_edge_mapping();

private:
	// the connection scheme
	ConnectScheme m_conn_scheme;
	// data vector for st. vts positions
	vector<TriPoint> m_stVts;
	// the start index and # consecutive st vts for each edge
	vector<pair<unsigned, unsigned>> m_triEdge_stVtsRange_list;
	// st vert id -> the residing tri edge idx & the idx within that edge
	// (-1 if 0 st vts on that edge)
	vector<pair<unsigned, int>> m_stV_triEdge_map;
	// total num of st edges
	size_t m_n_stEdges;
	// TODO: optimize away
	// list of all st. edges after subdivision
	vector<TriEdge> m_stEdges;
	// TODO: optimize away
	// st edge -> its index
	map<TriEdge, unsigned> m_stEdge_idx_map;
	// tri edge -> the st vts on it. OPTIMIZE_AWAY_THIS
	map<TriEdge, pair<unsigned, unsigned>> m_triEdge_stVts_map;
	
	// current subdiv scheme
	SubdivideScheme m_subd_scheme;
	// param for fixed scheme # sample on each edge
	unsigned m_nFixed;
	// param for adaptive scheme Neighboring st vts' distance not exceed this value.
	float m_d;

	std::shared_ptr<MyGraph> m_origG_ptr;
};

#endif

#include "steinerSubdivision_inline.h"