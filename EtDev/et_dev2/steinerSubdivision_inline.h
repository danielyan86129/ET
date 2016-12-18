#ifndef STEINER_SUBDIV_INLINE_H
#define STEINER_SUBDIV_INLINE_H

#include <algorithm>
#include <type_traits>

#include "steinerSubdivision.h"
#include "origGraph.h"

/* getters */
inline SteinerSubdivision::ConnectScheme SteinerSubdivision::CurrentConnectScheme() const
{
	return m_conn_scheme;
}

inline TriPoint SteinerSubdivision::getVert(const unsigned _i) const
{
	return m_stVts[_i];
}

inline void SteinerSubdivision::setVert(const unsigned _i, TriPoint _pos)
{
	m_stVts[_i] = _pos;
}

inline TriPoint SteinerSubdivision::getStVert(const unsigned _id) const
{
	return m_stVts[getStartOfStVts() + _id];
}

inline void SteinerSubdivision::setStVert(const unsigned _id, TriPoint _pos)
{
	m_stVts[_id] = _pos;
}

inline TriEdge SteinerSubdivision::getResidingEdge(const unsigned _i) const
{
	return m_origG_ptr->edges[m_stV_triEdge_map[_i - m_origG_ptr->vts.size()].first];
}

inline unsigned SteinerSubdivision::getResidingEdgeIdx(const unsigned _i) const
{
	return m_stV_triEdge_map[_i - m_origG_ptr->vts.size()].first;
}

inline void SteinerSubdivision::getResidingEdgeIdx(const unsigned _vi, unsigned& _ei, int& _idx) const
{
	const auto& ei_idx_pr = m_stV_triEdge_map[_vi - m_origG_ptr->vts.size()];
	_ei = ei_idx_pr.first;
	_idx = ei_idx_pr.second;
}

inline bool SteinerSubdivision::isSteinerVert(const unsigned _i) const
{
	return _i >= m_origG_ptr->vts.size();
}

inline unsigned SteinerSubdivision::getStartOfStVts() const
{
	return (unsigned)m_origG_ptr->vts.size();
}

inline unsigned SteinerSubdivision::getEndOfStVts() const
{
	return (unsigned)m_stVts.size();
}

inline TriEdge SteinerSubdivision::getStEdge(const unsigned _i) const
{
	return m_stEdges[_i];
}

inline TriEdge SteinerSubdivision::getOrigTriEdge(const unsigned _ei) const
{
	return m_origG_ptr->edges[_ei];
}

inline int SteinerSubdivision::getStEdgeIndex(const TriEdge _st_e) const
{
	auto find_it = m_stEdge_idx_map.find(_st_e);
	return find_it != m_stEdge_idx_map.end() ? find_it->second : -1;
}

inline bool SteinerSubdivision::isTriEdge(const unsigned _ei) const
{
	return _ei < m_triEdge_stVtsRange_list.size();
}

inline bool SteinerSubdivision::triEdgeSubdivided(const unsigned _i) const
{
	return m_triEdge_stVtsRange_list[_i].second > 0;
}

/*
<start index, # st vts> for the query edge
*/
inline void SteinerSubdivision::getStVertIndicesOnTriEdge(
	const unsigned _ei, vector<unsigned>& _stVts) const
{
	const auto& pr = m_triEdge_stVtsRange_list[_ei];
	_stVts.reserve(pr.second);
	for (unsigned i = pr.first; i < pr.first+pr.second; ++i)
	{
		_stVts.push_back(i);
	}
}
// TODO: OPTIMIZE_AWAY_THIS
inline void SteinerSubdivision::getStVertIndicesOnTriEdge(
	const TriEdge& _e, vector<unsigned>& _stVts) const
{
	const auto& pr = m_triEdge_stVts_map.find(_e)->second;
	_stVts.reserve(pr.second);
	for (unsigned i = pr.first; i < pr.first+pr.second; ++i)
	{
		_stVts.push_back(i);
	}
}

inline void SteinerSubdivision::getVtsOnTriEdge(const unsigned _ei, vector<unsigned>& _vts) const
{
	auto& e = m_origG_ptr->edges[_ei];
	_vts.push_back(e[0]);
	getStVertIndicesOnTriEdge(_ei, _vts);
	_vts.push_back(e[1]);
}

/*
the position of the st vts on the query edge
*/
inline void SteinerSubdivision::getStVtsPosOnTriEdge(
	const unsigned _ei, vector<TriPoint>& _stVts) const
{
	const auto& pr = m_triEdge_stVtsRange_list[_ei];
	_stVts.resize(pr.second);
	if (!_stVts.empty())
	{
		std::copy(m_stVts.begin()+pr.first, m_stVts.begin()+pr.second, _stVts.end());
	}
}

// TODO: OPTIMIZE_AWAY_THIS
inline void SteinerSubdivision::getStVtsPosOnTriEdge(
	const TriEdge& _e, vector<TriPoint>& _stVts) const
{
	const auto& pr = m_triEdge_stVts_map.find(_e)->second;
	_stVts.reserve(pr.second);
	for (unsigned i = pr.first; i < pr.second; ++i)
	{
		_stVts.push_back(m_stVts[i]);
	}
}

inline unsigned SteinerSubdivision::sizeOfVts() const
{
	return (unsigned)m_stVts.size();
}
inline unsigned SteinerSubdivision::sizeOfStEdges() const
{
	return m_n_stEdges;
	//return (unsigned)m_stEdges.size();
}

inline const vector<TriPoint>& SteinerSubdivision::getAllVts() const
{
	return m_stVts;
}

inline const vector<TriEdge>& SteinerSubdivision::getStEdges() const
{
	//static_assert(1==0, (std::string("deprecated: ")+std::string(__func__)).c_str() );
	std::cout << __FUNCTION__ << " is deprecated!" <<std::endl;
	exit(-1);
	return vector<TriEdge>();
}

inline float SteinerSubdivision::getIntervalSize() const
{
	return m_d;
}

#endif