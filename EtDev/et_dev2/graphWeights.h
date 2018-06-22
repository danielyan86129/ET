// Different types of edge weight on our graph structure
// 
// Copyright (C) 2018 Yajie Yan <danielyan86129@hotmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef GRAPH_WEIGHTS_H
#define GRAPH_WEIGHTS_H

#include "allTypes.h"
#include "commonUtils.h"
#include "origGraph.h"

struct GraphWeightFunc
{
public:
	virtual float operator() ( unsigned _vi, unsigned _vj )
	{
		return 1.0f;
	}
};

///
/// mean-value-coordinates function
///
struct MVCWeight : GraphWeightFunc
{
public:
	/* set the graph where weights will be computed later */
	inline void setup( std::shared_ptr<MyGraph> _g )
	{
		m_g = _g;
	}

	///
	/// re-implement parent functions
	///
	/* compute the mean-value coordinate weight for given vi with edge <vi, vj> */
	inline float operator() ( unsigned _vi, unsigned _vj )
	{
		auto e = util::makeEdge(_vi, _vj);
		const auto& nb_faces = m_g->getNbFaces( e );
		float theta[2] = {0.0f, 0.0f};
		TriPoint p0 = m_g->vts[_vi];
		TriPoint p1 = m_g->vts[_vj];
		for ( int i = 0; i < 2; ++i )
		{
			const auto& f = m_g->faces[nb_faces[i]];
			auto vk = util::oppositeVert( f, e );

			TriPoint p2 = m_g->vts[vk];
			theta[i] = acos( 
				trimesh::normalize(vec3(p1 - p0)) DOT trimesh::normalize(vec3(p2 - p0)) 
				);
		}

		return ( std::tan(theta[0] / 2) + tan(theta[1] / 2) ) / trimesh::dist(p0, p1);
	}

private:
	std::shared_ptr<MyGraph> m_g;
};

///
/// uniform weighting function
///
struct UniformWeight : GraphWeightFunc
{
public:
	///
	/// re-implement parent functions
	///
	/* return constant weight */
	inline float operator() ( unsigned _vi, unsigned _vj )
	{
		return 1.0f;
	}
};

#endif