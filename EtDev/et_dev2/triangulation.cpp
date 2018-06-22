// Triangulation of each MA face by the burn paths
// 
// Copyright (C) 2018 Yajie Yan <danielyan86129@hotmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.

#include "triangulation.h"

void Triangulation::triangulateChordBased( 
	const vector<unsigned>& _v_sorted, 
	// connect seq storing pairs: <vi, and relative pos in _v_indices
	const vector< std::pair<unsigned,unsigned> >& _conn_seq, 
	vector<TriFace>& _tri_faces 
	)
{
	assert(_v_sorted.size() >= 3);
	
	// the flag indicating whether a vert has been connected to a tri or not
	vector<bool> connected(_conn_seq.size(), false);
	// connect the first 3 vts up to form the base triangle
	_tri_faces.push_back( util::makeFace( _conn_seq[0].first,_conn_seq[1].first,_conn_seq[2].first) );
	connected[_conn_seq[0].second] = true;
	connected[_conn_seq[1].second] = true;
	connected[_conn_seq[2].second] = true;

	// iteratively grab the next vert in conn-seq to connect
	for ( int ii = 3; ii < _conn_seq.size(); ++ii )
	{
		unsigned vi = _conn_seq[ii].first;
		auto i = _conn_seq[ii].second;

		// we are at ith vert in of the sorted list
		// now walk to its left and right until we
		// find the two verts ui, wi, that bound vi
		int j = i;
		do // walk to left
		{
			j = ( j - 1 + _v_sorted.size() ) % _v_sorted.size();
		} while ( !connected[j] );
		unsigned ui = _v_sorted[j];
		j = i;
		do // walk to right
		{
			j = ( j + 1 ) % _v_sorted.size();
		} while ( !connected[j] );
		unsigned wi = _v_sorted[j];

		// connect vi, ui, wi
		connected[i] = true;
		_tri_faces.push_back( util::makeFace(ui, vi, wi) );
	}
}