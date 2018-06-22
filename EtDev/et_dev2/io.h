// I/O related functions
// 
// Copyright (C) 2018 Yajie Yan <danielyan86129@hotmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#pragma once
#include "commonUtils.h"
#include "allTypes.h"
namespace et
{
	bool loadMesh(
		const char * _ma_file, shared_ptr<MyMesh>& _m_MA,
		const char * _r_file, vector<float>* _radii,
		const char * _orig3d_file, shared_ptr<TriMesh>& _m_orig3d,
		trimesh::XForm<double> &_trans_mat,
		float _eps, float _pert, bool _remove_dup_faces );

	// _radii will be of size 0 if loading fails
	bool loadRadiiFromSAT( const char *_radii_file, vector<float>& _radii, vector<TriPoint>& _centers );
	bool loadRadiiFromRfile( const char *_radii_file, vector<float>& _radii, vector<TriPoint>& _centers );

}