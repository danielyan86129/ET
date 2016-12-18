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
		float _eps, float _pert );

	// _radii will be of size 0 if loading fails
	bool loadRadiiFromSAT( const char *_radii_file, vector<float>& _radii, vector<TriPoint>& _centers );
	bool loadRadiiFromRfile( const char *_radii_file, vector<float>& _radii, vector<TriPoint>& _centers );

}