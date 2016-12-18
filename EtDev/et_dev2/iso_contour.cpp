#include "iso_contour.h"

void IsoContourExtractor::marchTriangle( 
	const vector<TriPoint>& _vts_of_faces, 
	const vector<float>& _vals, 
	const float _iso_value,
	vector<TriPoint>& _iso_vts, 
	//vector<TriEdge>& _iso_edges, 
	vector<TriPoint>& _vts_for_newFaces, 
	vector<float>& _scalar_for_newFaces_vts/*,
	vector<TriFace>& _newFaces*/
	)
{
	TriFace f;
	float face_vals[3];
	bool e_has_iso[3]; // whether each edge of face has iso vert or not
	TriPoint iso_vts_on_f[3]; // the iso point on each edge
	for ( size_t fi = 0; fi < _vts_of_faces.size()/3; ++fi )
	{
		f = util::makeFace(3*fi, 3*fi+1, 3*fi+2);
		for ( size_t vi = 0; vi < 3; ++vi )
		{
			face_vals[vi] = _vals[f[vi]];
			e_has_iso[vi] = false;
		}

		// for each edge, locate the potential iso vert on it
		for ( size_t vi = 0; vi < 3; ++vi )
		{
			bool same_sign = 
				( above(face_vals[vi], _iso_value) && above(face_vals[(vi+1)%3], _iso_value) )
				|| ( below(face_vals[vi], _iso_value) && below(face_vals[(vi+1)%3], _iso_value) );
			if ( !same_sign )
			{
				new_iso_vert(_vts_of_faces[f[vi]], _vts_of_faces[f[(vi+1)%3]],
					face_vals[vi], face_vals[(vi+1)%3], _iso_value, iso_vts_on_f[vi]
				);
				e_has_iso[vi] = true;
			}
			else
			{
				e_has_iso[vi] = false;
			}
		}

		// There should be only two edges that contain iso verts
		int iso_cnt = e_has_iso[0]+e_has_iso[1]+e_has_iso[2];
		assert(iso_cnt == 2);

		// locate those two edges that contain iso verts 
		// and the shared vertex of the two edges
		TriEdge e1, e2; // the two edges
		int e1_i, e2_i; // index for the two edges
		int shared_v; // the shared vert between the edges
		if (e_has_iso[0] && e_has_iso[1]) // 0 & 1th edges of face 
		{
			e1 = TriEdge(f[0], f[1]);
			e2 = TriEdge(f[1], f[2]);
			e1_i = 0; e2_i = 1;
			shared_v = 1;
		}
		else if (e_has_iso[1] && e_has_iso[2]) // 1 & 2th edges of face 
		{
			e1 = TriEdge(f[1], f[2]);
			e2 = TriEdge(f[2], f[0]);
			e1_i = 1; e2_i = 2;
			shared_v = 2;
		}
		else // 0 & 2th edges of face 
		{
			e1 = TriEdge(f[0], f[1]);
			e2 = TriEdge(f[0], f[2]);
			e1_i = 0; e2_i = 2;
			shared_v = 0;
		}

		// save an iso-edge
		_iso_vts.push_back(iso_vts_on_f[e1_i]);
		_iso_vts.push_back(iso_vts_on_f[e2_i]);

		// the iso values at the other end for each edge
		int e1_other_end = util::otherEnd(e1, f[shared_v]);
		int e2_other_end = util::otherEnd(e2, f[shared_v]);
		float scalar_e_other_end = _vals[e1_other_end];
		float scalar_f_other_end = _vals[e2_other_end];

		// if the shared vert is ahead of iso contour, 
		// then a triangle remains
		// else, a quad (2 triangles) remains
		if ( above(_vals[f[shared_v]], _iso_value) )
		{
			_vts_for_newFaces.push_back( _vts_of_faces[f[shared_v]] );
			_vts_for_newFaces.push_back( iso_vts_on_f[e1_i] );
			_vts_for_newFaces.push_back( iso_vts_on_f[e2_i] );
			_scalar_for_newFaces_vts.push_back( face_vals[shared_v] );
			_scalar_for_newFaces_vts.push_back( _iso_value );
			_scalar_for_newFaces_vts.push_back( _iso_value );
		}
		else
		{
			_vts_for_newFaces.push_back( _vts_of_faces[e1_other_end] );
			_vts_for_newFaces.push_back( iso_vts_on_f[e1_i] );
			_vts_for_newFaces.push_back( iso_vts_on_f[e2_i] );
			_scalar_for_newFaces_vts.push_back( _vals[e1_other_end] );
			_scalar_for_newFaces_vts.push_back( _iso_value );
			_scalar_for_newFaces_vts.push_back( _iso_value );

			_vts_for_newFaces.push_back( _vts_of_faces[e2_other_end] );
			_vts_for_newFaces.push_back( _vts_of_faces[e1_other_end] );
			_vts_for_newFaces.push_back( iso_vts_on_f[e2_i] );
			_scalar_for_newFaces_vts.push_back( _vals[e2_other_end] );
			_scalar_for_newFaces_vts.push_back( _vals[e1_other_end] );
			_scalar_for_newFaces_vts.push_back( _iso_value );
		}
	}
}