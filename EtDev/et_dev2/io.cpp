#include "io.h"

namespace et
{
	void keepValidRadii( const vector<TriPoint>& _centers, vector<float>& _radii, const shared_ptr<MyMesh>& _MA )
	{
		vector<float> good_radii;
		KNNTree *tree = new KNNTree();
		for ( unsigned vi = 0; vi < _centers.size(); ++vi )
		{
			const auto& p = _centers[ vi ];
			tree->insert( Point_and_uint( P3( p[ 0 ], p[ 1 ], p[ 2 ] ), vi ) );
		}

		// keep the radii of the Nearest centers for MA vertices 
		for ( unsigned vi = 0; vi < _MA->vertices.size(); ++vi )
		{
			const auto& query = _MA->vertices[ vi ];
			K_neighbor_search search( *tree, P3( query[ 0 ], query[ 1 ], query[ 2 ] ), 1 );
			auto it = search.begin();
			auto vj = boost::get<1>( it->first );
			good_radii.push_back( _radii[ vj ] );
		}

		_radii = good_radii;
	}

	bool loadMesh(
		const char * _ma_file, shared_ptr<MyMesh>& _m_MA,
		const char * _r_file, vector<float>* _radii,
		const char * _orig3d_file, shared_ptr<TriMesh>& _m_orig3d,
		trimesh::XForm<double> &_trans_mat,
		float _eps, float _pert, bool _remove_dup_faces )
	{
		int ret = true;
		bool radii_ready = false;

		cout << "-----------Loading MA Mesh File-----------" << endl;
		string ma_file = _ma_file;
		_m_MA = shared_ptr<MyMesh>( MyMesh::read( _ma_file ) );
		if ( !_m_MA )
		{
			cout << "Error reading file " << _ma_file << endl;
			return false;
		}

		_m_MA->need_bbox();
		cout << "MA mesh bbox size: " << _m_MA->bbox.size().max() << endl;

		cout << "-----------Loading Orig 3d Mesh File-----------" << endl;
		_m_orig3d = shared_ptr<TriMesh>( MyMesh::read( _orig3d_file ) );
		if ( !_m_orig3d )
		{
			cout << "Error reading file " << _orig3d_file << endl;
			return false;
		}
		else
		{
			_m_orig3d->need_bbox();
			cout << "Original surface mesh bbox size: " << _m_orig3d->bbox.size().max() << endl;

			// perturb degenerate triangle faces
			_m_MA->need_bbox();
			_pert = _m_MA->bbox.radius() * _pert;
			_eps = _m_MA->bbox.radius() * _eps;
			TriFace f;
			unsigned pert_cnt = 0;
			cout << "pert. amount: " << _pert << endl;
			for ( unsigned fi = 0; fi < _m_MA->faces.size(); ++fi )
			{
				f = _m_MA->faces[ fi ];
				if ( util::is_degenerate( _m_MA->vertices[ f[ 0 ] ], _m_MA->vertices[ f[ 1 ] ], _m_MA->vertices[ f[ 2 ] ], _pert, _eps ) )
				{
					_m_MA->faces[ fi ] = f;
					pert_cnt++;
				}
			}
			cout << endl;
			cout << pert_cnt << " degenerate faces detected." << endl;

			_m_MA->need_normals();
			_m_MA->need_faces();
			_m_MA->need_bbox();
			_m_orig3d->need_normals();
			_m_orig3d->need_faces();
			_m_orig3d->need_bbox();

			// re-order indices in each face
			for ( auto & f : _m_MA->faces )
				f = util::makeFace( f );
			for ( auto & f : _m_orig3d->faces )
				f = util::makeFace( f );

			cout << endl << "Loading meshes done." << endl << endl;
		}

		if ( _remove_dup_faces )
		{
			cout << "------------Clean Duplicate Faces------------" << endl;

			/*struct FaceCompare
			{
				bool operator () ( const TriFace& a, const TriFace& b )
				{
					if ( a[ 0 ] == b[ 0 ] )
						if ( a[ 1 ] == b[ 1 ] )
							return a[ 2 ] < b[ 2 ];
						else
							return a[ 1 ] < b[ 1 ];
					else
						return a[ 0 ] < b[ 0 ];
				}
			};*/
			set<TriFace, FaceCompare> faces_map;
			vector<bool> faces_to_remove( _m_MA->faces.size(), false );
			int remove_cnt = 0;

			for ( unsigned fi = 0; fi < _m_MA->faces.size(); ++fi )
			{
				TriFace& f = _m_MA->faces[ fi ];
				//std::sort( f.v, f.v + 3 );
				f = util::makeFace( f );
				int pre_sz = faces_map.size();
				faces_map.insert( f );
				if ( faces_map.size() == pre_sz )
				{
					remove_cnt++;
					faces_to_remove[ fi ] = true;
				}
			}
			trimesh::remove_faces( _m_MA.get(), faces_to_remove );
			cout << "# duplicate faces removed: " << remove_cnt << endl;
			faces_map.clear();

			cout << endl << "Clean Duplicate Faces Done." << endl << endl;
		}

		/*cout << "Cleaning Unused Vertices... " << endl << endl;
		trimesh::remove_unused_vertices( _m_MA.get() );
		cout << endl << "Clean Unused Vertices Done." << endl << endl;*/

		// center both original surface and MA mesh to world origin
		// NOTE: must use the same translation vector
		_m_MA->need_bbox();
		_m_orig3d->need_bbox();
		auto trans_vec = -_m_orig3d->bbox.center();
		auto translate_mat = trimesh::xform::trans( trans_vec[ 0 ], trans_vec[ 1 ], trans_vec[ 2 ] );
		_trans_mat = translate_mat;
		// only translate MA mesh if it's not cleaned
		/*bool is_clean_MA = std::string(_ma_file).find(".clean.") != std::string::npos;
		if (!is_clean_MA)
		trimesh::apply_xform(_m_MA.get(), translate_mat);
		else
		{
		cout << "A cleaned MA loaded." << endl;
		}*/
		// translate the MA mesh anyway
		trimesh::apply_xform( _m_MA.get(), translate_mat );
		trimesh::apply_xform( _m_orig3d.get(), translate_mat );

		if ( _radii )
		{
			cout << "-----------Loading accompanied radii-----------" << endl;
			vector<TriPoint> centers;
			radii_ready = loadRadiiFromRfile( _r_file, *_radii, centers );
			if ( radii_ready )
			{
				// need to translate centers in the radii file to keep things consistent
				for ( size_t i = 0; i < centers.size(); ++i )
					centers[ i ] = translate_mat * centers[ i ];
				// keep only valid radii (there could be more radii than # vertices in MA)
				keepValidRadii( centers, *_radii, _m_MA );
				cout << "Loading accompanied radii done." << "(" << _radii->size() << " loaded)" << endl;
			}
			else
			{
				cout << "Error: Failed to load radii." << endl;
			}
		}

		// 	cout << "-----------preprocessing mesh-----------" << endl;
		// 	preprocess(_m_MA, *_radii, _eps);
		// 	cout << "preprocessing done" << endl;

		return true;
	}

	bool loadRadiiFromSAT( const char *_radii_file, vector<float>& _radii, vector<TriPoint>& _centers )
	{
		_centers.clear();
		_radii.clear();

		ifstream infile( _radii_file );
		if ( !infile.is_open() )
		{
			cout << "Error reading file " << _radii_file << endl;
			return false;
		}

		// get # radii to read from first line
		string dump_s;
		float dump_f;
		unsigned n_radii;

		infile >> dump_s;
		infile >> n_radii;
		infile >> dump_f >> dump_f;
		_centers.reserve( n_radii );
		_radii.reserve( n_radii );

		float x, y, z, r;
		unsigned i = 0;
		while ( i < n_radii )
		{
			infile >> x >> y >> z >> r;
			_centers.push_back( TriPoint( x, y, z ) );
			_radii.push_back( r );
			++i;
		}

		infile.close();

		return true;
	}

	bool loadRadiiFromRfile( const char *_radii_file, vector<float>& _radii, vector<TriPoint>& _centers )
	{
		ifstream infile( _radii_file );
		if ( !infile.is_open() )
		{
			cout << "Error couldn't open radii file " << _radii_file << endl;
			return false;
		}

		cout << "Loading radii field from " << _radii_file << endl;
		int n_radii;
		float x, y, z, r;
		infile >> n_radii;
		_radii.resize( n_radii );
		_centers.resize( n_radii );
		for ( int i = 0; i < n_radii; ++i )
		{
			infile >> x >> y >> z >> r;
			_radii[ i ] = r;
			_centers[ i ] = TriPoint( x, y, z );
		}

		infile.close();
		cout << "Loading radii completed!" << endl;

		return true;
	}
}