#pragma once
#include <map>
#include <iostream>
#include <vector>
#include <string>
#include <ply.h>

namespace ply
{
	using std::map;
	using std::vector;
	using std::string;
	enum status { FILE_NOT_OPEN, FAILURE, SUCCESS };
	struct Vertex
	{
		float x; float y; float z;
		unsigned char r, g, b;
		float s; // measure
	};
	struct Edge
	{
		int v1; int v2;
		unsigned char r, g, b;
		float s; // measure
	};
	struct Face
	{
		unsigned char nvts;
		int verts[ 3 ];
		unsigned char r, g, b;
		float s; // measure
	};
	class PLYreader
	{
	public:
		PLYreader() {};
		~PLYreader() {};

		//
		// open the file for reading
		bool open( const char* _fname )
		{
			FILE* fp = fopen( _fname, "r" );
			if ( !fp )
				return false;
			m_ply = read_ply( fp );
			return true;
		}
		//
		// try to read in a list of expected elements (e.g. "vertex", "edge", "face", etc.)
		// @param _v/e/f_props specifies the format of an element
		// @return list of each type of element, resp. empty if that element is missing.
		template<typename V, typename E, typename F>
		status read( const char* _fname,
			const map<string, PlyProperty>& _v_props,
			const map<string, PlyProperty>& _e_props,
			const map<string, PlyProperty>& _f_props,
			vector<V>& _vts,
			vector<E>& _edges,
			vector<F>& _faces )
		{
			if ( !open( _fname ) )
				return FILE_NOT_OPEN;
			int elem_size = 0;
			// find the element to setup
			char* elem_name = nullptr;
			for ( auto i = 0; i < m_ply->num_elem_types; ++i )
			{
				elem_name = setup_element_read_ply( m_ply, i, &elem_size );
				if ( equal_strings( "vertex", elem_name ) )
				{
					if ( elem_size == 0 )
						continue;
					std::cout << "about to read element: " << elem_name
						<< " (" << elem_size << ")" << std::endl;
					// prepare the format of the element to be read
					for ( auto j = 0; j < m_ply->elems[ i ]->nprops; ++j )
					{
						auto prop = m_ply->elems[ i ]->props[ j ];
						auto find_it = _v_props.find( prop->name );
						if ( find_it != _v_props.end() )
						{
							std::cout << "found property: " << prop->name << std::endl;
							setup_property_ply( m_ply, &( find_it->second ) );
						}
					}

					// read in all elements
					V v;
					_vts.clear();
					_vts.reserve( elem_size );
					for ( auto i = 0; i < elem_size; ++i )
					{
						get_element_ply( m_ply, (void *)&v );
						_vts.push_back( v );
					}
				}
				else if ( equal_strings( "edge", elem_name ) )
				{
					if ( elem_size == 0 )
						continue;
					std::cout << "about to read element: " << elem_name
						<< " (" << elem_size << ")" << std::endl;
					// prepare the format of the element to be read
					for ( auto j = 0; j < m_ply->elems[ i ]->nprops; ++j )
					{
						auto prop = m_ply->elems[ i ]->props[ j ];
						auto find_it = _e_props.find( prop->name );
						if ( find_it != _e_props.end() )
						{
							std::cout << "found property: " << prop->name << std::endl;
							setup_property_ply( m_ply, &( find_it->second ) );
						}
					}
					// read in all elements
					E e;
					_edges.clear();
					_edges.reserve( elem_size );
					for ( auto i = 0; i < elem_size; ++i )
					{
						get_element_ply( m_ply, (void *)&e );
						_edges.push_back( e );
					}
				}
				else if ( equal_strings( "face", elem_name ) )
				{
					if ( elem_size == 0 )
						continue;
					std::cout << "about to read element: " << elem_name
						<< " (" << elem_size << ")" << std::endl;
					// prepare the format of the element to be read
					for ( auto j = 0; j < m_ply->elems[ i ]->nprops; ++j )
					{
						auto prop = m_ply->elems[ i ]->props[ j ];
						auto find_it = _f_props.find( prop->name );
						if ( find_it != _f_props.end() )
						{
							std::cout << "found property: " << prop->name << std::endl;
							setup_property_ply( m_ply, &( find_it->second ) );
						}
						else
						{
							// when we don't want the prop, greg turk's ply reader didn't give internal type
							// a valid value. We make it the same as external type to avoid crashing.
							prop->internal_type = prop->external_type; 
						}
					}
					// read in all elements
					F f;
					_faces.clear();
					_faces.reserve( elem_size );
					for ( auto i = 0; i < elem_size; ++i )
					{
						get_element_ply( m_ply, (void *)&f );
						_faces.push_back( f );
					}
				}
			}
			close();
			return SUCCESS;
		}
		//
		// closes the reader.
		void close()
		{
			if ( m_ply )
			{
				close_ply( m_ply );
				free_ply( m_ply );
			}
		}
	private:
		PlyFile* m_ply;
	};

	class PLYwriter
	{
	public:
		PLYwriter() {}
		~PLYwriter() {}

		//
		// open a file for writing
		bool open( const char * _fname,
			int _n_elems, char** _elem_names, int _mode )
		{
			FILE* fp = fopen( _fname, "w" );
			if ( !fp )
				return false;
			m_ply = write_ply( fp, _n_elems, _elem_names, _mode );
			return true;
		}
		//
		// write a list of elements to file.
		// @param _export_v/e/f if v/e/f will be written
		// @param _v/e/f_prop specifies format for v/e/f elements.
		//        provide dummy if corresponding elements will not be written.
		// @param _vts/edges/faces_to_write contains the elements to write. 
		//        provide dummy if corresponding elements will not be written.
		template<typename V, typename E, typename F>
		status write( const char* _fname,
			bool _export_v,
			bool _export_e,
			bool _export_f,
			const map<string, PlyProperty>& _v_prop,
			const map<string, PlyProperty>& _e_prop,
			const map<string, PlyProperty>& _f_prop,
			const vector<V>& _vts_to_write,
			const vector<E>& _edges_to_write,
			const vector<F>& _faces_to_write )
		{
			char* names[] = {
				_export_v ? "vertex" : nullptr,
				_export_e ? "edge" : nullptr,
				_export_f ? "face" : nullptr
			};
			int cnt = _export_v + _export_e + _export_f;

			if ( !open( _fname, cnt, names, PLY_ASCII ) )
				return FILE_NOT_OPEN;

			const map<string, PlyProperty>* props[] = { &_v_prop, &_e_prop, &_f_prop };
			const int sizes[] = { _vts_to_write.size(), _edges_to_write.size(), _faces_to_write.size() };
			// get header ready
			for ( auto i = 0; i < 3; ++i )
			{
				if ( !names[ i ] )
					continue;
				describe_element_ply( m_ply, names[ i ], sizes[ i ] );
				for ( auto it = props[ i ]->begin(); it != props[ i ]->end(); ++it )
				{
					describe_property_ply( m_ply, &it->second );
				}
			}
			header_complete_ply( m_ply );
			// write data
			if ( _export_v )
			{
				put_element_setup_ply( m_ply, names[ 0 ] );
				for ( auto j = 0; j < _vts_to_write.size(); ++j )
					put_element_ply( m_ply, (void*)&_vts_to_write[ j ] );
			}
			if ( _export_e )
			{
				put_element_setup_ply( m_ply, names[ 1 ] );
				for ( auto j = 0; j < _edges_to_write.size(); ++j )
					put_element_ply( m_ply, (void*)&_edges_to_write[ j ] );
			}
			if ( _export_f )
			{
				put_element_setup_ply( m_ply, names[ 2 ] );
				for ( auto j = 0; j < _faces_to_write.size(); ++j )
					put_element_ply( m_ply, (void*)&_faces_to_write[ j ] );
			}

			close();
			return SUCCESS;
		}
		// 
		// simply write vts, edges, and tri angle faces to a file
		/*status write(
			const char* _fname,
			const vector<ply::Vertex>& output_vts,
			const vector<ply::Edge>& output_edges,
			const vector<ply::Face>& output_faces
		);*/
		//
		// closes the writer.
		void close()
		{
			if ( m_ply )
			{
				close_ply( m_ply );
				free_ply( m_ply );
			}
		}
	private:
		PlyFile* m_ply;
	};
}
//
//namespace ply
//{

	/*ply::PLYreader::PLYreader()
	{
	}

	ply::PLYreader::~PLYreader()
	{
	}*/

	/*bool ply::PLYreader::open( const char * _fname )
	{
		FILE* fp = fopen( _fname, "r" );
		if ( !fp )
			return false;
		m_ply = read_ply( fp );
		return true;
	}*/

	/*void ply::PLYreader::close()
	{
		if ( m_ply )
		{
			close_ply( m_ply );
			free_ply( m_ply );
		}
	}*/

	//template<typename V, typename E, typename F>
	//ply::status ply::PLYreader::read( const char* _fname,
	//	const map<string, PlyProperty>& _v_props,
	//	const map<string, PlyProperty>& _e_props,
	//	const map<string, PlyProperty>& _f_props,
	//	vector<V>& _vts,
	//	vector<E>& _edges,
	//	vector<F>& _faces )
	//{
	//	if ( !open( _fname ) )
	//		return FILE_NOT_OPEN;
	//	int elem_size = 0;
	//	// find the element to setup
	//	char* elem_name = nullptr;
	//	for ( auto i = 0; i < m_ply->num_elem_types; ++i )
	//	{
	//		elem_name = setup_element_read_ply( m_ply, i, &elem_size );
	//		if ( equal_strings( "vertex", elem_name ) )
	//		{
	//			if ( elem_size == 0 )
	//				continue;
	//			std::cout << "about to read element: " << elem_name
	//				<< " (" << elem_size << ")" << std::endl;
	//			// prepare the format of the element to be read
	//			for ( auto j = 0; j < m_ply->elems[ i ]->nprops; ++j )
	//			{
	//				auto prop = m_ply->elems[ i ]->props[ j ];
	//				auto find_it = _v_props.find( prop->name );
	//				if ( find_it != _v_props.end() )
	//				{
	//					std::cout << "found property: " << prop->name << std::endl;
	//					setup_property_ply( m_ply, &( find_it->second ) );
	//				}
	//			}
	//			// read in all elements
	//			V v;
	//			_vts.clear();
	//			_vts.reserve( elem_size );
	//			for ( auto i = 0; i < elem_size; ++i )
	//			{
	//				get_element_ply( m_ply, (void *)&v );
	//				_vts.push_back( v );
	//			}
	//		}
	//		else if ( equal_strings( "edge", elem_name ) )
	//		{
	//			if ( elem_size == 0 )
	//				continue;
	//			std::cout << "about to read element: " << elem_name
	//				<< " (" << elem_size << ")" << std::endl;
	//			// prepare the format of the element to be read
	//			for ( auto j = 0; j < m_ply->elems[ i ]->nprops; ++j )
	//			{
	//				auto prop = m_ply->elems[ i ]->props[ j ];
	//				auto find_it = _e_props.find( prop->name );
	//				if ( find_it != _e_props.end() )
	//				{
	//					std::cout << "found property: " << prop->name << std::endl;
	//					setup_property_ply( m_ply, &( find_it->second ) );
	//				}
	//			}
	//			// read in all elements
	//			E e;
	//			_edges.clear();
	//			_edges.reserve( elem_size );
	//			for ( auto i = 0; i < elem_size; ++i )
	//			{
	//				get_element_ply( m_ply, (void *)&e );
	//				_edges.push_back( e );
	//			}
	//		}
	//		else if ( equal_strings( "face", elem_name ) )
	//		{
	//			if ( elem_size == 0 )
	//				continue;
	//			std::cout << "about to read element: " << elem_name
	//				<< " (" << elem_size << ")" << std::endl;
	//			// prepare the format of the element to be read
	//			for ( auto j = 0; j < m_ply->elems[ i ]->nprops; ++j )
	//			{
	//				auto prop = m_ply->elems[ i ]->props[ j ];
	//				auto find_it = _f_props.find( prop->name );
	//				if ( find_it != _f_props.end() )
	//				{
	//					std::cout << "found property: " << prop->name << std::endl;
	//					setup_property_ply( m_ply, &( find_it->second ) );
	//				}
	//			}
	//			// read in all elements
	//			F f;
	//			_faces.clear();
	//			_faces.reserve( elem_size );
	//			for ( auto i = 0; i < elem_size; ++i )
	//			{
	//				get_element_ply( m_ply, (void *)&f );
	//				_faces.push_back( f );
	//			}
	//		}
	//	}
	//	close();
	//	return SUCCESS;
	//}

	/*ply::PLYwriter::PLYwriter()
	{
	}

	ply::PLYwriter::~PLYwriter()
	{
	}*/

	/*bool ply::PLYwriter::open( const char * _fname,
		int _n_elems, char** _elem_names, int _mode )
	{
		FILE* fp = fopen( _fname, "w" );
		if ( !fp )
			return false;
		m_ply = write_ply( fp, _n_elems, _elem_names, _mode );
		return true;
	}*/

	/*void ply::PLYwriter::close()
	{
		if ( m_ply )
		{
			close_ply( m_ply );
			free_ply( m_ply );
		}
	}*/

	/*template<typename V, typename E, typename F>
	ply::status ply::PLYwriter::write( const char* _fname,
		bool _export_v,
		bool _export_e,
		bool _export_f,
		const map<string, PlyProperty>& _v_prop,
		const map<string, PlyProperty>& _e_prop,
		const map<string, PlyProperty>& _f_prop,
		const vector<V>& _vts_to_write,
		const vector<E>& _edges_to_write,
		const vector<F>& _faces_to_write )
	{
		char* names[] = {
			_export_v ? "vertex" : nullptr,
			_export_e ? "edge" : nullptr,
			_export_f ? "face" : nullptr
		};
		int cnt = _export_v + _export_e + _export_f;

		if ( !open( _fname, cnt, names, PLY_ASCII ) )
			return FILE_NOT_OPEN;

		const map<string, PlyProperty>* props[] = { &_v_prop, &_e_prop, &_f_prop };
		const int sizes[] = { _vts_to_write.size(), _edges_to_write.size(), _faces_to_write.size() };
		// get header ready
		for ( auto i = 0; i < 3; ++i )
		{
			if ( !names[ i ] )
				continue;
			describe_element_ply( m_ply, names[ i ], sizes[ i ] );
			for ( auto it = props[ i ]->begin(); it != props[ i ]->end(); ++it )
			{
				describe_property_ply( m_ply, &it->second );
			}
		}
		header_complete_ply( m_ply );
		// write data
		if ( _export_v )
		{
			put_element_setup_ply( m_ply, names[ 0 ] );
			for ( auto j = 0; j < _vts_to_write.size(); ++j )
				put_element_ply( m_ply, (void*)&_vts_to_write[ j ] );
		}
		if ( _export_e )
		{
			put_element_setup_ply( m_ply, names[ 1 ] );
			for ( auto j = 0; j < _edges_to_write.size(); ++j )
				put_element_ply( m_ply, (void*)&_edges_to_write[ j ] );
		}
		if ( _export_f )
		{
			put_element_setup_ply( m_ply, names[ 2 ] );
			for ( auto j = 0; j < _faces_to_write.size(); ++j )
				put_element_ply( m_ply, (void*)&_faces_to_write[ j ] );
		}

		close();
		return SUCCESS;
	}*/

	/*status PLYwriter::write(
		const char* _fname,
		const vector<ply::Vertex>& output_vts,
		const vector<ply::Edge>& output_edges,
		const vector<ply::Face>& output_faces
	)
	{
		std::map<std::string, PlyProperty> vert_props;
		vert_props[ "x" ] = { "x", Float32, Float32, offsetof( Vertex, x ), PLY_SCALAR, 0, 0, 0 };
		vert_props[ "y" ] = { "y", Float32, Float32, offsetof( Vertex, y ), PLY_SCALAR, 0, 0, 0 };
		vert_props[ "z" ] = { "z", Float32, Float32, offsetof( Vertex, z ), PLY_SCALAR, 0, 0, 0 };
		std::map<std::string, PlyProperty> edge_props;
		edge_props[ "vertex1" ] = { "vertex1", Int32, Int32, offsetof( Edge, v1 ), PLY_SCALAR, 0, 0, 0 };
		edge_props[ "vertex2" ] = { "vertex2", Int32, Int32, offsetof( Edge, v2 ), PLY_SCALAR, 0, 0, 0 };
		std::map<std::string, PlyProperty> face_props;
		face_props[ "vertex_indices" ] = {
			"vertex_indices", Int32, Int32, offsetof( Face, verts ),
			PLY_LIST, Uint8, Uint8, offsetof( Face,nvts ) };

		std::cout << "writing to ply file... " << std::endl;
		std::cout << "# vts/edges/faces: "
			<< output_vts.size() << "/"
			<< output_edges.size() << "/"
			<< output_faces.size() << std::endl;

		return write( _fname, true, true, true,
			vert_props, edge_props, face_props,
			output_vts, output_edges, output_faces );
	}
*/
//}