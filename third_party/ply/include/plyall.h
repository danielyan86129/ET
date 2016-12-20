#pragma once
#include <map>
#include <vector>
#include <string>
#include <ply.h>

namespace ply 
{
	using std::map;
	using std::vector;
	using std::string;
	enum status { FILE_NOT_OPEN, FAILURE, SUCCESS };
	class PLYreader
	{
	public:
		PLYreader();
		~PLYreader();

		//
		// open the file for reading
		bool open(const char* _fname);
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
			vector<F>& _faces );
		//
		// closes the reader.
		void close();
	private:
		PlyFile* m_ply;
	};

	class PLYwriter
	{
	public:
		PLYwriter();
		~PLYwriter();

		//
		// open a file for writing
		bool open( const char * _fname,
			int _n_elems, char** _elem_names, int _mode );
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
			const vector<F>& _faces_to_write );
		//
		// closes the writer.
		void close();
	private:
		PlyFile* m_ply;
	};
}

namespace ply
{
	template<typename V, typename E, typename F>
	status PLYreader::read( const char* _fname,
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

	template<typename V, typename E, typename F>
	status PLYwriter::write( const char* _fname,
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

		if (!open( _fname, cnt, names, PLY_ASCII ))
			return FILE_NOT_OPEN;

		const map<string, PlyProperty>* props[] = { &_v_prop, &_e_prop, &_f_prop };
		const int sizes[] = { _vts_to_write.size(), _edges_to_write.size(), _faces_to_write.size() };
		// get header ready
		for ( auto i = 0; i < 3; ++i )
		{
			if (!names[i])
				continue;
			describe_element_ply( m_ply, names[i], sizes[i] );
			for ( auto it = props[ i ]->begin(); it != props[ i ]->end(); ++it )
			{
				describe_property_ply( m_ply, &it->second );
			}
		}
		header_complete_ply( m_ply );
		// write data
		if (_export_v )
		{
			put_element_setup_ply( m_ply, names[ 0 ] );
			for ( auto j = 0; j < _vts_to_write.size(); ++j )
				put_element_ply( m_ply, (void*)&_vts_to_write[ j ] );
		}
		if (_export_e)
		{
			put_element_setup_ply( m_ply, names[1] );
			for ( auto j = 0; j < _edges_to_write.size(); ++j )
				put_element_ply( m_ply, (void*)&_edges_to_write[ j ] );
		}
		if (_export_f )
		{
			put_element_setup_ply( m_ply, names[2] );
			for ( auto j = 0; j < _faces_to_write.size(); ++j )
				put_element_ply( m_ply, (void*)&_faces_to_write[ j ] );
		}

		close();
		return SUCCESS;
	}
}