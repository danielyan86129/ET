#pragma once
#include "plyall.h"
#include <cstdio>

ply::PLYreader::PLYreader()
{
}

ply::PLYreader::~PLYreader()
{
}

bool ply::PLYreader::open( const char * _fname )
{
	FILE* fp = fopen( _fname, "r" );
	if ( !fp )
		return false;
	m_ply = read_ply( fp );
	return true;
}

void ply::PLYreader::close()
{
	if ( m_ply )
	{
		close_ply( m_ply );
		free_ply( m_ply );
	}
}

ply::PLYwriter::PLYwriter()
{
}

ply::PLYwriter::~PLYwriter()
{
}

bool ply::PLYwriter::open( const char * _fname, 
	int _n_elems, char** _elem_names, int _mode )
{
	FILE* fp = fopen( _fname, "w" );
	if ( !fp )
		return false;
	m_ply = write_ply( fp, _n_elems, _elem_names, _mode );
	return true;
}

void ply::PLYwriter::close()
{
}
