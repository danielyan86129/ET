#ifndef COMMON_UTILS_H
#define COMMON_UTILS_H

/*------------------------- other utilities -------------------------------*/
#include <limits>
#include <assert.h>
#include <vector>
#include <cstdlib>
#include <Eigen/eigen>
#include <CGAL/Exact_spherical_kernel_3.h> // use sphere intersection routine

#include "allTypes.h"

using namespace std;

namespace
{
	bool is_little_endian()
	{
		union {
			uint32_t word;
			uint8_t bytes[4];
		} test_struct;
		test_struct.word = 0x1;
		if (test_struct.bytes[0] != 0)
			return true;
		else
			return false;
	}
	template<typename T>
	bool is_infinite( const T &value )
	{
		// Since we're a template, it's wise to use std::numeric_limits<T>
		//
		// Note: std::numeric_limits<T>::min() behaves like DBL_MIN, and is the smallest absolute value possible.
		//

		T max_value = std::numeric_limits<T>::max();
		T min_value = - max_value;

		return ! ( min_value <= value && value <= max_value );
	}

	template<typename T>
	bool is_nan( const T &value )
	{
		// True if NAN
		return value != value;
	}

	template<typename T>
	bool is_valid( const T &value )
	{
		return ! is_infinite(value) && ! is_nan(value);
	}

	template<typename T>
	bool almost_zero( const T &v, const T& eps)
	{
		return std::abs(v) < eps;
	}

};

/***** geom utilities *****/
namespace util
{
	inline int otherEnd(const TriEdge& _e, const int _v)
	{
		return _e[0] != _v ? _e[0] : _e[1];
	}

	inline TriEdge makeEdge(int _i, int _j)
	{
		return TriEdge(std::min(_i, _j), std::max(_i, _j));
	}

	inline TriFace makeFace(unsigned _i, unsigned _j, unsigned _k)
	{
		_i > _j ? std::swap(_i, _j) : NULL;
		_i > _k ? std::swap(_i, _k) : NULL;
		_j > _k ? std::swap(_j, _k) : NULL;
		return TriFace(_i, _j, _k);
	}

	inline TriEdge pureLoopEdge(void)
	{
		return TriEdge(-1, -1);
	}

	inline bool notPureLoop(const TriEdge& _e)
	{
		return _e[0] != -1 && _e[1] != -1;
	}

	inline bool notLoop(const TriEdge& _e)
	{
		return _e[0] != _e[1];
	}

	inline vector<TriEdge> edgesFromFace(const TriFace& _f)
	{
		vector<TriEdge> l;
		l.push_back(makeEdge(_f.v[0], _f.v[1]));
		l.push_back(makeEdge(_f.v[1], _f.v[2]));
		l.push_back(makeEdge(_f.v[2], _f.v[0]));

		return l;
	}

	inline TriEdge oppositeEdge(const TriFace& _f, int _v)
	{
		return _f[0] == _v ? makeEdge(_f[1], _f[2]) : 
			_f[1] == _v ? makeEdge(_f[0], _f[2]) : 
			makeEdge(_f[0], _f[1]);
	}

	// return the opposite vert of the given edge in a face.
	// Pre-condition: e has to be an edge in f.
	inline int oppositeVert(const TriFace& _f, const TriEdge& _e)
	{
		assert(
			(_e[0] == _f[0] || _e[0] == _f[1] || _e[0] == _f[2]) && 
			(_e[1] == _f[0] || _e[1] == _f[1] || _e[1] == _f[2])
			);
		for (int i = 0; i < 3; ++i)
		{
			if (_f[i] != _e[0] && _f[i] != _e[1])
				return _f[i];
		}

		return -1;
	}

	/* re-scale v to [0, 1], assuming it's from the [vmin, vmax] */
	float rescale(float v, float vmin, float vmax, float _exp, float _new_min, float _new_max);

	/*
	Return a RGB colour value (hot-cold scheme) given a scalar v in the range [vmin,vmax]
	In this case each colour component ranges from 0 (no contribution) to
	1 (fully saturated), modifications for other ranges is trivial.
	The colour is clipped at the end of the scales if v is outside
	the range [vmin,vmax]
	*/
	TriColor GetColour(float v,float vmin,float vmax);

	// find the barycenter of the given _p in the given triangle abc
	void find_barycentric(const TriPoint& _p, 
		const TriPoint& _a, const TriPoint& _b, const TriPoint& _c,
		float& _s, float& _t, float& _w);

	// test if the given three vts could form a valid triangle
	template <class T>
	bool valid_triangle(trimesh::Vec<3,T>& _a, trimesh::Vec<3,T>& _b, trimesh::Vec<3,T>& _c)
	{
		// deal with more-than-1-same-coordinates case specifically
		if ( ((_a[0]==_b[0]&&_b[0]==_c[0]) + (_a[1]==_b[1]&&_b[1]==_c[1]) + (_a[2]==_b[2]&&_b[2]==_c[2])) > 1 )
			return false;

		typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
		typedef K::Point_3 Point3;

		return !CGAL::collinear(
			Point3(_a[0], _a[1], _a[2]), 
			Point3(_b[0], _b[1], _b[2]), 
			Point3(_c[0], _c[1], _c[2])
			);
	}

	inline void longestEdge(
		// index of the triangle vertices
		int _i, int _j, int _k, 
		// their coordinates
		const TriPoint& _a, const TriPoint& _b, const TriPoint& _c,
		TriEdge& _longest_e, float& _len)
	{
		float len[3] = {
			trimesh::dist(_a, _b),
			trimesh::dist(_a, _c),
			trimesh::dist(_b, _c)
		};
		TriEdge edges[3] = {
			makeEdge(_i, _j),
			makeEdge(_i, _k),
			makeEdge(_j, _k)
		};
		float max_len = 0.0f;
		int max_ei = 0;
		for (int ii = 0; ii < 3; ++ii)
		{
			if (len[ii] > max_len)
			{
				max_len = len[ii];
				max_ei = ii;
			}
		}
		_longest_e = edges[max_ei];
		_len = max_len;
	}

	// return whether the two vertices have almost the same coordinates
	template <class T>
	bool is_equal(trimesh::Vec<3,T>& _x, trimesh::Vec<3,T>& _y, T _eps = (T)0.0)
	{
		T max_d = (T)0.0;
		for (int i = 0; i < 3; ++i)
		{
			max_d = std::max(std::abs(_x[i] - _y[i]), max_d);
		}
		if (max_d > (T)0.0)
			return trimesh::dist(_x, _y) < _eps;
		else
			return true;
	}

	// trace face boundary (a seq of vertices) into a seq of edges
	void traceFace( const vector<int>& _f_of_vts, vector<ivec2>& _f_of_edges );

	// take in a polygon's vertices and return
	// wether it's a degenerate polygon
	template <class T>
	bool is_degenerate(vector<trimesh::Vec<3,T>>& _vts, T _eps = (T)0.0)
	{
		_eps = std::max((T)0.0, _eps);
		trimesh::Box<3, T> bbox;
		for (unsigned vi = 0; vi < _vts.size(); ++vi)
		{
			bbox += _vts[vi];
		}
		if (bbox.radius() <= _eps)
			return true;

		// the line between min and max
		auto v_baseLine = bbox.max - bbox.min;
		trimesh::normalize(v_baseLine);

		// find the height of each vert wrt to this line
		T max_h = (T)0;
		for (unsigned vi = 0; vi < _vts.size(); ++vi)
		{
			auto v = _vts[vi];
			auto to_v = v - bbox.min;
			auto proj_len = to_v.dot(v_baseLine);
			// height
			auto h = std::sqrt( 
				std::max(
				trimesh::len2(to_v) - proj_len * proj_len,
				(T)0 ) 
				);
			max_h = std::max(h, max_h);
		}

		return max_h <= _eps;
	}

	// take in a triangle's 3 vertices and return
	// wether it's a degenerate triangle
	// also apply perturbation make it better shaped if _pert > 0.0
	template <class T>
	bool is_degenerate(trimesh::Vec<3,T>& _a, trimesh::Vec<3,T>& _b, trimesh::Vec<3,T>& _c, 
		T _perturb_amount = (T)0.0, T _eps = (T)0.0)
	{
		_eps = std::max((T)0.0, _eps);
		bool same_ab = is_equal(_a, _b, _eps);
		bool same_bc = is_equal(_b, _c, _eps);
		bool same_ac = is_equal(_a, _c, _eps);
		int same = same_ab + same_bc + same_ac;

		// if a, b, c are not identical:
		// u, v: a non-degenerate edge; w: its opposite vertex
		trimesh::Vec<3,T>* u, *v, *w; 
		trimesh::Vec<3,T> pert_v;// the vector along which perturbation is applied
		if ( same >= 2 ) 
		{
			// geometrically 3 vertices are the same point
			if (_perturb_amount <= (T)0.0) // return right away if no pert. is wanted
				return true;

			// else, we eliminate this by perturbing any vertex along the pert_v by pert_amount
			pert_v = trimesh::Vec<3,T>( 
				(T)(rand() % RAND_MAX), (T)(rand() % RAND_MAX), (T)(rand() % RAND_MAX) );
			pert_v = _perturb_amount * trimesh::normalize(pert_v);
			_a += pert_v;

			w = &_a;
			u = &_b;
			v = &_c;
		}

		// it still could be a sliver triangle. anyway set u v w properly first.
		if ( same_ab )
		{
			w = &_a;
			u = &_b;
			v = &_c;
		}
		else
		{
			w = &_c;
			u = &_a;
			v = &_b;
		}

		// now the triangle can only degenerate into one case, i.e. a line
		// since point degeneracy is eliminated already

		// compute dist from w to line defined by u, v (height of tri u,v,w)
		auto uw = (*w) - (*u);
		auto uv = (*v) - (*u);
		T proj_len = uw.dot(trimesh::normalize(uv));
		T h_sqrd = std::max(uw.dot(uw) - proj_len*proj_len, (T)0.0);
		T h = h_sqrd;
		if (h_sqrd > (T)0.0)
			h = std::sqrt(h_sqrd);

		//if (h <= _eps)
		//{
		//	cout << "sliver triangle detected! " <<*w<<","<<*u<<","<<*v<< endl; // debug
		//}

		if (_perturb_amount > (T)0.0)
		{
			// perturb w
			if ( h <= _eps )
			{
				pert_v = trimesh::Vec<3,T>( 
					(T)(rand() % RAND_MAX), (T)(rand() % RAND_MAX), (T)(rand() % RAND_MAX) );
				pert_v = _perturb_amount * trimesh::normalize(pert_v);
				*w += pert_v;

				//cout << "after perturbation: " <<*w<<","<<*u<<","<<*v<< endl; // debug
			}
		}

		return h <= _eps;
	}

	// take in a triangle's three (or any three) vertices and return 
	// whether at least two of them have the same coordinates
	template <class T>
	bool anySame(trimesh::Vec<3,T>& _a, trimesh::Vec<3,T>& _b, trimesh::Vec<3,T>& _c, 
		T _perturb_amount, T _eps = (T)0.0)
	{
		_eps = std::max((T)0.0, _eps);
		auto vert_equal = [_eps] (trimesh::Vec<3,T>& _x, trimesh::Vec<3,T>& _y) -> bool
		{
			return trimesh::dist(_x, _y) < _eps;
		};

		int same = vert_equal(_a, _b) + vert_equal(_a, _c) + vert_equal(_b, _c);
		if (same == 0 || _perturb_amount <= 0)
			return same > 0;

		trimesh::Vec<3, T> perturb;
		if ( vert_equal(_a, _b) ) 
		{
			perturb = trimesh::Vec<3, T>(
				(float)(rand() % 100000), (float)(rand() % 100000), (float)(rand() % 100000));
			perturb = _perturb_amount * trimesh::normalize(perturb);
			_b += perturb;
		}
		if ( vert_equal(_b, _c) )
		{
			perturb = trimesh::Vec<3, T>(
				(float)(rand() % 100000), (float)(rand() % 100000), (float)(rand() % 100000));
			perturb = _perturb_amount * trimesh::normalize(perturb);
			_c += perturb;
		}
		if ( vert_equal(_c, _a) )
		{
			perturb = trimesh::Vec<3, T>(
				(float)(rand() % 100000), (float)(rand() % 100000), (float)(rand() % 100000));
			perturb = _perturb_amount * trimesh::normalize(perturb);
			_a += perturb;
		}

		return same > 0;
	}

	enum IsectType{
		NONE_ISECT=0, 
		POINT_ISECT=1, 
		CIRCLE_ISECT=2, 
		SPHERE_ISECT=3, 
		INVALID_ISECT=4
	};

	typedef CGAL::Exact_spherical_kernel_3	Spherical_k;
	typedef CGAL::Point_3<Spherical_k>		Point3;
	typedef CGAL::Circle_3<Spherical_k>		Circle3;
	typedef CGAL::Sphere_3<Spherical_k>		Sphere3;
	typedef CGAL::Circular_arc_point_3<Spherical_k>	CircPoint3;
	typedef std::pair<CGAL::Circular_arc_point_3<Spherical_k>, unsigned> PointRetType;
	template <class T>
	IsectType cgalSphereIntersect(
		const Vec<3,T>& _a, const Vec<3,T>& _b, const Vec<3,T>& _c, 
		const T _ra, const T _rb, const T _rc, 
		vector<CGAL::Object>& _isecs
		)
	{
		Sphere3 s1(Point3(_a[0], _a[1], _a[2]), _ra*_ra);
		Sphere3 s2(Point3(_b[0], _b[1], _b[2]), _rb*_rb);
		Sphere3 s3(Point3(_c[0], _c[1], _c[2]), _rc*_rc);
		CGAL::intersection(
			s1,
			s2,
			s3,
			std::back_inserter(_isecs)
			);
		CGAL::set_pretty_mode(std::cout);
		/*cout << "-- s1:"<<s1<<endl
		<< "-- s2:"<<s2<<endl
		<< "-- s3:"<<s3<<endl;*/

		if (_isecs.empty())
			return NONE_ISECT;

		if ( const PointRetType* point_ret = CGAL::object_cast<PointRetType>(&_isecs[0]) )
		{
			//cout << "isect point: " << point_ret->first <<endl;
			return POINT_ISECT;
		}
		if ( CGAL::object_cast<Circle3>(&_isecs[0]) )
		{
			return CIRCLE_ISECT;
		}
		if ( CGAL::object_cast<Sphere3>(&_isecs[0]))
		{
			return SPHERE_ISECT;
		}
		return INVALID_ISECT;
	}

	IsectType computeFootPoints(
		Vec<3,double> _v1, 
		Vec<3,double> _v2, 
		Vec<3,double> _v3, 
		double _r1,
		double _r2,
		double _r3,
		Vec<3,double>* _foot_pts );

	template <class T>
	IsectType threeSphereIntersect(
		Vec<3,T> _a, Vec<3,T> _b, Vec<3,T> _c, 
		T _ra, T _rb, T _rc, 
		std::array<Vec<3,T>, 2>& _isecs
		)
	{
		T c1, c2, c3;
		Vec<3,T> ba = _b - _a;
		Vec<3,T> ba_n = ba;
		trimesh::normalize(ba_n);
		Vec<3,T> bc = _b - _c;
		Vec<3,T> bc_n = bc;
		trimesh::normalize(bc_n);
		Vec<3,T> ca = _c - _a;
		Vec<3,T> cb = _c - _b;

		Vec<3,T> n = ba_n CROSS bc_n;
		trimesh::normalize(n);

		c1 = ((T)(_ra*_ra) - (T)(_rb*_rb) - (_a DOT _a) + (_b DOT _b)) / T(2);
		c2 = (_ra*_ra - _rc*_rc - (_a DOT _a) + (_c DOT _c)) / T(2);
		c3 = (_rb*_rb - _rc*_rc - (_b DOT _b) + (_c DOT _c)) / T(2);

		Vec<3,T> u(ba DOT _a, ba DOT _b, ba DOT _c);
		Vec<3,T> v(ca DOT _a, ca DOT _b, ca DOT _c);
		Vec<3,T> w(cb DOT _a, cb DOT _b, cb DOT _c);

		// alternative: use Eigen
		// m1: (u, v, (1,1,1) )^t	cc1: (c1, c2, 1)^t
		// m2: composed of u, w, (1,1,1)	cc1: (c1, c3, 1)^t
		// m3: composed of v, w, (1,1,1)	cc1: (c2, c3, 1)^t
		Eigen::Matrix3d m1, m2, m3;
		Eigen::Vector3d cc1, cc2, cc3;
		m1 << u[0], u[1], u[2], 
			v[0], v[1], v[2], 
			1, 1, 1;
		m2 << u[0], u[1], u[2], 
			w[0], w[1], w[2], 
			1, 1, 1;
		m3 << v[0], v[1], v[2], 
			w[0], w[1], w[2], 
			1, 1, 1;
		cc1 << c1, c2, 1;
		cc2 << c1, c3, 1;
		cc3 << c2, c3, 1;

		Eigen::Vector3d bary; // barycentric coordinates
		double relative_err;  // error used to see if result good enough
		double err_t = 1e-4;  // error threshold relative error compares against
		// try solving m1*bary = cc1 first
		Eigen::FullPivLU<Eigen::Matrix3d> dec(m1);
		bary = dec.solve(cc1);
		// then verify if good enough
		relative_err = (m1*bary - cc1).norm() / cc1.norm();
		if (relative_err > err_t)
		{
			// try solving m2*bary = cc2 then
			dec.compute(m2);
			bary = dec.solve(cc2);
			relative_err = (m2*bary - cc2).norm() / cc2.norm();
		}
		if (relative_err > err_t)
		{
			// try solving m3*bary = cc3 then ...
			dec.compute(m3);
			bary = dec.solve(cc3);
			relative_err = (m3*bary - cc3).norm() / cc3.norm();
		}
		if (relative_err > err_t)
		{
			return NONE_ISECT;
		}

		Vec<3,T> q(	bary[0]*_a[0]+bary[1]*_b[0]+bary[2]*_c[0],
			bary[0]*_a[1]+bary[1]*_b[1]+bary[2]*_c[1], 
			bary[0]*_a[2]+bary[1]*_b[2]+bary[2]*_c[2] );
		//cout << "q: " << q << endl;

		Vec<3, T> qa = q - _a;
		T h = sqrt(_ra*_ra - (qa DOT qa));
		if ( !is_valid(h) )
		{
			Vec<3, T> qb = q - _b;
			h = sqrt(_rb*_rb - (qb DOT qb));
		}
		if ( !is_valid(h) )
		{
			Vec<3, T> qc = q - _c;
			h = sqrt(_rc*_rc - (qc DOT qc));
		}

		if ( !is_valid(h) )	
			return NONE_ISECT;

		_isecs[0] = (q + h*n);
		_isecs[1] = (q - h*n);
		return POINT_ISECT;
	}
};


#endif
