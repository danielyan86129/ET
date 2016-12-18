#include "commonUtils.h"

namespace util 
{
	void traceFace( const vector<int>& _f_of_vts, vector<ivec2>& _f_of_edges )
	{
		_f_of_edges.resize( _f_of_vts.size() );
		for ( size_t i = 0; i < _f_of_vts.size(); ++i )
		{
			int u = _f_of_vts[ i ];
			int v = _f_of_vts[ ( i + 1 ) % _f_of_vts.size() ];
			_f_of_edges[ i ] = util::makeEdge( u, v );
		}
	}

	Vec<3,double> computeFootPointsHelper(
		Vec<3,double> _v1, 
		Vec<3,double> _v2, 
		double _r1, 
		double _r2
		)
	{
		/* we want to compute the projection of q onto v1v2 */
		// a little renaming: v1, v2, q = C, A, q
		// therefore, v1v2, v1q, v2q = b, a, c
		double a = _r1, b = trimesh::len(_v1-_v2), c = _r2;
		double cos_c = (a*a + b*b - c*c)/(2.0*a*b);
		double cos_a = (b*b + c*c - a*a)/(2.0*b*c);
		double r2diff = _r1*_r1 - _r2*_r2;
		double l = b;

		double s1, s2;
		Vec<3,double> q;
		if ( cos_c <= 0.0 )
		{
			s1 = -0.5 * (l + r2diff/l);
			s2 = s1 + l;
			q = (s2*_v1 - s1*_v2) / l;
		}
		else if ( cos_a <= 0.0 )
		{
			s1 = 0.5 * (l + r2diff/l);
			s2 = s1 - l;
			q = (-s2*_v1 + s1*_v2) / l;
		}
		else
		{
			s1 = 0.5 * (l + r2diff/l);
			s2 = l - s1;
			q = (s2*_v1 + s1*_v2) / l;
		}

		return q;
	}

#include "linearAlg.h"
	IsectType computeFootPoints( 
		Vec<3,double> _v1, 
		Vec<3,double> _v2, 
		Vec<3,double> _v3, 
		double _r1, 
		double _r2, 
		double _r3,
		Vec<3,double>* _foot_pts )
	{
		typedef Vec<3,double> vec3d;
		bool is_bad_face = is_degenerate(_v1, _v2, _v3, 0.0, 0.0);
		if (is_bad_face)
		{
			auto max_double = std::numeric_limits<double>::max();
			_foot_pts[0] = vec3d(max_double);
			_foot_pts[1] = vec3d(max_double);
			return NONE_ISECT;
		}

		IsectType isect_type;

		vec3d q_temp;
		Eigen::Vector3d all_q[3];
		q_temp = computeFootPointsHelper(_v1, _v2, _r1, _r2);
		all_q[0] = Eigen::Vector3d(q_temp.data());
		q_temp = computeFootPointsHelper(_v2, _v3, _r2, _r3);
		all_q[1] = Eigen::Vector3d(q_temp.data());
		q_temp = computeFootPointsHelper(_v1, _v3, _r1, _r3);
		all_q[2] = Eigen::Vector3d(q_temp.data());

		Eigen::Vector3d v1(_v1[0], _v1[1], _v1[2]);
		Eigen::Vector3d v2(_v2[0], _v2[1], _v2[2]);
		Eigen::Vector3d v3(_v3[0], _v3[1], _v3[2]);

		Eigen::Vector3d normals[3] = {
			v1-v2, v2-v3, v1-v3
		};
		normals[0].normalize();
		normals[1].normalize();
		normals[2].normalize();
		Eigen::Vector3d c(0.0,0.0,0.0); 
		c.setZero();
		Eigen::Matrix3d matN; 
		matN.setZero();
		for (int i = 0; i < 3; ++i)
		{
			auto nni = normals[i] * normals[i].transpose();
			matN += nni;
			c += nni * (all_q[i] - v3);
		}

		Eigen::Vector3d x_alpha = v1 - v3;
		Eigen::Vector3d x_beta = v2 - v3;
		Eigen::MatrixXd matM(3, 2);
		matM << x_alpha, x_beta;
		Eigen::Matrix2d matA;
		auto matA_temp = matN * matM;
		matA << x_alpha.transpose() * matA_temp,
			x_beta.transpose() * matA_temp;

		Eigen::Vector2d b;
		b << x_alpha.transpose() * c,
			x_beta.transpose() * c;

		/*solve for alpha & beta: {alpha, beta}' = inverse-matA * b */
		auto inv_matA = util::pseudoInverse(matA);
		auto alpha_beta_vec = inv_matA * b;

		/* point x is the barycentric combination of v1v2v3 */
		auto x = 
			alpha_beta_vec(0) * v1 + 
			alpha_beta_vec(1) * v2 +
			(1.0 - alpha_beta_vec(0) - alpha_beta_vec(1)) * v3;

		/* compute height d */
		double avg_d2 = 0.0;
		Eigen::Vector3d x_v = x - v1;
		avg_d2 += std::max(_r1*_r1 - x_v.squaredNorm(), 0.0);
		x_v = x - v2;
		avg_d2 += std::max(_r2*_r2 - x_v.squaredNorm(), 0.0);
		x_v = x - v3;
		avg_d2 += std::max(_r3*_r3 - x_v.squaredNorm(), 0.0);
		avg_d2 *= 1.0 / 3.0;
		double d;
		//if ( avg_d2 >= 0.0 )
		//{
		//	/* compute two foot points */
		//	// the height of the imaginary tetrahedron
		//	d = std::sqrt(avg_d2);

		//	Eigen::Vector3d face_nrml = (v1 - v2).cross( v1 - v3 );
		//	face_nrml.normalize();
		//	Eigen::Vector3d footpt = x + d*face_nrml;
		//	_foot_pts[0] = vec3d( footpt[0], footpt[1], footpt[2] );
		//	footpt = x - d*face_nrml;
		//	_foot_pts[1] = vec3d( footpt[0], footpt[1], footpt[2] );

		//	isect_type = POINT_ISECT;
		//}
		//else
		//{
		//	// end up being here? triangle is badly shaped. 
		//	// maybe should return "no intersection".

		//	//d = (_r1 + _r2 + _r3) / 3.0;
		//	isect_type = POINT_ISECT;
		//	auto max_double = std::numeric_limits<double>::max();
		//	_foot_pts[0] = vec3d(max_double);
		//	_foot_pts[1] = vec3d(max_double);
		//}

		isect_type = POINT_ISECT;
		if ( avg_d2 >= 0.0 )
		{
			/* compute two foot points */
			// the height of the imaginary tetrahedron
			d = std::sqrt(avg_d2);
		}
		else
		{
			// end up being here? triangle is badly shaped. 
			// maybe should return "no intersection".

			d = (_r1 + _r2 + _r3) / 3.0;
		}
		Eigen::Vector3d face_nrml = (v1 - v2).cross( v1 - v3 );
		face_nrml.normalize();
		Eigen::Vector3d footpt = x + d*face_nrml;
		_foot_pts[0] = vec3d( footpt[0], footpt[1], footpt[2] );
		footpt = x - d*face_nrml;
		_foot_pts[1] = vec3d( footpt[0], footpt[1], footpt[2] );

		return isect_type;
	}
}
