#ifndef SURF_FUNC_H
#define SURF_FUNC_H

#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseLU>

#include "allTypes.h"
#include "commonUtils.h"
#include "graphWeights.h"
#include "graphDiffusion.h"

class SurfaceFunc
{
public:
	enum SurfFuncType {SHAPE_DIAM, SHAPE_WIDTH, SHAPE_WIDTH_EXTREMITY, SHAPE_LENGTH_EXTREMITY};
public:
	static const float INVALID_VAL;
public:
	SurfaceFunc ();
	~SurfaceFunc();

	void reset();
	/* init the surface function by passing in 
	   the steiner graph struct, and the dual line structure */
	void setup(
		const std::shared_ptr<SteinerGraph>& _stg_ptr, 
		const std::shared_ptr<TriMesh>& _3d_surf);
	void setCurrentFunction();
	void getSurfaceFunction(
		SurfFuncType _func, 
		vector<float>& _scalar_vals, 
		float _smooth_ratio, 
		bool _do_diffuse/* = false*/);
	// vert id _i on orig surf -> closest vert id of MC _i mapped to
	// return -1 if _i out of range or correspondence is not computed
	int getOrigToMC(int _i);

	/* compute gradient of bt2 per vertex */
	void computeNextVertAlongBT2Gradient();

	/* output diffuse system to *mathematica* */
	void outputDiffuseSystem( const vector<float>& _scalars );

	
	/* fill in m_closest_vert the correspondence 
	   between dual vert and its closest vert on orig surf */
	void computeCorrespondenceWithKNNSearch(
		int _n
		);
	void computeCorrespondenceWithRangeSearch(
		float _eps, int _k = -1
		);
	/* build & solve linear system used for surface scalar function diffusion */
	typedef Eigen::SparseMatrix<double> SpMat;
	typedef Eigen::Triplet<double> SpTriplet;
	void setupDiffusionSystem();

private:
	/* find local minimum from given vert _vi following bt2 gradient */
	int find_local_maximum(unsigned _vi, float _eps);
	/* compute the mean-value coordinate weight for given vi with edge <vi, vj> */
	float compute_mvc_weight(unsigned _vi, unsigned _vj, std::shared_ptr<MyGraph> _g);
	/* setup diffuser */
	void setup_diffuser();
	/* build & solve linear system used for surface scalar function diffusion */
	void solve_diffusion_system( 
		const vector<float>& _known,
		vector<float>& _unknown );
	/* carry out diffuse on given scalar function */
	void diffuse_scalar_field( vector<float>& _scalar_vals );

	void convert_to_shape_diam_eps_range(const float _smooth_ratio, float& _eps);
	void convert_to_shape_width_eps_range(const float _smooth_ratio, float& _eps);
	void convert_to_shape_width_extremity_eps_range(const float _smooth_ratio, float& _eps);
	void get_shape_diam(vector<float>& _scalar_vals, float _eps);
	void get_shape_width(vector<float>& _scalar_vals, float _eps);
	void get_shape_width_extremity(vector<float>& _scalar_vals, float _eps);
	
private:
	// compute gaussian with given params
	inline float gaussian( float _x, float _mu, float _sigma )
	{
		return 1.0f / (_sigma*std::sqrt(2.0f*M_PI)) * 
			std::exp( -0.5f * (_x-_mu)*(_x-_mu) / (_sigma*_sigma) );
	}
	// obtain smoothed field _scalars using the medial curve
	void obtain_smooth_field_using_mc(
		SurfFuncType _type, 
		const vector<float>& _mc_field, 
		vector<float>& _scalar_vals, 
		float _eps);

	std::shared_ptr<SteinerGraph> m_stg_ptr;
	std::shared_ptr<TriMesh> m_origSurf_mesh;
	std::shared_ptr<MyGraph> m_origSurf_graph;

	// the correspondence between MC vert and closest vert on orig surf
	vector<int> m_closest_vert_for_MC;
	// vert on orig surf -> vert on MC (inverse of the above mapping)
	vector<int> m_closest_vert_for_origSurf;
	// "close" MC points for each point on original surface
	// how close? depends on implementation of the corresp. building routine
	vector< vector<std::pair<unsigned, float>> > m_close_vts_for_origSurf;
	// the next vert of the local neighborhood along gradient of bt2
	//vector<int> m_next_vert_by_bt2;

	GraphDiffuser m_diffuser;

	/// the diffusion system components
	// static part
	SpMat A_sm, C_sm;
	//std::shared_ptr<Eigen::SimplicialCholesky<SpMat>> decomp_of_A;
	//Eigen::SimplicialCholesky<SpMat> decomp_of_A;
	Eigen::SparseLU<SpMat> decomp_of_A;
	// dynamic part 
	Eigen::VectorXd f;
};

#endif