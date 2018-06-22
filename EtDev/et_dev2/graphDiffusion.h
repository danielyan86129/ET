// Diffusion on our graph structure
// 
// Copyright (C) 2018 Yajie Yan <danielyan86129@hotmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef GRAPH_DIFFUSION_H
#define GRAPH_DIFFUSION_H

#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
//#include <Eigen/SparseLU>
//#include <Eigen/SparseQR>
//#include <Eigen/SparseCholesky>
#include <Eigen/OrderingMethods>

#include "allTypes.h"
#include "commonUtils.h"
#include "graphWeights.h"

class GraphDiffuser
{
public:
	enum Errcode {SUCCESS = 0, FAILED = 1, ZERO_UNKNOWNS = 2};

	GraphDiffuser();
	~GraphDiffuser();

	///
	/// public interfaces
	///

	/* setup diffusion so that later new diffusions can be solved quickly
	** provided that unknowns do not change */
	Errcode setupDiffuseSystem(
		const MyGraph* const _g,
		GraphWeightFunc& _wFunc, /*the weighting functor*/
		const vector<bool>& _is_known /*the indicating list*/
		);
	Errcode setupDiffuseSystem(
		const vector<TriPoint>& _vts,
		const vector<vector<int>>& _vts_adjacency,
		GraphWeightFunc& _wFunc, /*the weighting functor*/
		const vector<bool>& _is_known /*the indicating list*/
		);
	/* quickly solve for unknowns given new scalar values 
	** Note: only works when diffusion sytem has been setup 
	** and the system structure (i.e. unknowns) doesn't change*/
	Errcode solve( vector<float>& _scalar_vals );

	/* output diffuse system to *mathematica* */
	void outputDiffuseSystem( 
		const vector<TriPoint>& _g_vts,
		const vector<TriEdge>& _g_edges,
		const vector<float>& _scalars 
		);

	/* perform diffusion on input graph */
	Errcode diffuseOnGraph(
		const MyGraph* const _g, /*the graph to operate on*/
		GraphWeightFunc& _wFunc, /*the weighting functor*/
		const vector<bool>& _is_known, /*the indicating list*/
		vector<float>& _scalars /*the list of known and unknown scalars. unknowns will be filled when function returns*/
		);

	/* output diffuse system to *mathematica* */
	void outputDiffuseSystem( const vector<float>& _scalars );
	void outputMatA( /*const vector<float>& _scalars*/ );

private:
	typedef Eigen::SparseMatrix<double> SpMat;
	typedef Eigen::Triplet<double> SpTriplet;

	/* build linear system used for surface scalar function diffusion */
	//void setup_diffuse_system(
	//	const MyGraph* const _g, /*the graph to operate on*/
	//	GraphWeightFunc& _wFunc, /*the weighting functor*/
	//	const vector<bool>& _flags /*the indicating list*/
	//	);
	void solve_diffusion_system( 
		const vector<float>& _known,
		vector<float>& _unknown );

private:
	/// the diffusion system components
	// static part
	SpMat A_sm, C_sm;
	//std::shared_ptr<Eigen::SimplicialCholesky<SpMat>> decomp_of_A;
	//Eigen::SimplicialCholesky<SpMat> decomp_of_A;
	//// Very slow decomposition! (but works on any rectangular mat)
	//Eigen::SparseQR<SpMat, Eigen::COLAMDOrdering<int>> decomp_of_A;
	// fast decomposition (for square)
	Eigen::SparseLU<SpMat> decomp_of_A;
	/*Eigen::LeastSquaresConjugateGradient<SpMat> decomp_of_A;*/
	// dynamic part 
	Eigen::VectorXd f;
	// a list of graph edges (debug purpose)
	vector<TriPoint> m_gVts;
	vector<TriEdge> m_gEdges;
	// indicator: whether corresponding node has known value or not
	vector<bool> m_isKnown;
	// Often, input graph contains nodes w/o incident edges.
	// These nodes should be excluded from the diffusion system, 
	// by not including corresponding entries to A, C, and f.
	// In this case each node may needs a new id/order (stored in neworder_map)
	// 
	// new order map: nodes' new id -> its original idx
	vector<unsigned> m_new_orig_vertIdx_map;
};

#endif