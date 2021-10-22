//============================================================================
//                                  I B E X
// File        : ibex_QibexOptimizer.h
// Author      :  Bertrand Neveu
// Copyright   : IMT Atlantique (France)
// License     : See the LICENSE file
// Created     : Jul 23, 2021
// Last Update : Jul 23, 2021
//============================================================================

#ifndef __IBEX_QIBEXOPTIMIZER_H__
#define __IBEX_QIBEXOPTIMIZER_H__

#include "ibex_Optimizer.h"
#include <utility>
#include <vector>

using namespace std;
namespace ibex {


  class QibexOptimizer : public Optimizer {

  public:

        QibexOptimizer(int n, Ctc& ctc, Bsc& bsc, LoupFinder& finder, CellBufferOptim& buffer,
		       int goal_var, double minwidth, 
			double eps_x=OptimizerConfig::default_eps_x,
			double rel_eps_f=OptimizerConfig::default_rel_eps_f,
			double abs_eps_f=OptimizerConfig::default_abs_eps_f);

        void init();
        /* when using quadratic convex relaxation, contracts the objective and may return a new loup  */

        void qibex_contract(Cell& c);

        /* when using quadratic convex relaxation, calls the external solver
           returns 
           the best point v computed by the relaxation that may be a new loup, and the lower bound*/
    std::pair<Vector,double> qibex_newbounds(IntervalVector & box, int& var_to_bisect, double& ratio, double& gap0);
    int compute_var_to_bisect(const IntervalVector & box,int n_y, const Vector& v, const Vector & w, double& gap0);
    double compute_ratio(const IntervalVector & box,  const Vector& v, int i);
    vector<int> associ;
    vector<int> assocj;
    vector<int> ref_diag_coefs;
    vector<int> hessian_diag_coefs;
    vector<vector<int>> ref_coefs;
    // the width limit for taking into account the variable in the qibex bisection heuristic
    double minwidth;
  };
}
#endif // __IBEX_QIBEXOPTIMIZER_H__
