//============================================================================
//                                  I B E X
// File        : ibex_QibexOptimizer.h
// Author      :  Bertrand Neveu
// Copyright   : IMT Atlantique (France)
// License     : See the LICENSE file
// Created     : Jul 23, 2021
// Last Update : Oct 20, 2022
//============================================================================

#ifndef __IBEX_QIBEXOPTIMIZER_H__
#define __IBEX_QIBEXOPTIMIZER_H__

#include "ibex_Optimizer.h"
#include <utility>
#include <vector>

//using namespace std;
namespace ibex {

  class QibexOptimizer : public Optimizer {

  public:

    QibexOptimizer(int n, Ctc& ctc, Bsc& bsc, LoupFinder& finder, CellBufferOptim& buffer,
		   int goal_var, double minwidth, double tolerance=0.0, 
		   double eps_x=OptimizerConfig::default_eps_x,
		   double rel_eps_f=OptimizerConfig::default_rel_eps_f,
		   double abs_eps_f=OptimizerConfig::default_abs_eps_f);

    void init();
    void contract(Cell & c);
    /* when using quadratic convex relaxation, contracts the objective and may return a new loup  */
    void qibex_contract_and_bound(Cell& c);

    /* the external call through ampl */
    void qibex_relaxation_call(const IntervalVector & box, double objlb);

    /* results of the quadratic convex relaxation 
       when using quadratic convex relaxation,  the external solver  returns 
       the status of the relaxation returned by the solver (solved, solved?, infeasible, unbounded)
       the best point v computed by the relaxation that may be a new loup point,
       the point w of the auxiliary variables (used for the bisection strategy)
       and the lower bound newlb*/
    bool qibex_relaxation_results(std::string & status, double& newlb, Vector & v, Vector& w);

    /* the convex quadratic relaxation   returns a triple made of the best point, the point of the auxiliary variables  and  the new lower bound */
    std::tuple<Vector,Vector,double> qibex_relaxation_analysis (IntervalVector & box, std::string& status);
    /* the analysis of the relaxation gives the next variable to bisect ant its ratio*/
    void qibex_bisection_choice (Cell& c, IntervalVector & qcp_box, Vector&v, Vector& w);
    int compute_var_to_bisect(const IntervalVector & box, const Vector& v, const Vector & w, double& gap0);
    double compute_ratio(const IntervalVector & box,  const Vector& v, int i);
   
    // the width limit for taking into account the variable in the qibex bisection heuristic
    double minwidth;
    // the tolerance on the constraints for a feasible solution ; also used in qibex bisection heuristic
    double tolerance;
    // boolean to use the ibex loupfinder (the qibex loup finder is inside the relaxation)
    bool loupfinderp=true;
    // boolean indicating if the relaxation is made rigourous (by default ; not rigourous)
    bool rigor=false;
    // boolean indicating if the contractors are called after an improvement of the lower bound
    bool recontract=true;
    // redefining of update_loup to take into account the loupfinderp indicator
    bool update_loup(const IntervalVector& box, BoxProperties& prop);
    // checks if v is a new louppoint (after rounding it to integer in case of integer minlp variable) update the loup and louppoint  and returns the new ymax 
    double qibex_loupfinder(Vector& v);
    // redefining of compute_ymax to tka into account the integerobj indicator
    double compute_ymax();
    // redefining of compute_emptybuffer to tka into account the integerobj indicator
    double compute_emptybuffer_uplo();

  private :
    int n_x;
    int n_y;
    std::vector<int> associ;
    std::vector<int> assocj;
    std::vector<int> ref_diag_coefs;
    std::vector<int> hessian_diag_coefs;
    std::vector<std::vector<int>> ref_coefs;
  };
}
#endif // __IBEX_QIBEXOPTIMIZER_H__
