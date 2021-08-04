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
//#include "ibex_NormalizedSystem.h"

namespace ibex {


  class QibexOptimizer : public Optimizer {

  public:

        QibexOptimizer(int n, Ctc& ctc, Bsc& bsc, LoupFinder& finder, CellBufferOptim& buffer,
			int goal_var,
			double eps_x=OptimizerConfig::default_eps_x,
			double rel_eps_f=OptimizerConfig::default_rel_eps_f,
			double abs_eps_f=OptimizerConfig::default_abs_eps_f);


        /* when using quadratic convex relaxation, contracts the objective and may return a new loup  */
  
        void qibex_contract( Interval& y, Cell& c);

        /* when using quadratic convex relaxation, calls the external solver
           returns 
           the best point v computed by the relaxation that may be a new loup, and the lower bound*/
        std::pair<Vector,double> qibex_newbounds(const IntervalVector & box);

  };
}
#endif // __IBEX_QIBEXOPTIMIZER_H__
