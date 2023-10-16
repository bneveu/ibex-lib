//============================================================================
//                                  I B E X
// File        : ibex_QibexSmearSumRelative.cpp
// Author      : Bertrand Neveu
// License     : See the LICENSE file
// Created     : Sep, 1 2021
// Last Update : May, 3 2021
//============================================================================

#include "ibex_QibexSmearsumRelative.h"


namespace ibex {


  BisectionPoint QibexSmearSumRelative::choose_var(const Cell& cell){
   
    int var = cell.var_to_bisect;
    //       std::cout << " var to bisect " << var << std::endl;
    if (var==-1 || too_small(cell.box,var))
      {//std::cout << " call smear function " << std::endl;
	return SmearFunction::choose_var(cell);
	//	std::cout << " call lf " << std::endl;
	//	return lf->choose_var(cell);
      }
    else
      {//std::cout << " var " << var << std::endl;
	return BisectionPoint(var,cell.ratio,true); // output
      }
  }


} /* namespace ibex */
