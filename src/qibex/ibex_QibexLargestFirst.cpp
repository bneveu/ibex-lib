//============================================================================
//                                  I B E X                                   
// File        : ibex_OptimLargestFirst.cpp
// Author      : Bertrand Neveu, Gilles Chabert
// Copyright   : IMT Atlantique (France)
// License     : See the LICENSE file
// Created     : Dec 3, 2018
// Last Update : Apr 9, 2021
//============================================================================
#include "float.h"
#include "ibex_QibexLargestFirst.h"
#include "ibex_NoBisectableVariableException.h"

using namespace std;

namespace ibex {


  QibexLargestFirst::QibexLargestFirst(int goal_var,bool choose_obj, double prec,  double ratio) : OptimLargestFirst(goal_var, choose_obj, prec,ratio)   {
}

  QibexLargestFirst::QibexLargestFirst(int goal_var,bool choose_obj,const Vector& prec,double ratio) :OptimLargestFirst(goal_var, choose_obj, prec,ratio) {

}

BisectionPoint QibexLargestFirst::choose_var(const Cell& cell){
   
    int var = cell.var_to_bisect;
    //       std::cout << " var to bisect " << var << std::endl;
    if (var==-1 || too_small(cell.box,var))
      {
	return OptimLargestFirst::choose_var(cell);
      }
    else
      {//std::cout << " var " << var << std::endl;
	return BisectionPoint(var,cell.ratio,true); // output
	//	return BisectionPoint(var,ratio,true); // output
      }
  }


  

} // end namespace ibex
