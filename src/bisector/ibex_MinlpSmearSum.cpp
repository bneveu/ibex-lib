//============================================================================
//                                  I B E X                                   
// File        : ibex_MinlpSmearSum.cpp
// Author      : Bertrand Neveu
// Copyright   : IMT Atlantique (France)
// License     : See the LICENSE file
// Created     : Dec 3, 2018
// Last Update : Oct 17, 2019
//============================================================================
#include "float.h"
#include "ibex_BitSet.h"
#include "ibex_MinlpSmearSum.h"
#include "ibex_NoBisectableVariableException.h"

using namespace std;

namespace ibex {

 
  MinlpSmearSum::MinlpSmearSum(System& sys,  double prec,   LargestFirst& lf, bool gb) : SmearFunction(sys,prec, lf,gb)  {
}
    MinlpSmearSum::MinlpSmearSum(System& sys,  const Vector & prec,   LargestFirst& lf, bool gb) : SmearFunction(sys,prec, lf,gb)  {
}




   int MinlpSmearSum::var_to_bisect(IntervalMatrix& J, const IntervalVector& box) const {
    double max_magn = NEG_INFINITY;
    int var = -1;
    BitSet& b= *(sys.get_integer_variables());
    //    cout << "integer variables "  << b << nbvars << " goal_to_bisect " << goal_to_bisect << endl;
    for (int j=0; j<nbvars; j++) {
      if ((!too_small(box,j))&& (goal_to_bisect || (j!= goal_var())&& b[j])) { 
	double sum_smear=0;
	for (int i=0; i<sys.f_ctrs.image_dim(); i++) {
	  if (constraint_to_consider (i, box))
	    sum_smear+= J[i][j].mag() *box[j].diam();
	}
	if (sum_smear > max_magn) {
	  max_magn = sum_smear;
	  var = j;
	}
      }
    }
    if (var==-1)  // no integer variable was chosen
      {
	max_magn = NEG_INFINITY;
	for (int j=0; j<nbvars; j++) {
	  if ((!too_small(box,j))&& (goal_to_bisect || j!= goal_var())) { // && (box[j].mag() <1 ||  box[j].diam()/ box[j].mag() >= prec(j))) {
	    double sum_smear=0;
	    for (int i=0; i<sys.f_ctrs.image_dim(); i++) {
	      if (constraint_to_consider (i, box))
		sum_smear+= J[i][j].mag() *box[j].diam();
	    }
	    if (sum_smear > max_magn) {
	      max_magn = sum_smear;
	      var = j;
	    }
	  }
	}
	  
      }

    //    cout << " var " << var << endl;
    return var;
   }




 
} // end namespace ibex
