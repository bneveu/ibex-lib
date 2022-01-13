//============================================================================
//                                  I B E X                                   
// File        : ibex_MinlpLargestFirst.cpp
// Author      : Bertrand Neveu
// Copyright   : IMT Atlantique (France)
// License     : See the LICENSE file
// Created     : Dec 3, 2018
// Last Update : Oct 17, 2019
//============================================================================
#include "float.h"
#include "ibex_BitSet.h"
#include "ibex_MinlpSmearSumRelative.h"
#include "ibex_NoBisectableVariableException.h"

using namespace std;

namespace ibex {

 
  MinlpSmearSumRelative::MinlpSmearSumRelative(System& sys,  double prec,   LargestFirst& lf, bool gb) : SmearFunction(sys,prec, lf,gb)  {
}

  MinlpSmearSumRelative::MinlpSmearSumRelative(System& sys,const Vector& prec,LargestFirst& lf) : SmearFunction (sys,prec, lf)  {
}



   int MinlpSmearSumRelative::var_to_bisect(IntervalMatrix& J, const IntervalVector& box) const {
    double max_magn = NEG_INFINITY;
    int var = -1;
    double* ctrjsum = new double[sys.f_ctrs.image_dim()];
    BitSet b= sys.get_integer_variables();
    for (int i=0; i<sys.f_ctrs.image_dim(); i++) {
      ctrjsum[i]=0;
      // not an extended system or constraint is active or it is the objective 
      if (constraint_to_consider(i, box))
	for (int j=0; j<nbvars ; j++) {
	  if (b[j] ||  (goal_to_bisect && j== goal_var())){

	    ctrjsum[i]+= J[i][j].mag() * box[j].diam();
	    //	    cout << " j " << j << " " << ctrjsum[i] << endl;
	  }
	}
      //      cout << " i " << ctrjsum[i] << endl;
    }
    // computes the variable with the maximal sum of normalized impacts
    for (int j=0; j<nbvars; j++) {
      if (!too_small(box,j) && (goal_to_bisect || j!= goal_var()) && b[j]==1){
	double sum_smear=0;
	for (int i=0; i<sys.f_ctrs.image_dim(); i++) {
	  if (ctrjsum[i]!=0)
	    sum_smear+= J[i][j].mag() * box[j].diam() / ctrjsum[i];
	}

	//        cout << " j " << j << " " << sum_smear << " " << max_magn << endl;	  
	if (sum_smear > max_magn) {
	  max_magn = sum_smear;
	  var = j;
	}
      }
    }

    if (var==-1)
      {
	double max_magn = NEG_INFINITY;
	for (int i=0; i<sys.f_ctrs.image_dim(); i++) {
	  ctrjsum[i]=0;
      // not an extended system or constraint is active or it is the objective 
	  if (constraint_to_consider(i, box))
	    for (int j=0; j<nbvars ; j++) {
	      ctrjsum[i]+= J[i][j].mag() * box[j].diam();
	    }
	}
	
	for (int j=0; j<nbvars; j++) {
	  if (!too_small(box,j) && (goal_to_bisect || j!= goal_var())){ //&&  (box[j].mag() <1 ||  box[j].diam()/ box[j].mag() >= prec(j))) {
	    double sum_smear=0;
	    for (int i=0; i<sys.f_ctrs.image_dim(); i++) {
	      if (ctrjsum[i]!=0)
		sum_smear+= J[i][j].mag() * box[j].diam() / ctrjsum[i];
	    }

	  
	    if (sum_smear > max_magn) {
	      max_magn = sum_smear;
	      var = j;
	    }
	  }
	}	  
      }

    delete[] ctrjsum;
    //    cout << " var " << var << endl;
    return var;
   }




 
} // end namespace ibex
