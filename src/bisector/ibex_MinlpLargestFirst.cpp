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
#include "ibex_MinlpLargestFirst.h"
#include "ibex_NoBisectableVariableException.h"

using namespace std;

namespace ibex {

  // to not bisect very big objectives 
  double objectivebisect_ratiolimit0=1.e10;

  MinlpLargestFirst::MinlpLargestFirst(System& sys, int goal_var,bool choose_obj, double prec,  double ratio) : OptimLargestFirst(goal_var,choose_obj,prec, ratio), sys(sys)  {
}

  MinlpLargestFirst::MinlpLargestFirst(System& sys,int goal_var, bool choose_obj,const Vector& prec,double ratio) :OptimLargestFirst(goal_var,choose_obj,prec, ratio), sys(sys) {

}


BisectionPoint MinlpLargestFirst::choose_var(const Cell& cell) {
        const IntervalVector& box=cell.box;
	int var =-1;
	double l=0.0;
	BitSet& b= *(sys.get_integer_variables());
	//	cout << " b" << b << endl ;
	for (int i=0; i< box.size(); i++){

	  if (i!= goal_var && b[i]){
	    if ( ! nobisectable (box,i)){
	      if (var==-1) {
		var=i;
		l = uniform_prec()? box[i].diam() : (box[i].diam()/prec(i));
	      }
	      else {
		double l_tmp = uniform_prec()? box[i].diam() : (box[i].diam()/prec(i));
		if (l_tmp>l) {
		  var = i;
		  l = l_tmp;
		}
	      }
	    }
	  }
	}
	//	cout << " bisected var " << var  << " l " << l << endl ;
	if (var !=-1){
	  return BisectionPoint(var,ratio,true);
	}
	else{
	  for (int i=0; i< box.size(); i++){

	  if (i!= goal_var && !b[i]){
	    if ( ! nobisectable (box,i)){
	      if (var==-1) {
		var=i;
		l = uniform_prec()? box[i].diam() : (box[i].diam()/prec(i));
	      }
	      else {
		double l_tmp = uniform_prec()? box[i].diam() : (box[i].diam()/prec(i));
		if (l_tmp>l) {
		  var = i;
		  l = l_tmp;
		}
	      }
	    }
	  }
	}
	}
	if ((choose_obj == true)
	      &&  !(nobisectable (box,goal_var))
	      && (l < box[goal_var].diam())
	      && box[goal_var].diam()/l < objectivebisect_ratiolimit0)
	  var=goal_var;
	//	cout << " bisected var " << var  << " l " << l << endl ;
	if (var !=-1){
	  return BisectionPoint(var,ratio,true);
	}
	else {
	  throw NoBisectableVariableException();
	}
	
}
  
 
} // end namespace ibex
