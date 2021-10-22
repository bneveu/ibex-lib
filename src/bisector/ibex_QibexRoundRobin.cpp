//============================================================================
//                                  I B E X                                   
// File        : ibex_RoundRobin.cpp
// Author      : Gilles Chabert
// Copyright   : IMT Atlantique (France)
// License     : See the LICENSE file
// Created     : May 8, 2012
// Last Update : Jun 29, 2018
//============================================================================

#include "ibex_QibexRoundRobin.h"
#include "ibex_NoBisectableVariableException.h"

using namespace std;

namespace ibex {

QibexRoundRobin::QibexRoundRobin(double prec, double ratio) : Bsc(prec), ratio(ratio) {

}

QibexRoundRobin::QibexRoundRobin(const Vector& prec, double ratio) : Bsc(prec), ratio(ratio) {

}

BisectionPoint QibexRoundRobin::choose_var(const Cell& cell) {

  int proposed_var=cell.var_to_bisect;

	const IntervalVector& box=cell.box;

	int n = box.size();

	//	cout << "proposed var " << proposed_var << endl;
	// if no proposed var, choice becomes roundrobin  ??  seems random
	if (proposed_var == -1 || too_small(box,proposed_var))
	    proposed_var = rand()%n;

	int var = proposed_var;
	

	while (var != (proposed_var-1)%n && too_small(box,var))
		var = (var + 1)%n;
	// if no variable can be bisected an exception is thrown
	if (var==(proposed_var-1)%n && too_small(box,var))
		throw NoBisectableVariableException();

	else  {// cout << " var " << var << endl;
	  
		return BisectionPoint(var,ratio,true); // output
	}
}

} // end namespace ibex
