//============================================================================
//                                  I B E X                                   
// File        : ibex_CellStackOptim.cpp
// Author      : Gilles Chabert
// Copyright   : Ecole des Mines de Nantes (France)
// License     : See the LICENSE file
// Created     : May 12, 2012
// Last Update : May 12, 2012
//============================================================================

#include "ibex_CellStackOptim.h"
#include <cfloat>

using namespace std;


namespace ibex {

  CellStackOptim::CellStackOptim(const ExtendedSystem & sys): CellBufferOptim(),sys(sys){}
  
void CellStackOptim::flush() {
	while (!cstack.empty()) {
		delete cstack.top();
		cstack.pop();
	}
	costs.clear();
}

unsigned int CellStackOptim::size() const {
	return cstack.size();
}

bool CellStackOptim::empty() const {
	return cstack.empty();
}

void CellStackOptim::push(Cell* cell) {
	if (capacity>0 && size()==capacity) throw CellBufferOverflow();
	cstack.push(cell);
	int goal_var=sys.goal_var();
	costs.push_back(cell->box[goal_var].lb());
}

Cell* CellStackOptim::pop() {
	Cell* c = cstack.top();
	cstack.pop();
	int goal_var=sys.goal_var();
	suppress(c->box[goal_var].lb());
	return c;
}

  void CellStackOptim::suppress(double val){
    list<double>::iterator it;
    for (it=costs.begin() ; it != costs.end(); it++) 
      {if (*it==val) {
	  costs.erase (it); break;
	}
      }
  }

	  
Cell* CellStackOptim::top() const {
	return cstack.top();
}

  void CellStackOptim::contract(double loup) {;}

  double CellStackOptim::minimum() const{
    double minimum = DBL_MAX;
    list<double>::const_iterator it;
    for (it=costs.begin() ; it != costs.end(); it++)
      {if (*it < minimum)
	  minimum=*it;
      }
    return minimum;
  }

  

} // end namespace ibex
