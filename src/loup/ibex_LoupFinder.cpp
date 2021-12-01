//============================================================================
//                                  I B E X
// File        : ibex_LoupFinder.cpp
// Author      : Gilles Chabert, Ignacio Araya, Bertrand Neveu
// Copyright   : IMT Atlantique (France)
// License     : See the LICENSE file
// Created     : July 09, 2017
//============================================================================

#include "ibex_LoupFinder.h"

using namespace std;

namespace ibex {

LoupFinder::~LoupFinder() {

}

void LoupFinder::add_property(const IntervalVector& init_box, BoxProperties& prop) {

}
  
  /*
  bool LoupFinder::point_check(const System& sys, Vector & pt){
    return    sys.is_inner(pt);
  }
  */
  
bool LoupFinder::integer_and_bound_check(const System& sys, Vector & pt){
  if (sys.minlp){
      //      double eps=1.e-2;
      double eps=0.5;
      //      cout << " pt " << pt << endl;
      BitSet b = sys.get_integer_variables();
      for (int i=0; i< pt.size(); i++){
	//	cout << i << " " << integer_variables[i] << endl;
	if (b[i])
	  {Interval intvec =integer(Interval(pt[i]-eps,pt[i]+eps));
	    if (intvec.is_empty() || intvec.diam() >=1)
	      return false;
	    else{
	      //     cout << i << "  " << pt[i] << endl;
	      pt[i]= intvec.mid();
	      if (pt[i] < sys.box[i].lb() || pt[i] > sys.box[i].ub())
		return false;
	    }
	  }
	else bound_check_i(sys,pt,i);
      }
    }

    else bound_check(sys,pt);

    return true;

}
  void LoupFinder::bound_check(const System& sys,Vector & pt){
    for (int i=0; i< pt.size() ; i++)
      bound_check_i(sys,pt,i);
  }
    
  
  void LoupFinder::bound_check_i(const System& sys,Vector & pt, int i){
    if (pt[i] < sys.box[i].lb())
	  pt[i]=sys.box[i].lb();
      if (pt[i] > sys.box[i].ub())
	pt[i]=sys.box[i].ub();
  }

  void LoupFinder::bound_check(const System& sys,IntervalVector & vec){
    for (int i=0; i< vec.size() ; i++)
      bound_check_i(sys,vec,i);
  }

  void LoupFinder::bound_check_i(const System& sys,IntervalVector & vec,int i){
    if (vec[i].lb() < sys.box[i].lb())
      vec[i]=Interval(sys.box[i].lb(),vec[i].ub());
    if (vec[i].ub() > sys.box[i].ub())
      vec[i]=Interval(vec[i].lb(),sys.box[i].ub());
  }		      

  
bool LoupFinder::check(const System& sys,  Vector& pt, double& loup, bool _is_inner) {
	// "res" will contain an upper bound of the criterion
	double res = sys.goal_ub(pt);

	// check if f(x) is below the "loup" (the current upper bound).
	//
	// The "loup" and the corresponding "loup_point" (the current minimizer)
	// will be updated if the constraints are satisfied.
	// The test of the constraints is done only when the evaluation of the criterion
	// is better than the loup (a cheaper test).

	if (res<loup) {
	  if (sys.minlp){
	  //	  cout << " loup " << loup << " res " << res << " pt " << pt << endl;
	    if (integer_and_bound_check(sys,pt) && sys.is_inner(pt)) {
	      res = sys.goal_ub(pt); // integer_and_bound_check may modify the loup_point by making integer the values of integer variables
	      // and putting it in the initial box
	      // we have to recompute the criterion.
	      if (res<loup) {
		loup = res;
		return true;}
	      else return false;
	    }
	    return false;
	  }
	  else
	    if (_is_inner || sys.is_inner(pt)){
	      loup = res;
	      return true;
	    }
	  return false;
	}
	return false;
}

  bool LoupFinder::is_inner0(const System& sys, Vector& pt) {return sys.is_inner(pt);}
  double LoupFinder::goal_ub0(const System& sys, Vector& pt) {return sys.goal_ub(pt);}
  
  void LoupFinder::monotonicity_analysis(const System& sys, IntervalVector& box, bool is_inner) {

	size_t n=sys.nb_var;

	if (!is_inner && sys.f_ctrs.used_vars.size()==n)
		// if there is no inner box and all the variables appear
		// in the constraints, nothing can be done
		return;

	IntervalVector g(n);
	sys.goal->gradient(box,g);

	for (size_t j=0; j<n; j++) {
		if (is_inner || !sys.f_ctrs.used(j)) {
			if (g[j].lb()>=0 && box[j].lb()!=NEG_INFINITY) box[j]=box[j].lb();
			if (g[j].ub()<=0 && box[j].ub()!=POS_INFINITY) box[j]=box[j].ub();
		}
	}
}

} // namespace ibex
