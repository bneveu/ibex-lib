//============================================================================
//                                  I B E X
// File        : ibex_LoupFinderDefaultIpopt.cpp
// Author      : Bertrand Neveu
// Copyright   : IMT Atlantique (France)
// License     : See the LICENSE file
// Created     : Aug 03, 2023
// Last update : Aug 03, 2023
//============================================================================

#include "ibex_LoupFinderDefaultIpopt.h"
#include "ibex_LoupFinderInHC4.h"
#include "ibex_LoupFinderFwdBwd.h"
#include "ibex_BxpLinearRelaxArgMin.h"
#include "ibex_LoupFinderProbing.h"
#include "ibex_LoupFinderIpopt.h"
#include "ibex_CellHeap.h"
using namespace std;

namespace ibex {

  LoupFinderDefaultIpopt::LoupFinderDefaultIpopt(System& sys, const System& normsys, const ExtendedSystem& extsys, bool inHC4)  : sys(sys), normsys(normsys), extsys(extsys),									 //	finder_probing(inHC4? (LoupFinder&) *new LoupFinderInHC4(sys) : (LoupFinder&) *new LoupFinderFwdBwd(sys)),
        finder_probing(inHC4? (LoupFinder&) *new LoupFinderInHC4(normsys) : (LoupFinder&) *new LoupFinderProbing(normsys)),
        finder_x_taylor(normsys),
	finder_ipopt(sys,normsys,extsys){

}

  bool LoupFinderDefaultIpopt::integer_check(Vector& pt) {return integer_and_bound_check(normsys,pt);}
  bool LoupFinderDefaultIpopt::is_inner (Vector& pt) {return is_inner0(normsys,pt);}
  double LoupFinderDefaultIpopt::goal_ub (Vector& pt) {return goal_ub0(normsys,pt);}
  void LoupFinderDefaultIpopt::sysbound(Vector& pt) {bound_check(normsys,pt);}
  void LoupFinderDefaultIpopt::sysbound(IntervalVector& vec) {bound_check(normsys,vec);}
  
void LoupFinderDefaultIpopt::add_property(const IntervalVector& init_box, BoxProperties& prop) {
	finder_probing.add_property(init_box,prop);
	finder_x_taylor.add_property(init_box,prop);
	finder_ipopt.add_property(init_box,prop);

	//--------------------------------------------------------------------------
	/* Using line search from LP relaxation minimizer seems not interesting. */
//	if (!prop[BxpLinearRelaxArgMin::get_id(finder_x_taylor.sys)]) {
//		prop.add(new BxpLinearRelaxArgMin(finder_x_taylor.sys));
//	}
	//--------------------------------------------------------------------------

}

std::pair<IntervalVector, double> LoupFinderDefaultIpopt::find(const IntervalVector& box, const IntervalVector& old_loup_point, double old_loup, BoxProperties& prop) {

	pair<IntervalVector,double> p=make_pair(old_loup_point, old_loup);

	bool found=false;

        try { p=finder_ipopt.find(box,p.first,p.second);
	      found=true;
	} catch(NotFound&) { }

	try {
		p=finder_probing.find(box,p.first,p.second,prop);
		found=true;
	} catch(NotFound&) { }

	try {
		// TODO
		// in_x_taylor.set_inactive_ctr(entailed->norm_entailed);
		p=finder_x_taylor.find(box,p.first,p.second,prop);
		found=true;
	} catch(NotFound&) { }

	if (found) {
		//--------------------------------------------------------------------------
		/* Using line search from LP relaxation minimizer seems not interesting. */
		//	BxpLinearRelaxArgMin* argmin=(BxpLinearRelaxArgMin*) prop[BxpLinearRelaxArgMin::get_id(finder_x_taylor.sys)];
		BxpLinearRelaxArgMin* argmin=NULL;
		//--------------------------------------------------------------------------

		if (argmin && argmin->argmin()) {
			Vector loup_point = p.first.lb();
			double loup = p.second;
			LoupFinderProbing(finder_x_taylor.sys).dichotomic_line_search(loup_point, loup, *argmin->argmin(), false);
			//cout << "better loup found! " << loup << endl;
			p=make_pair(loup_point,loup);
		}
		return p;
	} else
		throw NotFound();
}

LoupFinderDefaultIpopt::~LoupFinderDefaultIpopt() {
	delete &finder_probing;
}

} /* namespace ibex */
