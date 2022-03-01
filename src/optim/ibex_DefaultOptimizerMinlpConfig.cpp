//============================================================================
//                                  I B E X
// File        : ibex_DefaultOptimizerMinlpConfig.cpp
// Author      : Gilles Chabert
// Copyright   : IMT Atlantique (France)
// License     : See the LICENSE file
// Created     : Dec 11, 2014
// Last Update : Oct 14, 2019
//============================================================================

#include "ibex_DefaultOptimizerMinlpConfig.h"

#include "ibex_CtcHC4.h"
#include "ibex_CtcInteger.h"
#include "ibex_CtcAcid.h"
#include "ibex_CtcCompo.h"
#include "ibex_CtcFixPoint.h"
#include "ibex_CtcLinearRelax.h"
#include "ibex_CellDoubleHeap.h"
#include "ibex_SmearFunction.h"
#include "ibex_MinlpSmearSumRelative.h"
#include "ibex_OptimLargestFirst.h"
#include "ibex_LoupFinderDefault.h"
#include "ibex_LoupFinderCertify.h"
#include "ibex_Array.h"
#include "ibex_Random.h"
#include "ibex_CellBeamSearch.h"
#include "ibex_CellHeap.h"
#include "ibex_CtcKuhnTucker.h"
#include "ibex_CtcKuhnTuckerLP.h"
#include "ibex_BitSet.h"

using namespace std;

namespace ibex {

namespace {

enum { 	NORMALIZED_SYSTEM_TAG,
		EXTENDED_SYSTEM_TAG,
		CTC_TAG,
		BSC_TAG,
		CELL_BUFFER_TAG,
		LOUP_FINDER_TAG };
}

DefaultOptimizerMinlpConfig::DefaultOptimizerMinlpConfig(const System& sys) : sys(sys), kkt(false) {
	set_eps_h(ExtendedSystem::default_eps_h);
	set_rigor(default_rigor);
	set_inHC4(default_inHC4);
	// by defaut, we apply KKT for unconstrained problems
	//	set_kkt(sys.nb_ctr==0);
	set_random_seed(default_random_seed);
}

// note:deprecated.
DefaultOptimizerMinlpConfig::DefaultOptimizerMinlpConfig(const System& sys, double rel_eps_f, double abs_eps_f,
							double eps_h, bool rigor, bool inHC4, bool kkt,
							int random_seed, double eps_x) : sys(sys) {

	set_rel_eps_f(rel_eps_f);
	set_abs_eps_f(abs_eps_f);
	set_eps_h(eps_h);
	set_rigor(rigor);
	set_inHC4(inHC4);
	set_kkt(kkt);
	set_random_seed(random_seed);
	set_eps_x(eps_x);
}

DefaultOptimizerMinlpConfig::~DefaultOptimizerMinlpConfig() {

}

void DefaultOptimizerMinlpConfig::set_eps_h(double _eps_h) {
	eps_h = _eps_h;
}

void DefaultOptimizerMinlpConfig::set_rigor(bool _rigor) {
	rigor = _rigor;

	if (!rigor && kkt) {
		for (int i=0; i<sys.nb_ctr; i++)
			if (sys.ctrs[i].op==EQ) {
				kkt=false;
				ibex_warning("[OptimizerMinlpConfig] KKT automatically disabled if rigor mode is switched off with equalities.");
				return;
			}
	}
}

void DefaultOptimizerMinlpConfig::set_inHC4(bool _inHC4) {
	// check if inHC4 can be applied if it
	// has not already been done.
	if (!inHC4 && _inHC4 && sys.nb_ctr>0) {
		// *******
		// Warning: generates components of f_ctrs!!
		//          (but LoupFinderInHC4 does anyway)
		// *******
		for (int i=0; i<sys.f_ctrs.image_dim(); i++) {
			if (!sys.f_ctrs[i].inhc4revise().implemented()) {
				return;
			}
		}
	}
	inHC4 = _inHC4;
}

void DefaultOptimizerMinlpConfig::set_kkt(bool _kkt) {
	kkt = _kkt;

	// if KKT is applied with equalities, rigor mode is forced.
	if (kkt && !rigor && sys.nb_ctr>1) {
		rigor=true;
		// mandatory with equalities and strongly preferable with inequalities only
		ibex_warning("[OptimizerMinlpConfig] Rigor mode automatically activated with KKT.");
	}
	
}

 

void DefaultOptimizerMinlpConfig::set_random_seed(int _random_seed) {
	random_seed = _random_seed;
	RNG::srand(random_seed);
}

  void DefaultOptimizerMinlpConfig::set_integer_variables(vector<int>& l){
    ExtendedSystem& ext_sys=get_ext_sys();
    BitSet b (ext_sys.nb_var);
    for (int i=0; i<l.size();i++)
      b.add(l[i]);
    ext_sys.set_integer_variables(b);
    ext_sys.minlp=true;
    cout << " integer variables " << ext_sys.get_integer_variables() << endl;
    NormalizedSystem& norm_sys=get_norm_sys();
    norm_sys.set_integer_variables(b);
    norm_sys.minlp=true;
  }

// The two next functions are necessary because we need
// the normalized and extended system to build
// arguments of the base class constructor (ctc, bsc, loup finder, etc.)
// and we don't know which argument is evaluated first

NormalizedSystem& DefaultOptimizerMinlpConfig::get_norm_sys() {
	if (found(NORMALIZED_SYSTEM_TAG)) {
		return get<NormalizedSystem>(NORMALIZED_SYSTEM_TAG);
	} else {
		return rec(new NormalizedSystem(sys,eps_h), NORMALIZED_SYSTEM_TAG);
	}
}

ExtendedSystem& DefaultOptimizerMinlpConfig::get_ext_sys() {
	if (found(EXTENDED_SYSTEM_TAG)) {
		return get<ExtendedSystem>(EXTENDED_SYSTEM_TAG);
	} else {
		return rec(new ExtendedSystem(sys,eps_h), EXTENDED_SYSTEM_TAG);
	}
}

unsigned int DefaultOptimizerMinlpConfig::nb_var() {
	return sys.nb_var;
}

Ctc& DefaultOptimizerMinlpConfig::get_ctc() {
	if (found(CTC_TAG)) // in practice, get_ctc() is only called once by Optimizer.
		return get<Ctc>(CTC_TAG);

	const ExtendedSystem& ext_sys = get_ext_sys();

	Array<Ctc> ctc_list( 3);
	
        cout << "minlp " << ext_sys.minlp << endl;
	cout << "integer_variables " << ext_sys.get_integer_variables() << endl;
	// first contractor on ext_sys : incremental HC4 (propag ratio=0.01)
	ctc_list.set_ref(0, rec(new CtcCompo (
						rec (new CtcInteger (ext_sys.nb_var, ext_sys.get_integer_variables())),
						rec (new CtcHC4 (ext_sys,0.01,true)),
						rec (new CtcInteger (ext_sys.nb_var, ext_sys.get_integer_variables())))));
	// second contractor on ext_sys : "Acid" with incremental HC4 (propag ratio=0.1)
        ctc_list.set_ref(1, rec(new CtcCompo (rec (new CtcAcid (ext_sys,rec(new CtcHC4 (ext_sys,0.1,true)),true)),
					      rec(new  CtcInteger (ext_sys.nb_var, ext_sys.get_integer_variables())))));

				;
	// the last contractor is "XNewton"

	if (ext_sys.nb_ctr > 1) {
	  ctc_list.set_ref(2,rec(new CtcFixPoint
				(rec(new CtcCompo(
						rec (new CtcInteger (ext_sys.nb_var, ext_sys.get_integer_variables())),

						rec(new CtcLinearRelax(ext_sys)),
						rec (new CtcInteger (ext_sys.nb_var, ext_sys.get_integer_variables())),
						rec(new CtcHC4(ext_sys,0.1)),
						rec (new CtcInteger (ext_sys.nb_var, ext_sys.get_integer_variables())))),
				 default_relax_ratio)));
	} else {
		ctc_list.set_ref(2,rec(new CtcLinearRelax(ext_sys)));
	}
	return rec(new CtcCompo(ctc_list), CTC_TAG);
}


Bsc& DefaultOptimizerMinlpConfig::get_bsc() {
	if (found(BSC_TAG)) // in practice, get_bsc() is only called once by Optimizer.
			return get<Bsc>(BSC_TAG);

	return rec(new MinlpSmearSumRelative(
			get_ext_sys(),eps_x,
			rec(new OptimLargestFirst(get_ext_sys().goal_var(),false,eps_x,default_bisect_ratio))),
			BSC_TAG);
}

LoupFinder& DefaultOptimizerMinlpConfig::get_loup_finder() {
	if (found(LOUP_FINDER_TAG)) // in practice, get_loup_finder() is only called once by Optimizer.
			return get<LoupFinder>(LOUP_FINDER_TAG);

	const NormalizedSystem& norm_sys = get_norm_sys();

	return rec(rigor? (LoupFinder*) new LoupFinderCertify(sys,rec(new LoupFinderDefault(norm_sys, inHC4))) :
			(LoupFinder*) new LoupFinderDefault(norm_sys, inHC4), LOUP_FINDER_TAG);
}

CellBufferOptim& DefaultOptimizerMinlpConfig::get_cell_buffer() {
	if (found(CELL_BUFFER_TAG)) // in practice, get_cell_buffer() is only called once by Optimizer.
			return get<CellBufferOptim>(CELL_BUFFER_TAG);

	const ExtendedSystem& ext_sys = get_ext_sys();

//			  (CellBufferOptim&) rec(new CellDoubleHeap(get_ext_sys()))
	return (CellBufferOptim&) rec(new CellHeap(get_ext_sys()));
	/*
	return (CellBufferOptim&) rec (new  CellBeamSearch (
			(CellHeap&) rec (new CellHeap (ext_sys)),
			(CellHeap&) rec (new CellHeap (ext_sys)),
			ext_sys), CELL_BUFFER_TAG);
	*/
}

int DefaultOptimizerMinlpConfig::goal_var() {
	return get_ext_sys().goal_var();
}

} /* namespace ibex */
