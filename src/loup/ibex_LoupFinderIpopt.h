//============================================================================
//                                  I B E X
// File        : ibex_LoupFinderProbing.h
// Author      : Gilles Chabert
// Copyright   : IMT Atlantique (France)
// License     : See the LICENSE file
// Created     : May 14, 2012
// Last Update : Jul 09, 2017
//============================================================================

#ifndef __IBEX_LOUP_FINDER_IPOPT_H__
#define __IBEX_LOUP_FINDER_IPOPT_H__

#include "ibex_LoupFinder.h"
#include "ibex_System.h"
#include "ibex_Optimizer.h"
#include "ibex_CellHeap.h"
namespace ibex {

/**
 * \ingroup optim
 *
 * \brief 
 *
 */
class LoupFinderIpopt : public LoupFinder {
public:

	/**
	 * \brief Create the algorithm for a given system.
	 *
	 * \param sys         - The NLP problem.
	 * \param sample_size - Number of sample points inside the box.
	 */
  LoupFinderIpopt(System& sys, const System& normsys, const ExtendedSystem& extsys);

	/**
	 * \brief Find a new loup in a given box.
	 *
	 * \see comments in LoupFinder.
	 */
  	virtual std::pair<IntervalVector, double> find(const IntervalVector& box, const IntervalVector& loup_point, double loup);


	/**
	 * \brief Second method (line probing).
	 *
	 * Performs a dichotomic search between the current loup-point and its projection on the
	 * facet of the input box in the opposite direction of its gradient.
	 *
	 * return true if the loup has been modified.
	 */

	/**
	 * \brief The NLP problem.
	 */
  System& sys; // original system used for generating .mod file for ampl for ipopt solver ; not const before box is set to current box during the generation  of .mod file. 
  const System& normsys; // normalized system used for the checking the constraints.
  const ExtendedSystem& extsys; // extended system used for the recursing call of optimizer
  Optimizer* optimizer=nullptr;
  bool recursive_call=true;
  int correction_nodes=0;
  double correction_time=0.0;
  void correct_ipopt_sol (Vector&v, double& loup);
  void write_system_ampl(const IntervalVector& box);
  bool ipopt_ampl_results(int n,std::string& status, double& obj, Vector & v);
  bool integer_check(Vector& pt);
  bool is_inner(Vector& pt);
  double goal_ub(Vector& pt);
  void sysbound(Vector& pt);
  void sysbound(IntervalVector& vec);  

 
protected:
//
//	/**
//	 * Current loup-point
//	 */
//	Vector loup_point;
//
//	/**
//	 * Current loup
//	 */
//	double loup;
};

} /* namespace ibex */

#endif /* __IBEX_LOUP_FINDER_PROBING_H__ */
