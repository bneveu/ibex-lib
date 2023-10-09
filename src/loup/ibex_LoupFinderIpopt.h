//============================================================================
//                                  I B E X
// File        : ibex_LoupFinderIpopt.h
// Author      : Bertrand Neveu
// Copyright   : IMT Atlantique (France)
// License     : See the LICENSE file
// Created     : Aug 2020
// Last Update : Sep 28, 2020
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
 * \brief Calls Ipopt to find a feasible point, which is a local minimum.
 *
 */
class LoupFinderIpopt : public LoupFinder {
public:

	/**
	 * \brief Create the algorithm for a given system.
	 *
	 * \param sys         - The NLP problem.
         * \param normsys     - The normalized NLP problem.
         * \param extsys      - The extended NLP problem
	 */
  LoupFinderIpopt(System& sys, const System& normsys, const ExtendedSystem& extsys);

	/**
	 * \brief Find a new loup in a given box.
	 *
	 * \see comments in LoupFinder.
	 */
  	virtual std::pair<IntervalVector, double> find(const IntervalVector& box, const IntervalVector& loup_point, double loup);


	/**
	 * \brief The NLP problem.
	 */
  System& sys; // original system used for generating .mod file for ampl for ipopt solver ; not const because box is set to current box during the generation  of .mod file.
  
  const System& normsys; // normalized system used for checking the constraints.
  const ExtendedSystem& extsys; // extended system used for the recursing call of optimizer
  Optimizer* optimizer=nullptr; // optimizer data are used for building an optimizer for correcting the point returned by ipopt if it does not verify the constraints
  bool recursive_call=true; // boolean to prevent double recursion of optimizer
  int correction_nodes=0;  // additional nodes for correcting the point given by ipopt
  double correction_time=0.0; // additional time for correcting the point given by ipopt
  
  bool integer_check(Vector& pt);
  bool is_inner(Vector& pt);
  double goal_ub(Vector& pt);
  void sysbound(Vector& pt);
  void sysbound(IntervalVector& vec);  
  std::string ipopt_ampl_run = "/libre/neveu/ampl/ampl model_ipopt.run > amplout.txt";
  void write_system_ampl(const IntervalVector& box); // write the system.mod with box as variable bounds.
  

  void correct_ipopt_sol (Vector&v, double& loup);

  bool ipopt_ampl_results(int n,std::string& status, double& obj, Vector & v);

};

} /* namespace ibex */

#endif /* __IBEX_LOUP_FINDER_PROBING_H__ */
