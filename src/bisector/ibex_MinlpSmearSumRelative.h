//============================================================================
//                                  I B E X                                   
// File        : Smearsumrelative bisector variant for MINLP optimization 
// Author      : Bertrand Neveu
// Copyright   : IMT Atlantique (France)
// License     : See the LICENSE file
// Created     : Aug 19, 2021
// Last Update : Aug 19, 2021
//============================================================================

#ifndef __IBEX_MINLP_SMEARSUMREL_H__
#define __IBEX_MINLP_SMEARSUMREL_H__

#include "ibex_Bsc.h"
#include "ibex_SmearFunction.h"
#include "ibex_System.h"
#include "ibex_LargestFirst.h"

namespace ibex {

/**
 * \ingroup bisector
 *
 * \brief smear sum relative bisector for minlp
 *
 */
class MinlpSmearSumRelative : public SmearFunction {
public:

	/**
	 * \brief Create a bisector with SmearSumRelative heuristic choosing first among the integer variables
	 *
	 */
  MinlpSmearSumRelative(System& sys,  double prec,   LargestFirst& lf, bool gb=true);

	/**
	 * \brief Create a bisector with SmearSumRelative heuristic (choosing first among the integer variables, and if
         * no integer variable could be chosen, among the real variables) see ibex_SmearFunction.h for explanation
         * of this heuristic.
	 *
	 * \param prec             - see #Bsc::Bsc(double).
	 * \param lf : a largest first bisector to be used when the Smear based heuristic could not choose any variable
         * \param gb : boolean indicating if the goal variable can be bisected : default true.
         * TODO . reintroduce param ratio . It is now set to its default value 0.45
	 */

  MinlpSmearSumRelative(System& sys,const Vector& prec,LargestFirst& lf, bool gb=true);

	/**
	 * \brief Return next variable to be bisected.
	 *
	 * called by SmearFunction::choose_var
	 */
  int var_to_bisect(IntervalMatrix& J, const IntervalVector& box) const;
	


	  
};

}// end namespace ibex

#endif // __IBEX_MINLP_SMEARSUMREL_H__
