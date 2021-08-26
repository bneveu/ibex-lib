//============================================================================
//                                  I B E X                                   
// File        : Smearsumr bisector variant for MINLP optimization 
// Author      : Bertrand Neveu
// Copyright   : IMT Atlantique (France)
// License     : See the LICENSE file
// Created     : Aug 19, 2021
// Last Update : Aug 19, 2021
//============================================================================

#ifndef __IBEX_MINLP_SMEARSUM_H__
#define __IBEX_MINLP_SMEARSUM_H__

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
class MinlpSmearSum : public SmearFunction {
public:

	/**
	 * \brief Create a bisector with smearsum heuristic for minlp
	 *
	 */
  MinlpSmearSum(System& sys,  double prec,   LargestFirst& lf, bool gb=true);

	/**
	 * \brief Create a bisector with largest first heuristic.
	 *
	 * \param prec             - see #Bsc::Bsc(double).
	 * \param ratio (optional) - the ratio between the diameters of the left and the right parts of the
	 *                           bisected interval. Default value is 0.45.
	 */
  MinlpSmearSum(System& sys,const Vector& prec,LargestFirst& lf);

	/**
	 * \brief Return next variable to be bisected.
	 *
	 * called by Bsc::bisect(...)
	 */
  int var_to_bisect(IntervalMatrix& J, const IntervalVector& box) const;
	


	  
};

}// end namespace ibex

#endif // __IBEX_MINLP_SMEARSUM_H__
