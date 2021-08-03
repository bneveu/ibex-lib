//============================================================================
//                                  I B E X                                   
// File        : Largest First bisector variant for MINLP optimization 
// Author      : Bertrand Neveu
// Copyright   : IMT Atlantique (France)
// License     : See the LICENSE file
// Created     : Dec 3, 2018
// Last Update : Oct 17, 2019
//============================================================================

#ifndef __IBEX_MINLP_FIRST_H__
#define __IBEX_MINLP_LARGEST_FIRST_H__

#include "ibex_Bsc.h"
#include "ibex_OptimLargestFirst.h"
#include "ibex_System.h"

namespace ibex {

/**
 * \ingroup bisector
 *
 * \brief largest-first bisector.
 *
 */
class MinlpLargestFirst : public OptimLargestFirst {
public:

	/**
	 * \brief Create a bisector with largest-first heuristic.
	 *
	 * \param prec             - see #Bsc::Bsc(double). By default, 0 which means an endless uniform bisection process.
	 * \param ratio (optional) - the ratio between the diameters of the left and the right parts of the
	 *                           bisected interval. Default value is 0.45.
	 */
  MinlpLargestFirst( System& sys, int goal_var,bool choose_obj,double prec=0, double ratio=Bsc::default_ratio());

	/**
	 * \brief Create a bisector with largest first heuristic.
	 *
	 * \param prec             - see #Bsc::Bsc(double).
	 * \param ratio (optional) - the ratio between the diameters of the left and the right parts of the
	 *                           bisected interval. Default value is 0.45.
	 */
  MinlpLargestFirst( System& sys, int goal_var,bool choose_obj,const Vector& prec, double ratio=Bsc::default_ratio());

	/**
	 * \brief Return next variable to be bisected.
	 *
	 * called by Bsc::bisect(...)
	 */
	virtual BisectionPoint choose_var(const Cell& cell);


	

        


 protected :
    System& sys;
/** 
	  * \brief  Condition for not bisecting the variable i
	  */
//	bool nobisectable(const IntervalVector& box, int i) const;

  //private :
  //	double max_diam_nobj;
  //        int max_diam_var;
	  
};

} // end namespace ibex

#endif // __IBEX_OPTIM_LARGEST_FIRST_H__
