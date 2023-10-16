//============================================================================
//                                  I B E X
// File        : ibex_LSmear.h
// Author      : Ignacio Araya
// License     : See the LICENSE file
// Last update : May, 25 2017
//============================================================================

#ifndef __IBEX_QIBEX_LARGEST_FIRST__
#define __IBEX_QIBEX_LARGEST_FIRST__


#include "ibex_OptimLargestFirst.h"
namespace ibex {




/**
 * \ingroup bisector
 *
 * \brief Qibex bisection strategy 
 *
 */
class QibexLargestFirst : public OptimLargestFirst {

public :
	/**
	 * \brief Create the QibexLargestFirst bisector. 
	
	 */
  QibexLargestFirst (int goal_var,bool choose_obj,  double prec, double ratio=Bsc::default_ratio());
	
	/*
	 * \brief Variant with a vector of precisions.
	 *
	 * \see #QibexLargestFirst(System&, double, double, lsmear_mode)
	 */
  QibexLargestFirst (int goal_var,bool choose_obj, const Vector& prec, double ratio=Bsc::default_ratio());
		

	/**
	 * \brief Returns the variable to bisect.
	 */
  BisectionPoint choose_var(const Cell & cell) ;

};

} /* namespace ibex */

#endif /* __IBEX_QIBEXSMEARSUMREL__ */
