//============================================================================
//                                  I B E X
// File        : ibex_QibexSmearSumRelative.h
// Author      : Bertrand Neveu
// License     : See the LICENSE file
// Last update : May, 25 2021
//============================================================================

#ifndef __IBEX_QIBEXSMEARSUMREL__
#define __IBEX_QIBEXSMEARSUMREL__

#include "ibex_SmearFunction.h"

#include "ibex_OptimLargestFirst.h"

#include "ibex_ExtendedSystem.h"
namespace ibex {




/**
 * \ingroup bisector
 *
 * \brief Qibex bisection strategy 
 *
 */
class QibexSmearSumRelative : public SmearSumRelative {

public :
	/**
	 * \brief Create the QibexSmearSumRelative bisector. 
	
	 */
  QibexSmearSumRelative (ExtendedSystem& sys,  double prec, double ratio=Bsc::default_ratio());
	
	/*
	 * \brief Variant with a vector of precisions.
	 *
	 * \see #QibexSmearSumRelative(System&, double, double, lsmear_mode)
	 */
  QibexSmearSumRelative (ExtendedSystem& sys, const Vector& prec, double ratio=Bsc::default_ratio());
		
	/*
	 * \brief Variant with an OptimLargestFirst bisector as default bisector
	 * used by default optimizer
	 */
  QibexSmearSumRelative (ExtendedSystem& sys, double prec, OptimLargestFirst& lf);
/*
	 * \brief Variant with an OptimLargestFirst bisector as default bisector and a vector of precisions

	 */
	QibexSmearSumRelative (ExtendedSystem& sys,const Vector& prec, OptimLargestFirst& lf);
	
	/**
	 * \brief Returns the variable to bisect.
	 */
  BisectionPoint choose_var(const Cell & cell) ;

};
/*============================================ inline implementation ============================================ */

  inline QibexSmearSumRelative::QibexSmearSumRelative(ExtendedSystem& sys,  double prec, double ratio) : SmearSumRelative(sys,prec,ratio)
		 {;}

  inline QibexSmearSumRelative::QibexSmearSumRelative(ExtendedSystem& sys, double prec, OptimLargestFirst& lf) : SmearSumRelative(sys,prec,lf)
  {;}

  inline QibexSmearSumRelative::QibexSmearSumRelative(ExtendedSystem& sys, const Vector& prec, OptimLargestFirst& lf) : SmearSumRelative(sys,prec,lf)
  {;}


  inline QibexSmearSumRelative::QibexSmearSumRelative(ExtendedSystem& sys, const Vector& prec, double ratio) : SmearSumRelative(sys,prec,ratio)
  {;}

} /* namespace ibex */

#endif /* __IBEX_QIBEXSMEARSUMREL__ */
