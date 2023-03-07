//============================================================================
//                                  I B E X
// File        : ibex_QibexSmearSum.h
// Author      : Bertrand Neveu
// License     : See the LICENSE file
// Last update : May, 25 2021
//============================================================================

#ifndef __IBEX_QIBEXSMEARSUM__
#define __IBEX_QIBEXSMEARSUM__

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
class QibexSmearSum : public SmearSum {

public :
	/**
	 * \brief Create the QibexSmearSum bisector. 
	
	 */
  QibexSmearSum (ExtendedSystem& sys,  double prec, double ratio=Bsc::default_ratio());
	
	/*
	 * \brief Variant with a vector of precisions.
	 *
	 * \see #QibexSmearSum(System&, double, double, lsmear_mode)
	 */
  QibexSmearSum (ExtendedSystem& sys, const Vector& prec, double ratio=Bsc::default_ratio());
		
	/*
	 * \brief Variant with an OptimLargestFirst bisector as default bisector
	 * used by default optimizer
	 */
  QibexSmearSum (ExtendedSystem& sys, double prec, OptimLargestFirst& lf);
/*
	 * \brief Variant with an OptimLargestFirst bisector as default bisector and a vector of precisions

	 */
	QibexSmearSum (ExtendedSystem& sys,const Vector& prec, OptimLargestFirst& lf);
	
	/**
	 * \brief Returns the variable to bisect.
	 */
  BisectionPoint choose_var(const Cell & cell) ;

};
/*============================================ inline implementation ============================================ */

  inline QibexSmearSum::QibexSmearSum(ExtendedSystem& sys,  double prec, double ratio) : SmearSum(sys,prec,ratio)
		 {;}

  inline QibexSmearSum::QibexSmearSum(ExtendedSystem& sys, double prec, OptimLargestFirst& lf) : SmearSum(sys,prec,lf)
  {;}

  inline QibexSmearSum::QibexSmearSum(ExtendedSystem& sys, const Vector& prec, OptimLargestFirst& lf) : SmearSum(sys,prec,lf)
  {;}


  inline QibexSmearSum::QibexSmearSum(ExtendedSystem& sys, const Vector& prec, double ratio) : SmearSum(sys,prec,ratio)
  {;}

} /* namespace ibex */

#endif /* __IBEX_QIBEXSMEARSUMREL__ */
