//============================================================================
//                                  I B E X
// File        : ibex_MinlpLSmear.h
// Author      : Ignacio Araya Bertrand Neveu
// License     : See the LICENSE file
// Last update : April, 26 2023
//============================================================================

#ifndef __IBEX_MINLPLSMEAR__
#define __IBEX_MINLPLSMEAR__

#include "ibex_SmearFunction.h"
#include "ibex_MinlpSmearSumRelative.h"
#include "ibex_LargestFirst.h"
#include "ibex_LPSolver.h"
#include "ibex_ExtendedSystem.h"
namespace ibex {

typedef enum{ MINLPLSMEAR=0, MINLPLSMEAR_MG } minlplsmear_mode;


/**
 * \ingroup bisector
 *
 * \brief bisector which first ponderates the constraints by using the dual
 * solution of a linear programming relaxation and then computes the impact
 * by using the Smear function heuristic.
 *
 */
class MinlpLSmear : public MinlpSmearSumRelative {

public :
	/**
	 * \brief Create the MinlpLSmear bisector. See the article:
	 * Araya, I., Neveu, B. lsmear: a variable selection strategy for interval
	 * branch and bound solvers (2017)
	 *
	 * \param lsmode	- variant (see the paper)
	 *
	 * For the parameters, see #SmearFunction::SmearFunction(ExtendedSystem&, double, double).

	MinlpLSmear (ExtendedSystem& sys,  double prec, double ratio=Bsc::default_ratio(),
			minlplsmear_mode lsmode=MINLPLSMEAR_MG);
	*/
	/*
	 * \brief Variant with a vector of precisions.
	 *
	 * \see #MinlpLSmear(System&, double, double, lsmear_mode)

  	MinlpLSmear (ExtendedSystem& sys, const Vector& prec, double ratio=Bsc::default_ratio(),
  			minlplsmear_mode lsmode=MINLPLSMEAR_MG);
	*/
	/*
	 * \brief Variant with an LargestFirst bisector as default bisector
	 * used by default optimizer
	 */
  MinlpLSmear (ExtendedSystem& sys, double prec, LargestFirst& lf,bool gb, minlplsmear_mode lsmode=MINLPLSMEAR_MG);
/*
	 * \brief Variant with an OptimLargestFirst bisector as default bisector and a vector of precisions

	 */
  MinlpLSmear (ExtendedSystem& sys,const Vector& prec, LargestFirst& lf, bool gb, minlplsmear_mode lsmode=MINLPLSMEAR_MG);
	/**
	 * \brief Delete this.
	 */
	~MinlpLSmear();

	/**
	 * \brief Returns the variable to bisect.
	 */
	virtual int var_to_bisect(IntervalMatrix& J,const IntervalVector& box) const;

	/**
	 * \brief Computes the dual solution of the linear program mid(J).x<=0
	 *
	 * Returns OPTIMAL if the dual solution was computed successfully and UNKNOWN otherwise.
	 *
	 * \param J 	- the jacobian matrix
	 * \param x 	- the current box
	 * \param dual 	- the dual solution that will be returned
	 */
	LPSolver::Status getdual(IntervalMatrix& J,const IntervalVector& x, Vector& dual) const;

	/**
	 * \brief The linear solver
	 */
	LPSolver* mylinearsolver;

	/**
	 * \brief The lsmear variant (MINLPLSMEAR or MINLPLSMEAR_MG)
	 */
	minlplsmear_mode lsmode;


};

/*============================================ inline implementation ============================================ */
/*
inline MinlpLSmear::MinlpLSmear(ExtendedSystem& sys,  double prec, double ratio, minlplsmear_mode lsmode) : MinlpSmearSumRelative(sys,prec,ratio),
		lsmode(lsmode) {
		  mylinearsolver = new LPSolver(sys.nb_var);
}
*/
  inline MinlpLSmear::MinlpLSmear(ExtendedSystem& sys, double prec, LargestFirst& lf, bool gb, minlplsmear_mode lsmode) : MinlpSmearSumRelative(sys,prec,lf,gb),
		lsmode(lsmode) {
	mylinearsolver = new LPSolver(sys.nb_var);
}

  inline MinlpLSmear::MinlpLSmear(ExtendedSystem& sys, const Vector& prec, LargestFirst& lf, bool gb, minlplsmear_mode lsmode) : MinlpSmearSumRelative(sys,prec,lf,gb),
		lsmode(lsmode) {
	mylinearsolver = new LPSolver(sys.nb_var);
}

  /*
inline MinlpLSmear::MinlpLSmear(ExtendedSystem& sys, const Vector& prec, double ratio,minlplsmear_mode lsmode) : MinlpSmearSumRelative(sys,prec,ratio),
		lsmode(lsmode) {
	mylinearsolver = new LPSolver(sys.nb_var);
}
  */
} /* namespace ibex */

#endif /* __IBEX_LSMEAR__ */
