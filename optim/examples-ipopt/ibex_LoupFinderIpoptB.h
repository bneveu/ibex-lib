//============================================================================
//                                  I B E X
// File        : ibex_LoupFinderIpoptA.h
// Author      : Bertrand Neveu,Ignacio Araya
// Copyright   : IMT Atlantique (France)
// License     : See the LICENSE file
// Created     : Oct 26 2023
// Last Update : Nov 3 2023

//============================================================================

#ifndef __IBEX_LOUP_FINDER_IPOPTB_H__
#define __IBEX_LOUP_FINDER_IPOPTB_H__

#include <cassert>
#include <iostream>
#include "ibex_System.h"
#include "ibex_LoupFinder.h"
#include <sstream>
#include <fstream>
#include <stdio.h>
#include "/usr/local/include/coin-or/IpTNLP.hpp"
#include "/usr/local/include/coin-or/IpIpoptApplication.hpp"
#include "/usr/local/include/coin-or/IpSolveStatistics.hpp"
#include "/usr/local/include/coin-or/IpUtils.hpp"

#include "ibex_Optimizer.h"
#include "ibex_CellHeap.h"

using namespace Ipopt;
using namespace std;


namespace ibex {
 /**
 * \ingroup optim
 *
 * \brief Calls Ipopt to find a feasible point, which is a local minimum.
 *
 */

    class LoupFinderIpoptB: public TNLP, public LoupFinder  {

    private:
      SmartPtr<Ipopt::IpoptApplication> app;
      /** \brief the current box
       */
      IntervalVector the_box;
	/**
	 * \brief Symbolic gradient of the objective.
	 */
      Function* df=nullptr;;

	/**
	 * \brief Symbolic gradient of constraints.
	 */
      Function** dg=nullptr;;
    public:

      double optimalValue = POS_INFINITY;

        
      Vector solution;


      LoupFinderIpoptB(const System& sys, const System& normsys, const ExtendedSystem& extsys);

      virtual ~LoupFinderIpoptB();
      const System& sys;

      const System& normsys;
      const ExtendedSystem& extsys; // the extended system : used by recursive call of ipopt for correcting a solution 
      Optimizer* optimizer=nullptr; // optimizer data are used for building an optimizer for correcting the point returned by ipopt if it does not verify the constraints
      bool recursive_call=true; // boolean to prevent double recursion of optimizer ; when true, the optimizer can be recursevely called
      int ipopt_frequency=1; // frequency of ipopt calls in number of calls.
      int ipopt_calls=-1;
      bool force=false; // when true it forces Ipopt call (after a new loup is found by another loup finder
      int correction_nodes=0;  // additional nodes for correcting the point given by ipopt
      double correction_time=0.0; // additional time for correcting the point given by ipopt

      void correct_ipopt_sol (Vector&v, double& loup);
        /**
      	 * \brief Find a new loup in a given box.
      	 *
      	 * \see comments in LoupFinder.
      	 */
      	virtual std::pair<IntervalVector, double> find(const IntervalVector& box, const IntervalVector& loup_point, double loup);

        /************ methods required  by IpOpt ***********/

        /** Method to return some info about the nlp */
      virtual bool get_nlp_info(int& n, int& m, int& nnz_jac_g, int& nnz_h_lag, IndexStyleEnum& index_style);


        /** Method to return the bounds for my problem */
        virtual bool get_bounds_info(int n, Number* x_l, Number* x_u,
            int  m, Number* g_l, Number* g_u);

        /** Method to return the starting point for the algorithm */
        virtual bool get_starting_point(int n, bool init_x, Number* x,
            bool init_z, Number* z_L, Number* z_U,
            int m, bool init_lambda, Number* lambda);

        /** Method to return the objective value */
        virtual bool eval_f(int n, const Number* x, bool new_x, Number& obj_value);

        /** Method to return the gradient of the objective */
        virtual bool eval_grad_f(int n, const Number* x, bool new_x, Number* grad_f);

        /** Method to return the constraint residuals */
        virtual bool eval_g(int n, const Number* x, bool new_x, int m, Number* g);

        /** Method to return:
         *   1) The structure of the jacobian (if "values" is NULL)
         *   2) The values of the jacobian (if "values" is not NULL)
         */
        virtual bool eval_jac_g(int n, const Number* x, bool new_x,
            int m, int nele_jac, int* iRow, int* jCol, Number* values);

        /** Method to return:
         *   1) The structure of the hessian of the lagrangian (if "values" is NULL)
         *   2) The values of the hessian of the lagrangian (if "values" is not NULL)
         */
        virtual bool eval_h(int n, const Number* x, bool new_x, Number obj_factor,
			    int m, const Number* lambda, bool new_lambda,
			    int  nele_hess, int* iRow, int* jCol, Number* values);


        /** @name Solution Methods */
        //@{
        /** This method is called when the algorithm is complete so the TNLP can store/write the solution */
        virtual void finalize_solution(SolverReturn status,
            int n, const Number* x, const Number* z_L, const Number* z_U,
            int m, const Number* g, const Number* lambda, Number obj_value, const IpoptData* ip_data, IpoptCalculatedQuantities* ip_cq);

        /************************************************************************/
        bool integer_check(Vector& pt);
        bool is_inner(Vector& pt);
        double goal_ub(Vector& pt);
        void sysbound(Vector& pt);
        void sysbound(IntervalVector& vec);  
        void set_quadratic(bool quadratic);

    private:
        /**@name Methods to block default compiler methods.
         * The compiler automatically generates the following three methods.
         *  Since the default compiler implementation is generally not what
         *  you want (for all but the most simple classes), we usually
         *  put the declarations of these methods in the private section
         *  and never implement them. This prevents the compiler from
         *  implementing an incorrect "default" behavior without us
         *  knowing. (See Scott Meyers book, "Effective C++")
         *
         */
         //@{
        LoupFinderIpoptB();
        LoupFinderIpoptB(const LoupFinderIpoptB&);
        LoupFinderIpoptB& operator=(const LoupFinderIpoptB&);
        //@}


    };

}

#endif
