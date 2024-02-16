#include <iostream>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <limits>
#include "ibex_LoupFinderIpoptB.h"
#include "ibex_DefaultOptimizerConfig.h"
#include "ibex.h"

#include "ibex_NormalizedSystem.h"
#include <cassert>

using namespace Ipopt;
using namespace std;


namespace ibex {
  double expansion_precisionB=1.e-6;
  double ipopt_diam=1.e8;
  
  LoupFinderIpoptB::LoupFinderIpoptB(const System& sys,const System& normsys, const ExtendedSystem& extsys) : sys(sys), normsys(normsys), extsys(extsys), solution(sys.nb_var), the_box(sys.box), ipopt_box(sys.box) {
	try {
		df = new Function(*sys.goal,Function::DIFF);
		//                cout << " nb ctr " << sys.nb_ctr << endl;
		if (sys.nb_ctr>0) {
			dg = new Function*[sys.f_ctrs.image_dim()];

			for (int i=0; i<sys.f_ctrs.image_dim(); i++) {
				dg[i] = new Function(sys.f_ctrs[i], Function::DIFF);
			}
		} else {
			dg = NULL;
		}

	}
	catch(Exception&) {
		//TODO: replace with ExprDiffException.
		// Currently, DimException is also sometimes raised.
		cerr << "Warning: symbolic differentiation has failed ==> ipopt hessian disabled" << endl;
		df = NULL;
		dg = NULL;
	}
	//	cout << " fin initialisation gradients " << endl;
      //The IpOPT initialization should be here
      app = IpoptApplicationFactory();
      app->RethrowNonIpoptException(true);
     

      // Change some options
      // Note: The following choices are only examples, they might not be
      // suitable for your optimization problem.
      
      //      app->Options()->SetNumericValue("tol", 1e-8);
      //      app->Options()->SetNumericValue("constr_viol_tol", 1e-10);
      // app->Options()->SetStringValue("mu_strategy", "adaptive");
      app->Options()->SetStringValue("output_file", "ipopt.out");
      app->Options()->SetIntegerValue("print_level",0);
      app->Options()->SetIntegerValue("max_iter",1000);
      app->Options()->SetStringValue("linear_solver","mumps");
      // Initialize the IpoptApplication and process the options
      
      ApplicationReturnStatus status;
      status = app->Initialize();
      //      cout << "initialisation ipopt " << endl;
      if (status != Solve_Succeeded) {
      	std::cout << std::endl << std::endl << "*** Error during initialization!" << std::endl;
      	exit(0);
      }
      }

    LoupFinderIpoptB::~LoupFinderIpoptB()
    {}


    std::pair<IntervalVector, double> LoupFinderIpoptB::find(const IntervalVector& box, const IntervalVector& loup_point, double loup) {
      double loup0=loup;
      if (ipopt_calls== -1) ipopt_box=box;  // initialization of ipopt_box before first call.
      if (sys.minlp) ipopt_box=box;  // in case of Minlp, ipopt is run with the current box.
 
      the_box  = box;
      double newloup=false;
      IntervalVector loup_point0=loup_point;
      if (recursive_call){
	ipopt_calls++;
	if (
	    (!sys.minlp ||box.max_diam()<= ipopt_diam)
	    &&
	    (ipopt_calls%ipopt_frequency==0 || ipopt_calls==10 || ipopt_calls==20  || ipopt_calls==50 || (force && box.max_diam()<= ipopt_diam) )
	    ){
	
	  //	  cout << "nb_cells " <<  optimizer->get_nb_cells() << endl;
	  ApplicationReturnStatus status = app->OptimizeTNLP(this);
	  force=0;
	  //	  cout << " status " << status << endl;
	  //	  if (status == Solve_Succeeded) {
	  //        	std::cout << std::endl << std::endl << "*** The problem solved!" << std::endl;
	  //	  cout << " optimal_value ipopt "<< optimalValue << endl;
	  //	  cout << " solution ipopt " << solution << endl;
	  if( optimalValue< loup){
	    
	    for (int i=0; i< box.size(); i++) {
	      bound_check_i(sys,solution,i);
	    }
	    //	    cout << " pt avant check "  << solution << endl;
	    if (check(normsys,solution,loup,false))
		  { 
		    loup_point0=solution;
		    if (optimizer->trace) 
		    cout << "*** ipopt      " ;
		  }

	    if (sys.get_integer_variables()->size() < sys.nb_var){
	    //	    else if (optimalValue< loup && sys.get_integer_variables()->size() < sys.nb_var){
	      //	      cout << " solution " << solution << endl;
	      //	      cout << " optimal value " << optimalValue << " loup " << loup << endl;
		  double ipoptloup=POS_INFINITY;
		  correct_ipopt_sol(solution, ipoptloup);
		  if (ipoptloup < loup)
		    {
		      loup_point0=solution;
		      loup=ipoptloup;
		      if (optimizer->trace)
			cout << "*** ipopt+corr " ;
		      newloup=true;
		    }
	    }
	  }
	}
	
      }
        
      if (loup < loup0)
	return std::make_pair(loup_point0, loup);
      else  throw NotFound();


}



    bool LoupFinderIpoptB::get_nlp_info(int& n, int& m, int& nnz_jac_g,
					int& nnz_h_lag, IndexStyleEnum& index_style)
    {
        n = sys.nb_var;

        m = sys.nb_ctr;

        nnz_jac_g = n * m;

        nnz_h_lag = (n * (n+1))/2;

        index_style = TNLP::C_STYLE;

        return true;
    }



    // returns the variable bounds
    bool LoupFinderIpoptB::get_bounds_info(int n, Number* x_l, Number* x_u, int m, Number* g_l, Number* g_u)
    {
        // here, the n and m we gave IPOPT in get_nlp_info are passed back to us.
        // If desired, we could assert to make sure they are what we think they are.


        assert(n == sys.nb_var);
        assert(m == sys.nb_ctr);
	//        cout << "ipopt box " << ipopt_box << endl;
        Array<NumConstraint> contrains = sys.ctrs;
        for (int i = 0; i < n; i++) {
	  x_l[i] = ipopt_box[i].lb();   // Ipopt has for domain the initial box after contraction

        }


        for (int i = 0; i < n; i++) {
	  x_u[i]=ipopt_box[i].ub();
        }
	
	
        for (int i = 0; i < m; i++)
        {

            if (contrains[i].op == EQ)
            {
	      g_u[i] = 0;
	      g_l[i] = 0;
	    }

            else if (contrains[i].op == GEQ)
            {
                g_u[i] = 2.e19;
                g_l[i] = 0;

            }

            else if (contrains[i].op == LEQ)
            {
                g_u[i] = 0;
                g_l[i] = -2.e19;

            }
        }

        return true;
    }


    // returns the initial point for the problem
    bool LoupFinderIpoptB::get_starting_point(int n, bool init_x, Number* x,
        bool init_z, Number* z_L, Number* z_U,
        int
        m, bool init_lambda,
        Number* lambda)
    {
        // the starting point is the center of the current box, and in case of a
        // new loup found by another loup finder (force=true) , the corresponding point (solution)

        assert(init_x == true);

        assert(init_z == false);
        assert(init_lambda == false);

        // initialize to the given starting point
	//	ibex::Vector v = the_box.random();
	ibex::Vector v = the_box.mid();
	//	cout << " force " << force << endl;
	if (force)  v=solution;

	//	cout << " v " << v << endl;
	for (int i=0; i<n ; i++)
	  x[i]=v[i];
        return true;
    }

    // returns the value of the objective function
    bool LoupFinderIpoptB::eval_f(int n, const Number* x, bool new_x, Number& obj_value)
    {
        assert(n == the_box.size());
        ibex::Vector v(n);

        for (int i = 0; i < n; i++) {
            v[i] = x[i];
        }

        obj_value = sys.goal_ub(v);
     
        return true;
    }


    // return the gradient of the objective function grad_{x} f(x)

    bool LoupFinderIpoptB::eval_grad_f(int n, const Number* x, bool new_x, Number* grad_f)
    {
      assert(n == the_box.size());

      
        ibex::IntervalVector v(n);

        Function* funcion = sys.goal;
	//	cout << " goal " << *funcion << endl;
	
        for (int i = 0; i < v.size(); i++)
	  {
            v[i] = x[i];
	  }
	ibex::IntervalVector gradiente(n);
  	gradiente= funcion->gradient(v);
        //  gradiente= df->eval_vector(v);
        //	  cout << " gradient " << gradiente << endl;
	//        cout << " grad size " << gradiente.size() << endl;
	//	cout << " v " << v << endl;
        //        cout << "df[0] " << df[0] << endl;
        //  cout << "df[0] val" << *df[0].eval(v)[0].mid() <<  " " << gradiente[0].mid() << endl;
        for (int i = 0; i < n; i++)
        {
	  grad_f[i] = gradiente[i].mid();
	  //	  double val=df[i].eval(v).mid();
	  //	  cout << i << " df[0] val" << val << " "  << gradiente[i].mid()        << endl;
	  //	  grad_f[i] = df[i].eval(v).mid();

        }


        return true;
    }

    // return the value of the constraints: g(x)
    bool LoupFinderIpoptB::eval_g(int n, const Number* x, bool new_x, int m, Number* g)
    {
        assert(n == sys.nb_var);
        assert(m == sys.nb_ctr);

        ibex::IntervalVector v(n);

        for (int i = 0; i < v.size(); i++)
            v[i] = x[i];


        for (int i = 0; i < m; i++)
        {
	    g[i] = sys.ctrs[i].f.eval(v).mid();
        }
        return true;
    }

    // return the structure or values of the jacobian
    bool LoupFinderIpoptB::eval_jac_g(int n, const Number* x, bool new_x, int m, int nele_jac, int* iRow, int* jCol, Number* values)
    {

        ibex::IntervalVector v(n);

        int cont = 0;
        if (values == NULL) {
            // return the structure of the jacobian
            int p = 0;
            for (int i = 0; i < m; i++)
                for (int j = 0; j < n; j++) {
                    iRow[p] = i;
                    jCol[p] = j;
                    p++;
                }


        }
        else {
            // return the values of the jacobian of the constraints
            for (int i = 0; i < v.size(); i++) {
                v[i] = x[i];
               // cout << v[i] << endl;
            }

            Array<NumConstraint> contrains = sys.ctrs;
       
            for (int c = 0; c < m; c++)
            {

	      IntervalVector grad(n);
	      grad= sys.ctrs[c].f.gradient(v);
	       //               cout << "restriccion:" << c << endl;
	       //               cout <<grad<< endl;

	       for (int k = 0; k < n; k++)
                    {

                            values[cont] = grad[k].mid();
			    //		    cout << cont << " " << values[cont] << endl;
			    cont++;

                    }

            }

        }

        return true;
    }
 

  
    bool LoupFinderIpoptB::eval_h(int n, const Number* x, bool new_x, Number obj_factor, int m, const Number* lambda, bool new_lambda, int nele_hess, int* iRow, int* jCol, Number* values)
    {

      if( values == NULL )
	{
      // return the structure. This is a symmetric matrix, fill the lower left
      // triangle only.
 
      // the hessian  is actually dense
	  int idx = 0;
	  for( int row = 0; row < n; row++ )
	    {
	      for( int col = 0; col <= row; col++ )
		{
		  iRow[idx] = row;
		  jCol[idx] = col;
		  idx++;
		}
	    }
	  assert(idx == nele_hess);
        }
  
      else {
	ibex::IntervalVector v(n);
	for (int i = 0; i < n; i++) {
	  v[i] = x[i];

	}

	IntervalMatrix jacobianobj(n,n);
	df->jacobian(v,jacobianobj);

	int nn=(n*(n+1))/2;

	for (int i=0; i< nn ; i++)
	  values[i]=0;
	int l=0;
	for (int i = 0; i < n; i++){
          for (int j=0; j<=i; j++){
	    //	    cout << " jacobianobj " << i << "  " << j << " " << jacobianobj[i][j].mid() << endl;
	   
	    values[l]+=obj_factor*jacobianobj[i][j].mid();
	    l++;
	  }
	}
	for (int k=0 ; k< m; k++){
          l=0;
	  Function* f= dg[k];
	  //	  cout << "function " << *f << endl;
	  IntervalMatrix jacobianct(n,n);
	  f->jacobian(v, jacobianct);
	  for (int i = 0; i < n; i++){
	    for (int j=0; j<=i; j++){
	      // cout << k << " jacobianct " << i << "  " << j << " " << jacobianct[i][j].mid() << endl;
	      values[l]+=lambda[k]*jacobianct[i][j].mid();
	      l++;
	    }
	  }
	}
      }
      
      return true;
    }
    

    void LoupFinderIpoptB::finalize_solution(SolverReturn status,
        int n, const Number* x, const Number* z_L, const Number* z_U,
        int m, const Number* g, const Number* lambda,
        Number obj_value,
        const IpoptData* ip_data,
        IpoptCalculatedQuantities* ip_cq)
    {
        // here is where we would store the solution to variables, or write to a file, etc
        // so we could use the solution.

        // For this example, we write the solution to the console
        /*
        std::cout << std::endl << std::endl << "Solution of the primal variables, x" << std::endl;
        for (int
            i = 0; i < n; i++) {
            std::cout << "x[" << i << "] = " << x[i] << std::endl;
        }

        std::cout << std::endl << std::endl << "Solution of the bound multipliers, z_L and z_U" << std::endl;
        for (int
            i = 0; i < n; i++) {
            std::cout << "z_L[" << i << "] = " << z_L[i] << std::endl;
        }
        for (int
            i = 0; i < n; i++) {
            std::cout << "z_U[" << i << "] = " << z_U[i] << std::endl;
        }

        std::cout << std::endl << std::endl << "Objective value" << std::endl;
        std::cout << "f(x*) = " << obj_value << std::endl;

        std::cout << std::endl << "Final value of the constraints:" << std::endl;
        for (int i = 0; i < m; i++) {
            std::cout << "g(" << i << ") = " << g[i] << std::endl;
        }
        */

        optimalValue = obj_value;

        for (int i = 0; i < n; i++) {
            solution[i] = x[i];
        }

    }


  bool LoupFinderIpoptB::integer_check(Vector& pt)
           {return integer_and_bound_check(normsys,pt);}
  bool LoupFinderIpoptB::is_inner(Vector& pt)
           {return is_inner0(normsys,pt);}
  double LoupFinderIpoptB::goal_ub(Vector& pt)
           {return goal_ub0(normsys,pt);}
  void LoupFinderIpoptB::sysbound(Vector& pt) {bound_check(normsys,pt);}
    void LoupFinderIpoptB::sysbound(IntervalVector& vec) {bound_check(normsys,vec);}
  // OK only for QP problems (quadratic objective and linear constraints)
  void LoupFinderIpoptB::set_quadratic(bool quadra){
    if (quadra) app->Options()->SetStringValue("hessian_constant", "yes");
  }
  
  void LoupFinderIpoptB::correct_ipopt_sol (Vector&v, double& loup){
    if (recursive_call){
      
      recursive_call=false;
      IntervalVector box = sys.box;
      double eps=expansion_precisionB;
      //      if (sys.minlp) eps=0.01;
      IntervalVector boxsol(v.size());
      for ( int i=0; i< v.size() ; i++){
	double epsi = eps;
	if (fabs(v[i])>1) epsi= eps*fabs(v[i]);
	boxsol[i]= sys.box[i] & Interval(v[i]- epsi, v[i]+ epsi);
      }
      CellHeap buffer(extsys);

      Optimizer opt(sys.nb_var,optimizer->ctc,optimizer->bsc,optimizer->loup_finder,buffer,extsys.goal_var(),optimizer->eps_x[0],optimizer->rel_eps_f, optimizer->abs_eps_f);

      opt.integerobj=optimizer->integerobj;

      opt.set_uplo(optimizer->get_uplo());
      opt.set_loup(optimizer->get_loup());
      opt.timeout=1;
      //      cout << " boxsol " << boxsol << endl;
      opt.optimize(boxsol);
      recursive_call=true;
      correction_nodes+=opt.get_nb_cells();
      correction_time+=opt.get_time();
      if (optimizer->trace)
	cout << " correction nodes " << correction_nodes << " correction_time " << correction_time << endl;
      if (opt.get_loup() < optimizer->get_loup()){
	if (optimizer->trace)
	  cout << "new loup after correction " << opt.get_loup() << endl;
	loup= opt.get_loup();
	v = opt.get_loup_point().mid();
      }
    
    }
}


}
