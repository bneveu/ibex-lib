//                                  I B E X
// File        : ibex_IpoptPreprocessing.cpp
// Author      : Bertrand Neveu
// Copyright   : IMT Atlantique (France)
// License     : See the LICENSE file
// Created     : May 14, 2012
// Last Update : Oct 20, 2023
//============================================================================

#include "ibex_Optimizer.h"
#include "ibex_Timer.h"
#include "ibex_LoupFinderIpopt.h"

#include <float.h>
#include <stdlib.h>
#include <iomanip>
using namespace std;

namespace ibex {


// analysis of results of a direct call of ipopt without passing through ampl
bool ipopt_direct_results(int n,string& status, double& obj, Vector & val){
  ifstream fic1("ipopt_res.txt");
   if (fic1.good()){

     string a;
     string b;
     string c;
     double v;
     for (string a ; fic1 >>a; ){
       if (a=="Objective...............:")
	 {cout << a;
	   fic1 >> obj;
	   cout << obj << endl;}
       else if (a=="EXIT:"){
	 cout << "exit " << a << endl;
	 fic1 >> a; fic1 >>b ; fic1 >> c;
	 cout << a << " " << b << " " << c << endl;
	 if (a=="Optimal" && b=="Solution" && c=="Found.")
	   status = "solved";
       }
       else if (a=="_svar[1]"){
	 fic1 >> v;
	 val[0]=v;
	 for (int i=1;i<n;i++){
	   fic1 >> a;
	   fic1 >> v;
	   val[i]=v;
	 }
       }
     }
     fic1.close();
     cout.precision(12);
     for (int i=0; i<n ;i++)
       cout << "var"<<i<< " " << val[i] << endl;

     cout << "status " << status << endl;
     if (status=="solved")return 1;
     else return 0;
   }
   return 0;
}
     

  
  double Optimizer:: check_ipopt(LoupFinder& loup_finder, Vector& v){
  double newub=POS_INFINITY;
    double loup=POS_INFINITY;
    if (loup_finder.integer_check(v)){
      //cout << "v after integer check " << v << endl;
      newub=loup_finder.goal_ub(v);
      //      cout << "newub " << newub << endl;
      if ( loup_finder.is_inner(v)){
	//	cout << "is_inner " << endl;
	if (integerobj)
	  {
	    Interval loupint (newub-loup_finder.integer_tolerance, newub+loup_finder.integer_tolerance);
	    loupint=integer(loupint);
	    //	    cout << "loupint " << loupint << endl;
	    if (!(loupint.is_empty()))
	      newub=loupint.ub();
	    else
	      newub=POS_INFINITY;
	  }
			      
	  loup= newub;
	  if (trace) cout << " loup= " << loup << endl;
      }
      
    }
    return loup;
  }
 

  
  void  Optimizer::ipopt_preprocessing(System& sys, const System& normsys, const ExtendedSystem & extsys){
   double obj=POS_INFINITY;
	Vector v(sys.nb_var);
	string status= "undefined";
	bool res;
	int res0;
	/*
	if (found!=std::string::npos){
	  cout << "appel direct nl " << endl;
	  string str = "time /libre/neveu/ibex-lib/ampl/ampl.linux-intel64/ipopt";
	  str.append(string(" "));
	  str.append(string(argv[1]));
	  str.append(" 'print_level 3 wantsol 2' > ipopt_res.txt");
	  //	  cout << str << endl;
	  res0= system(str.c_str());
          res=ipopt_direct_results(sys.nb_var,status,obj,v);
	}
	else{
	*/
	LoupFinderIpopt loupfindipopt(sys,normsys,extsys);
	loupfindipopt.optimizer=this;
	loupfindipopt.integerobj=integerobj;
	loupfindipopt.write_system_ampl(sys.box);
	res0=system("rm results.txt");

	res0= system(loupfindipopt.ipopt_ampl_run.c_str());
	res=loupfindipopt.ipopt_ampl_results(sys.nb_var,status,obj,v);
	preprocampltime= loupfindipopt.ampltime;
	preprocipopttime= loupfindipopt.ipopttime;
	  //	}

	if (trace) cout << " resultats ipopt " << status << endl;
	double initloup=POS_INFINITY;
	if (status=="solved"){
	    	  cout << "v " << v ;
		  cout << " obj " << obj << endl;
	  initloup= check_ipopt(loup_finder,v);
	  if (initloup==POS_INFINITY && sys.get_integer_variables()->size() < sys.nb_var)
	    loupfindipopt.correct_ipopt_sol(v,initloup);
	  if (trace) {
	    cout << "initloup " << initloup ;
	    if (initloup < POS_INFINITY) cout <<  " feasible point " << v;
	    cout << endl;
	  }
	  set_loup(initloup);
	  IntervalVector lp(v);
	  if (initloup == POS_INFINITY) lp.set_empty();
	  set_loup_point(lp);
	  
	  if (trace) cout << "preprocessing correction time " << loupfindipopt.correction_time << endl;
	  if (trace) cout << "preprocessing correction nodes " << loupfindipopt.correction_nodes << endl;

	}
  }

}

