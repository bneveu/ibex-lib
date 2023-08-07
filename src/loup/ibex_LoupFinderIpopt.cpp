//============================================================================
//                                  I B E X
// File        : ibex_LoupFinderIpopt.cpp
// Author      : Bertrand Neveu
// Copyright   : IMT Atlantique (France)
// License     : See the LICENSE file
// Created     : May 14, 2012
// Last Update : Jul 09, 2017
//============================================================================

#include "ibex_LoupFinderIpopt.h"
#include <stdlib.h>
#include<fstream>

using namespace std;
namespace ibex {

  LoupFinderIpopt::LoupFinderIpopt(System& sys, const System& normsys) : sys(sys),normsys(normsys){

}


  void LoupFinderIpopt::write_system_ampl(const IntervalVector& box){
	ofstream ficsys ("system.mod", ofstream::trunc);
	
	//	cout << *sys  <<endl;
	ficsys.precision(15);
	IntervalVector initbox=sys.box;
	sys.box=box;
	ficsys << sys  << endl;
	sys.box=initbox;
	ficsys << "solve > toto.txt;" << endl ;
	ficsys << "option display_precision 20;" << endl;
	ficsys << "display solve_result > results.txt;" << endl;
	ficsys << "display obj > results.txt; "<< endl;
	ficsys << "display _total_solve_time > results.txt;" << endl;
	ficsys << "display _ampl_time > results.txt; " << endl;

	for (int i=0; i<sys.args.size(); i++) {
		const ExprSymbol& x = sys.args[i];
		ficsys << " display " << x << "> results.txt;" << endl;
	};
	ficsys.close();
}

  bool LoupFinderIpopt::ipopt_ampl_results(int n,string& status, double& obj, Vector & v){
    ifstream fic1 ("results.txt");
    if (fic1.good()){
      string a;
      int j;
      double otime=0;
      fic1 >> a;
      fic1 >> a;
      fic1 >> status; 
      //      cout << "status " << status << endl;
      fic1 >> a;
      fic1 >> a; 
      fic1 >> obj;
      //       cout << "obj " << obj << endl;
      
      fic1 >> a;
      fic1 >> a;
      fic1 >> otime;
      //cout << "ipopttime " << otime;
      ipopttime +=otime;

      fic1 >> a;
      fic1 >> a;
      fic1 >> otime;
      //      cout << "ampltime " << otime << endl;
      ampltime +=otime;

      //      cout << "ampltime " << ampltime << " ipopttime " << ipopttime << endl;
      
      for (int i =0 ; i< n; i++){
	fic1 >> a;
	fic1 >> a;
	fic1 >> v[i];
	//	cout << "j" << a << "  v["<<i<<"] " <<  v[i] << endl;
      }
    
      fic1 >> a;
      fic1.close();
      return true;
    }
    return false;
  }

 
std::pair<IntervalVector, double> LoupFinderIpopt::find(const IntervalVector& box, const IntervalVector& current_loup_point, double current_loup) {

	int n=sys.nb_var;
	Vector loup_point(n);
	double loup = current_loup;

	Vector pt(n);
	bool loup_changed=false;
	string status= "undefined";
        double obj=POS_INFINITY;

	write_system_ampl(box);
	
	int res0= system("/libre/neveu/ampl/ampl model_ipopt.run > toto.txt");
	int res=ipopt_ampl_results(sys.nb_var,status,obj,pt);

	//	cout << " resultats ipopt " << status;
	double ipoptloup=POS_INFINITY;
	if (status=="solved"){
	  //	  bool _is_inner = normsys.is_inner(box);
	  //	  cout << " obj " << obj;
	  bool _is_inner=false;
	  if( check(normsys,pt,loup,_is_inner))

	      
	  /*  A FAIRE ...
	  if (initloup==POS_INFINITY)
	    initloup=correct_ipopt_sol(v,o,*sys);
	  */
	    { loup_changed=true;
	      loup_point=pt;
	      cout << "*** ipopt " ;
	    }
	}
	if (loup_changed)
	  return std::make_pair(loup_point,loup);
	else
	  throw NotFound();
}


  
  bool LoupFinderIpopt::integer_check(Vector& pt)
           {return integer_and_bound_check(normsys,pt);}
  bool LoupFinderIpopt::is_inner(Vector& pt)
           {return is_inner0(normsys,pt);}
  double LoupFinderIpopt::goal_ub(Vector& pt)
           {return goal_ub0(normsys,pt);}
  void LoupFinderIpopt::sysbound(Vector& pt) {bound_check(normsys,pt);}
    void LoupFinderIpopt::sysbound(IntervalVector& vec) {bound_check(normsys,vec);}
} /* namespace ibex */
