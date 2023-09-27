//============================================================================
//                                  I B E X
// File        : ibex_LoupFinderIpopt.cpp
// Author      : Bertrand Neveu
// Copyright   : IMT Atlantique (France)
// License     : See the LICENSE file
// Created     : Aug 2023
// Last Update : Sep 27, 2023
//============================================================================

#include "ibex_LoupFinderIpopt.h"
#include <stdlib.h>
#include<fstream>

using namespace std;
namespace ibex {

  double expansion_precision=1.e-6;

  LoupFinderIpopt::LoupFinderIpopt(System& sys, const System& normsys, const ExtendedSystem& extsys) : sys(sys),normsys(normsys),extsys(extsys){

}

  void LoupFinderIpopt::correct_ipopt_sol (Vector&v, double& loup){
    if (recursive_call){
      recursive_call=false;
      IntervalVector box = sys.box;
      double eps=expansion_precision;
      IntervalVector boxsol(v.size());
      for ( int i=0; i< v.size() ; i++){
	double epsi = eps;
	if (fabs(v[i])>1) epsi= eps*fabs(v[i]);
	boxsol[i]= box[i] & Interval(v[i]- epsi, v[i]+ epsi);
      }
  //  cout << "boxsol " << boxsol << endl;
      CellHeap buffer(extsys);
      Optimizer opt(sys.nb_var,optimizer->ctc,optimizer->bsc,optimizer->loup_finder,buffer,extsys.goal_var(),optimizer->eps_x[0],optimizer->rel_eps_f, optimizer->abs_eps_f);
      opt.integerobj=optimizer->integerobj;
      opt.integer_tolerance=optimizer->integer_tolerance;
      opt.set_uplo(optimizer->get_uplo());
      opt.set_loup(optimizer->get_loup());
      opt.optimize(boxsol);
      recursive_call=true;
      correction_nodes+=opt.get_nb_cells();
      correction_time+=opt.get_time();
      if (opt.get_loup() < optimizer->get_loup()){
      //      cout << "new loup after correction " << opt.get_loup() << endl;
	loup= opt.get_loup();
	v = opt.get_loup_point().mid();
      }
    
    }



}


  void LoupFinderIpopt::write_system_ampl(const IntervalVector& box){
	ofstream ficsys ("system.mod", ofstream::trunc);
	
	ficsys.precision(15);
	IntervalVector initbox=sys.box;
	sys.box=box;
	ficsys << sys  << endl;
	sys.box=initbox;
	ficsys << "solve > amplout.txt;" << endl ;
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
  if (recursive_call){
	int n=sys.nb_var;
	Vector loup_point(n);
	double loup = current_loup;

	Vector pt(n);
	bool loup_changed=false;
	string status= "undefined";
        double obj=POS_INFINITY;

	write_system_ampl(box);
	int res0=system("rm results.txt");
	res0= system(ipopt_ampl_run.c_str());
	int res=ipopt_ampl_results(n,status,obj,pt);

	//	cout << " resultats ipopt " << status;
	double ipoptloup=POS_INFINITY;
	if (status=="solved"){
	  //	  bool _is_inner = normsys.is_inner(box);
	  //	  cout << " obj " << obj;
	  bool _is_inner=false;
	  if( check(normsys,pt,loup,_is_inner))
	    { loup_changed=true;
	      loup_point=pt;
	      if (optimizer->trace && !optimizer->integerobj)   cout << "*** ipopt      " ;
	    }
	  else{
	    
	    correct_ipopt_sol(pt, ipoptloup);
	    if (ipoptloup < current_loup)
	      {loup_changed=true;
		loup_point=pt;
		loup=ipoptloup;
		 if (optimizer->trace && !optimizer->integerobj) cout << "*** ipopt+corr " ;
	      }
	  }
	    
	  
	}
	if (loup_changed)
	  return std::make_pair(loup_point,loup);
	else
	  throw NotFound();
  }
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
