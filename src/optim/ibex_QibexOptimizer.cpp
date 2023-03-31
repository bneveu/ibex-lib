//                                  I B E X
// File        : ibex_QibexOptimizer.cpp
// Author      : Bertrand Neveu
// Copyright   : IMT Atlantique (France)
// License     : See the LICENSE file
// Created     : Sep 01, 2021
// Last Update : Oct 24, 2022
//============================================================================

#include "ibex_QibexOptimizer.h"

#include <float.h>
#include <stdlib.h>
#include <iomanip>


using namespace std;

namespace ibex {

  QibexOptimizer::QibexOptimizer(int n, Ctc& ctc, Bsc& bsc, LoupFinder& finder,
		CellBufferOptim& buffer,
				 int goal_var, double minwidth, double tolerance, double eps_x, double rel_eps_f, double abs_eps_f) : Optimizer (n,ctc,bsc,finder,buffer,goal_var,eps_x,rel_eps_f,abs_eps_f), minwidth(minwidth), tolerance(tolerance) {init();}

  void QibexOptimizer::init(){
    ampltime=0;
    solvertime=0;
    int res0=system("rm results.txt");
    int res=system("/libre/neveu/ampl/ampl model_quad_relax_init.run");
    ifstream fic4 ("associ.txt");
    string a;
    int j;
    fic4 >> a;
    fic4 >> a;
    fic4 >> n_y;
    for (int i=0 ; i< n_y; i++)
      associ.push_back(-1);
    fic4 >> a;
    fic4 >> a;
    fic4 >> a;
    for (int i =0 ; i< n_y; i++){
      fic4 >> j;
      fic4 >> associ[j-1];
    }
    fic4 >> a;
    fic4.close();

    ifstream fic5 ("assocj.txt");
    fic5 >> a;

    fic5 >> a;
    fic5 >> n_y;
    for (int i=0 ; i< n_y; i++)
      assocj.push_back(-1);
    
    fic5 >> a;
    fic5 >> a;
    fic5 >> a;
    for (int i =0 ; i< n_y; i++){
      fic5 >> j;
      fic5 >> assocj[j-1];
    }
    fic5 >> a;
    fic5.close();
    ifstream fic6("diag_ref.txt");
    fic6 >> a;
    fic6 >> a;
    fic6 >> n_x;
    fic6 >> a;
    fic6 >> a;
    fic6 >> a;
    for (int i=0 ; i< n_x; i++){
      ref_diag_coefs.push_back(-1);}
    for (int i=0 ; i< n_x; i++){
      fic6 >> j;
      fic6 >>ref_diag_coefs[j-1];
    }
    
    ifstream fic8("coef_ref.txt");
    for (int i=0 ; i< n_x; i++){
      vector<int> v ;
      for (int j=0 ; j< n_x; j++){
	v.push_back(-1);
      }
      ref_coefs.push_back(v);
    }
    
    for (int i=0 ; i< n_x; i++){
      for (int j=i ; j< n_x; j++)
	fic8 >> ref_coefs[i][j];
    }
  

    /*
    res=system("/libre/neveu/ampl/ampl model_quad_init.run");
    ifstream fic7("diag_hessian.txt");

    fic7 >> a;
    fic7 >> a;
    fic7 >> n_x;
    fic7 >> a;
    fic7 >> a;
    fic7 >> a;
    for (int i=0 ; i< n_x; i++){
      hessian_diag_coefs.push_back(-1);}
    for (int i=0 ; i< n_x; i++){
      fic7 >>j;
      fic7 >>hessian_diag_coefs[j-1];
    }
    */

 for (int i=0 ; i< n_x; i++){
      hessian_diag_coefs.push_back(-1);
 }
 ifstream fic7("diag_hessian.txt");
 for (int i=0; i < n_x; i++){
   fic7 >>  hessian_diag_coefs[i];
 }
 fic7.close();
 /*
 cout << " hessian diag coefs " ;
 for (int i=0; i < n_x; i++){
   cout <<  hessian_diag_coefs[i] << " " ;
 }
 cout << endl;
 */
 /*
    IntervalVector gradient1(n_x);
    IntervalVector boxn(n_x);
    for (int i=0; i < n_x; i++)
	   boxn[i]=Interval(0.0,0.0);
    gradient=  goal->gradient(boxn);
    for (int i=0 ; i< n_x; i++){
      hessian_diag_coefs.push_back(-1);
    
  }
 */
  }
  bool QibexOptimizer::update_loup(const IntervalVector& box, BoxProperties& prop) {
    if (loupfinderp)
      return Optimizer::update_loup(box,prop);
    else
      return false;
  }

  void QibexOptimizer::qibex_relaxation_call(const IntervalVector & box, double objlb){
    ofstream fic ("bound.ampl", ofstream::trunc);
    IntervalVector pt(n_x);
    for (int i=0; i< n_x; i++){
      pt[i]=Interval(box[i].lb()-1.e-9,box[i].ub()+1.e-9);
	  
    }
    loup_finder.sysbound(pt);  // to ensure that pt is inside the initial box.
	
    fic << "data;" << endl;
    fic << "param ub :=" << endl;

    for (int i=0; i< n_x; i++){
      fic << i+1 << " " << std::setprecision(9) << pt[i].ub() << endl;
    }
    fic << ";" <<endl;
    fic << "param lb :=" << endl; 
    for (int i=0; i< n_x; i++){
      
      fic << i+1 << " " << std::setprecision(9) << pt[i].lb() << endl;
    }
    fic << ";" << endl;
    fic << "param objlb :=" << endl;
    fic << objlb << ";" << endl;
    fic << "end ;";
    
    fic.close();

    int res;
    if (rigor)
      res=system("/libre/neveu/ampl/ampl model_quad_relax_rigueur.run > amplout");
    else
      res=system("/libre/neveu/ampl/ampl model_quad_relax.run > amplout");
      
  }

  bool QibexOptimizer::qibex_relaxation_results(string& status, double& newlb, Vector & v, Vector& w){
    ifstream fic1 ("results.txt");
    if (fic1.good()){
      string a;
      int j;
      double newobj;
      int lbused =0;
      double otime=0;
      fic1 >> a;
      fic1 >> a;
      fic1 >> status; 
      
      fic1 >> a;
      fic1 >> a;
      fic1 >> newlb;
      
      fic1 >> a;
      fic1 >> a; 
      fic1 >> newobj;

      fic1 >> a;
      fic1 >> a;
      fic1 >> lbused;

      fic1 >> a;
      fic1 >> a;
      fic1 >> otime;
      solvertime +=otime;

      fic1 >> a;
      fic1 >> a;
      fic1 >> otime;
      ampltime +=otime;

      
      
      
      //      cout << " status " << status << endl;
      //      cout << " newlb " << newlb << endl;
      //      cout << " newobj " << newobj << endl;
      //      cout << "lbused " << lbused << " rigor " << rigor << endl;
      //      cout.precision(12);
      if (trace && lbused==0)
	cout << " relaxation not useful " << endl;
      if (rigor && lbused==0)
	newlb=-1.e300;
      if (trace && newlb < newobj && status=="solved" && rigor && lbused ==1)
	cout << "CORRECTION " << newlb << " " << newobj <<  "  " << newobj-newlb << endl;
      
      fic1 >> a;
      fic1>> a;
      fic1>> a;
      for (int i =0 ; i< n_x; i++){
	fic1 >> j;
	fic1 >> v[j-1];
	//	    cout << "j" << j << "  v[j-1] " <<  v[j-1] << endl;
      }
      fic1 >> a;
      fic1 >> a;
      fic1 >> a;
      fic1 >> a;

      for (int i =0 ; i< n_y; i++){
	fic1 >> j;
	fic1 >> w[j-1];
	//cout << "j" << j << " w[j-1] " <<  w[j-1] << endl;
	
      }
      fic1 >> a;
      fic1.close();
      return true;
    }
    return false;
  }
	  

  
  tuple<Vector,Vector,double> QibexOptimizer::qibex_relaxation_analysis(IntervalVector & box, string & status){

	double newlb = NEG_INFINITY;

	int n_y_max=n_x+(n_x*(n_x+1))/2;

	Vector v (n_x);
	
	Vector w(n_y_max);
	
	if (qibex_relaxation_results(status, newlb, v, w)){
	  if (trace){
	   cout << "status " << status << endl;
	   cout << " nb cells " << nb_cells << endl;
	  }
	  //	  cout << " v " << v << endl;
	  //	  cout << "n_x " << n_x << " n_y " << n_y << " v " << v << endl;
	  //	 cout << " w " << w << endl;
	
	  /*
	    if (status=="infeasible")  //too strong for cplex (useful for ipopt ??)
	    box.set_empty();
	  */
	 
	   //	   if  (status!= "solved")
	   if  (status!= "solved" && status!= "'solved?'")
	    {
             //cout << b << "  " << newlb << endl;
             newlb = NEG_INFINITY;
           }

	}

	if (integerobj && newlb != NEG_INFINITY){
	  
	  if (rigor) {
	    if (newlb != std::floor(newlb))
		newlb= std::ceil(newlb);
	  }
	  
	  else{
	    double epsinteger=1.e-4;
	    if (newlb >= std::floor(newlb) + epsinteger)
	      newlb= std::ceil(newlb);
	    else
	      newlb=std::floor(newlb) ;
	  }
	}

        tuple <Vector,Vector,double> triple(v,w,newlb);
	return triple;
  }
    
  
  double  QibexOptimizer:: compute_ratio (const IntervalVector& box,const Vector& v, int i) {
    double bisectionpoint = 0.75* box[i].mid() + 0.25 *v[i];
    //    cout << " ratio " << (bisectionpoint-box[i].lb())/(box[i].ub()-box[i].lb()) << endl;
    return ((bisectionpoint-box[i].lb())/(box[i].ub()-box[i].lb()));
  }


  // code utilisé uniquement si bisector = qibex..."
  int QibexOptimizer:: compute_var_to_bisect(const IntervalVector& box,  const Vector& v, const Vector & w, double& gap0) {
    double mingap=tolerance;
    double epsbound=1.e-4;
    double gap=0.0;
    double width=0.0;
    int var=-1;

    for (int i =0; i < n_y; i++){
      if (associ[i]== assocj[i]){
	//	cout << i << "  " << associ[i] << "  " << assocj[i] << endl;
	//gap = fabs((w[i+n_x]-v[associ[i]-1]*v[assocj[i]-1]) * (hessian_diag_coefs[associ[i]-1] - ref_diag_coefs[associ[i]-1]));
	//gap = fabs((w[i+n_x]-v[associ[i]-1]*v[assocj[i]-1])*ref_diag_coefs[associ[i]-1]);
	gap = fabs(w[i+n_x]-v[associ[i]-1]*v[assocj[i]-1]);
	if (gap > gap0) gap0=gap;
	width =  box[associ[i]-1].ub() - box[associ[i]-1].lb() ;
	//	cout << "gap " << gap << " i " << i << " w[i+n_x] " << w[i+n_x] << " " << v[associ[i]-1] << " " << v[assocj[i]-1] << endl;
	if (gap > mingap  
	    && v[associ[i]-1] > box[associ[i]-1].lb() + epsbound
	    && v[associ[i]-1] < box[associ[i]-1].ub() - epsbound
	    && width > minwidth
	    )
	    
	    {
	      //cout << " gap " << gap << " maxgap " << maxgap << " var " <<  associ[i]-1 << endl;
	      mingap=gap;
	      var=associ[i]-1;
	    }
      }
    }
    //    cout << " variable " << var << endl;
    
    if (var==-1){
      //          cout << " cas produit " << endl;
      mingap=tolerance;
      for (int i =0; i < n_y; i++){
	if (associ[i]!= assocj[i]){
	  gap = fabs(w[i+n_x]-v[associ[i]-1]*v[assocj[i]-1]);
	  //	  gap = fabs((w[i+n_x]-v[associ[i]-1]*v[assocj[i]-1])* ref_coefs[associ[i]-1][assocj[i]-1]);
	if (gap > gap0) gap0=gap;
	//	cout << "gap " << gap << " i " << i << " w[i] " << w[i] << " " << v[associ[i]-1] << " " << v[assocj[i]-1] << endl;
	if (gap > mingap
	    && v[associ[i]-1] > box[associ[i]-1].lb() + epsbound
	    && v[associ[i]-1] < box[associ[i]-1].ub() - epsbound
	    && box[associ[i]-1].ub() - box[associ[i]-1].lb() > minwidth
	    )
	    {
	      // cout << " gap " << gap << " maxgap " << maxgap << " var " <<  associ[i]-1 << endl;
	      mingap=gap;
	      var=associ[i]-1;
	    }
	}
      }
    }
    //    cout << " computed var " << var << endl;
    return var;
  
  }

  double QibexOptimizer::qibex_loupfinder(Vector& v){
    double newub=POS_INFINITY;
    double ymax=POS_INFINITY;
    if (loup_finder.integer_check(v)){
      newub=loup_finder.goal_ub(v); 
      if (newub < loup && loup_finder.is_inner(v)){
	loup= newub;
	if (trace) {
	  cout << "                    ";
	  cout << "*** \033[32m loup= " << loup << "\033[0m" << endl;
	  //	  cout << "*** \033[32m loup_point= " << v << "\033[0m" << endl;
	}
	loup_point = v;
	loup_changed=true;

	ymax=compute_ymax();
      }
    }
    return ymax;
  }


  void QibexOptimizer::qibex_bisection_choice (Cell& c, IntervalVector & qcp_box, Vector& v, Vector& w){
   int var_to_bisect=-1;
 
   double ratio=0.45;
   // double ratio=0.5;
   double gap0=0.0;

   var_to_bisect=compute_var_to_bisect(qcp_box,  v,  w, gap0);
   if (var_to_bisect != -1)
     ratio=compute_ratio(qcp_box, v, var_to_bisect);
   //  cout << " var_to_bisect " << var_to_bisect << endl;
   c.var_to_bisect=var_to_bisect;
   c.ratio=ratio;
  
  }

  void QibexOptimizer::qibex_contract_and_bound(Cell & c){
    Interval& y=c.box[goal_var];
    if (integerobj){
      y=integer(y);
      if (y.is_empty()){
	c.box.set_empty(); return;
      }
    }
    IntervalVector qcp_box(n);
    read_ext_box(c.box,qcp_box);

    double objlb= y.lb();
    if (objlb < -1.e300) objlb=-1.e300;
    //  cout << "c.box " << c.box << endl;
    //    cout << "objlb " << objlb << endl;
    qibex_relaxation_call(qcp_box,objlb);
  
    string status;
  
  
    int n_y_max=n_x+(n_x*(n_x+1))/2;

    Vector v (n_x);
	
    Vector w(n_y_max);
	

    double newlb;
  
    tie (v,w,newlb) = qibex_relaxation_analysis (qcp_box,status);
  
    if (qcp_box.is_empty()) {c.box.set_empty(); return;}
    if  (status== "solved" || status== "'solved?'")
      qibex_bisection_choice(c, qcp_box, v,w );
  
  
    Interval y0=y;
    double ymax=POS_INFINITY;
    if (newlb > NEG_INFINITY){
      ymax=qibex_loupfinder(v);
      if (ymax < POS_INFINITY){
      //cout << " ymax " << ymax <<" newlb " << newlb << endl;
	if (newlb<= ymax){
	  y &= Interval(newlb,ymax);
	//	buffer.contract(ymax);
	}
	else{
	//	y &= Interval(newlb,newlb);}
	  c.box.set_empty(); return;}
      
	//  semble inutile (si le pt faisable est l'optimum, newlb=newub et l'arret de la branche est automatique 
	/*
	if (c.var_to_bisect==-1 && gap0 <= 1.e-10 && loup_finder.integer_check(v))
	  {cout << " leaf found " << endl;
	    y.set_empty();
	  }
	*/

      
	/*
        if(newub > newlb){
	  cout << " error : relaxed bound greater than original objective " << newlb << "  >  " << newub << endl;
	  y.set_empty();}
	*/
	c.prop.update(BoxEvent(c.box,BoxEvent::CONTRACT,BitSet::singleton(n+1,goal_var)));
      }
      else{
	y &= Interval(newlb,POS_INFINITY);
	c.prop.update(BoxEvent(c.box,BoxEvent::CONTRACT,BitSet::singleton(n+1,goal_var)));
      //      cout << " y after newlb" << y << endl;
      }
    

      if (y.is_empty() ) {
	c.box.set_empty();
      }
      else if (recontract && y.diam() < y0.diam()){
	if (integerobj){
	  y=integer(y);
	  if (y.is_empty()){
	    c.box.set_empty(); return;}
	  if (y.lb()==y.ub()){
	    c.box.set_empty(); return;}
	}
	ctc.contract(c.box);
      }
    }
  }

  void QibexOptimizer:: contract(Cell & c)
  {Optimizer::contract(c);
    if (c.box.is_empty()) return;
    qibex_contract_and_bound(c);
  }
  
  double QibexOptimizer::compute_ymax(){
     if (integerobj)
       return loup-1;
     else
       return Optimizer::compute_ymax();
   }
  
  double QibexOptimizer::compute_emptybuffer_uplo(){
     if (integerobj)
       return loup;
     else
       return Optimizer::compute_emptybuffer_uplo();
   }


} // end namespace ibex

