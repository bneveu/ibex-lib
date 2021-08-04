//                                  I B E X
// File        : ibex_QibexOptimizer.cpp
// Author      : Gilles Chabert, Bertrand Neveu
// Copyright   : IMT Atlantique (France)
// License     : See the LICENSE file
// Created     : May 14, 2012
// Last Update : Apr 08, 2019
//============================================================================

#include "ibex_QibexOptimizer.h"

#include <float.h>
#include <stdlib.h>



using namespace std;

namespace ibex {

  QibexOptimizer::QibexOptimizer(int n, Ctc& ctc, Bsc& bsc, LoupFinder& finder,
		CellBufferOptim& buffer,
				 int goal_var, double eps_x, double rel_eps_f, double abs_eps_f) : Optimizer (n,ctc,bsc,finder,buffer,goal_var,eps_x,rel_eps_f,abs_eps_f){;}

  

  pair<Vector,double> QibexOptimizer::qibex_newbounds(const IntervalVector & box){
    //    cout << " box before qibex " << box << endl;
	ofstream fic ("bound.ampl", ofstream::trunc);
	fic << "data;" << endl;
	fic << "param ub :=" << endl;
	
	for (int i=0; i< box.size(); i++)
	  fic << i+1 << " " << box[i].ub() << endl;
	fic << ";" <<endl;
	fic << "param lb :=" << endl; 
	for (int i=0; i< box.size(); i++)
	  fic << i+1 << " " << box[i].lb() << endl;
	fic << ";" << endl;
	fic << "end ;";
	fic.close();

	double newlb = NEG_INFINITY;
	Vector v(box.size());

	system("~/imagine4/RECHERCHE/ampl/ampl model_quad_relax.run");


	ifstream fic1 ("lower_bound.txt");
	if (fic1.good()){
	    string a;
            fic1 >> a; fic1 >>a; fic1 >> a;
	    if (a=="solved"){
	      fic1 >> a;
	      fic1>> a;
	      fic1 >> newlb;

	      ifstream fic2 ("xval.txt");
	      fic2 >> a;
	      fic2>> a;
	      fic2>> a;
	      for (int i =0 ; i< box.size(); i++){
		fic2 >> a;
		fic2 >> v[i];
	      }
	      fic2 >> a;
	      fic2.close();
	    }
	    fic1.close();
	}
        pair<Vector,double> p(v,newlb);
	return p;
  }
  


 
void QibexOptimizer::qibex_contract(Interval& y, Cell & c){
  IntervalVector qcp_box(n);
  read_ext_box(c.box,qcp_box);

  pair<Vector,double> vecnewbounds  = qibex_newbounds(qcp_box);
  Vector v= vecnewbounds.first;
  double newlb=vecnewbounds.second;
  if (newlb > NEG_INFINITY){
    double newub=POS_INFINITY;
    if (loup_finder.integer_check(v)){
      newub=loup_finder.goal_ub(v); 
      if (newub < loup && loup_finder.is_inner(v)){
	loup= newub;
	if (trace) {
	  cout << "                    ";
	  cout << "*** \033[32m loup= " << loup << "\033[0m" << endl;
	      //	      cout << "*** \033[32m loup_point= " << v << "\033[0m" << endl;
	}
	loup_point = v;
	loup_changed=true;
	double ymax=compute_ymax();
	if (newlb<= ymax)
	  y &= Interval(newlb,ymax);
	else if (newlb <= newub)
	  //	  y &= Interval(ymax ,ymax);
	  y.set_empty();
	else{
	  cout << " error : relaxed bound greater than original objective " << newlb << "  >  " << newub << endl;
	  y.set_empty();}
	c.prop.update(BoxEvent(c.box,BoxEvent::CONTRACT,BitSet::singleton(n+1,goal_var)));
      }
      else{
	y &= Interval(newlb,POS_INFINITY);
	c.prop.update(BoxEvent(c.box,BoxEvent::CONTRACT,BitSet::singleton(n+1,goal_var)));
      }
    }
    
    else {
      y &= Interval(newlb,POS_INFINITY);
      c.prop.update(BoxEvent(c.box,BoxEvent::CONTRACT,BitSet::singleton(n+1,goal_var)));
    }
  }
}


} // end namespace ibex

