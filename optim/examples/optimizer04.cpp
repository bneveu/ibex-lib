//============================================================================
//                                  I B E X                                   
// File        : optimizer04int.cpp
// Author      : Bertrand Neveu
// Copyright   : Ecole des Mines de Nantes (France)
// License     : See the LICENSE file
// Created     : Jul 12, 2012
// Last Update : Jul 12, 2012
//============================================================================


#include "ibex.h"

#include "ibex_LinearizerAffine2.h"
#include "ibex_AmplInterface.h"



const double default_relax_ratio =0.2;
const double initbox_limit = 1.e8;  // TODO . parameter ??
const double eqeps= 1.e-8;  // TODO parameter ??

using namespace std;
using namespace ibex;

  

int main(int argc, char** argv){
	// ------------------------------------------------
	// Parameterized Optimizer (with a system loaded from a file, and choice of contractor, linearization , bisector, and search strategy)
        // Load a problem to optimize (in format ampl .nl or minibex (.mbx or .bch ))
	// --------------------------
	try {
	  
	if (argc<12) {
		cerr << "usage: optimizer04int filename filtering linear_relaxation bisection upperbounding strategy [beamsize] integerobj prec goal_prec timelimit randomseed"  << endl;
		exit(1);
	}
	  
	System * sys;
	#ifdef __IBEX_AMPL_INTERFACE_H__
	std::size_t found = string(argv[1]).find(".nl");
	if (found!=std::string::npos){
	  AmplInterface interface (argv[1]);
	  sys= new System(interface);
          vector<int> integers = sys->find_integer_variables (argv[1]);
	  sys->set_integer_variables(integers);
	
	}
	else{
	  sys = new System(argv[1]);
	  if (sys->minlp==false){
	    BitSet b (sys->nb_var);
	    sys->set_integer_variables(b);}
	}
	    
	    
	#else
	sys = new System(argv[1]);
	#endif

	
	for (int i=0; i< sys->box.size(); i++){
	  if (sys->box[i].lb() < -initbox_limit) 
	    sys->box[i]= Interval(-initbox_limit,sys->box[i].ub()) ;
	  if (sys->box[i].ub() >initbox_limit)
 	    sys->box[i] =  Interval(sys->box[i].lb(), initbox_limit);
	}


	string filtering = argv[2];
	string linearrelaxation= argv[3];
	string bisection= argv[4];
	string loupfindermethod=argv[5];
	string strategy= argv[6];
	int nbinput=6;
	int beamsize;
	if (strategy=="bs" || strategy== "beamsearch") {beamsize=atoi(argv[6]); nbinput++;}

	int integerobjective= atoi(argv[nbinput+1]);
	double prec= atof(argv[nbinput+2]);
	double goalprec= atof (argv[nbinput+3]);
	double timelimit= atof(argv[nbinput+4]);
	//	double eqeps= 1.e-6;
	double eqeps= 1.e-8;
	int randomseed = atoi(argv[nbinput+5]);
	//	double initloup=atof(argv[nbinput+5]);
	RNG::srand(randomseed);
	//        cout << "fin lecture parametres " << endl;
	// the extended system 
	ExtendedSystem ext_sys(*sys,eqeps);
        NormalizedSystem norm_sys(*sys,eqeps);

        if (sys->minlp)	cout << " integer variables " << *(sys->get_integer_variables()) << endl;

	
	ext_sys.tolerance=eqeps;
	norm_sys.tolerance=eqeps;
	sys->tolerance=eqeps;
	
	//	cout << "nor_sys" << norm_sys << endl;
	//	LoupFinderDefault loupfinder (norm_sys,true);
	LoupFinder* loupfinder;
	if (loupfindermethod=="xninhc4")
	  loupfinder = new LoupFinderDefault (norm_sys,true);
	else if (loupfindermethod=="xn")
	  loupfinder = new LoupFinderDefault (norm_sys,false);
	else if (loupfindermethod=="prob")
	  loupfinder = new LoupFinderProbing (norm_sys);
	else if (loupfindermethod=="inhc4")
	  loupfinder = new LoupFinderInHC4 (norm_sys);
	CellBufferOptim* buffer;
	CellHeap futurebuffer (ext_sys);
       	CellHeap currentbuffer (ext_sys);
	if (strategy=="bfs")
	  buffer = new CellHeap   (ext_sys);
	else if (strategy=="dh")
	  buffer = new CellDoubleHeap  (ext_sys);
       	else if (strategy=="bs")
	  buffer = new CellBeamSearch  (currentbuffer, futurebuffer, ext_sys, beamsize);

	cout << "file " << argv[1] << endl;
	

	cout << " filtering " << filtering; 
        cout << " linearrelaxation " << linearrelaxation;
	cout << " bisection " << bisection ;
	cout << " strategy " << strategy ;
	cout << " randomseed " << randomseed << endl;
	

	// Build the bisection heuristic
	// --------------------------

	Bsc * bs;
	OptimLargestFirst * bs1;

	if  (bisection=="lsmear" || bisection=="smearsum" || bisection=="smearmax" || bisection=="smearsumrel" || bisection=="smearmaxrel" || bisection == "minlpsmearsumrel" ||  bisection == "minlpsmearsum" || bisection=="lsmearmg" || bisection=="lsmearss" || bisection=="lsmearmgss")
	  bs1=  new OptimLargestFirst(ext_sys.goal_var(),true,prec);
        else if
	  (bisection=="lsmearnoobj" || bisection=="smearsumnoobj" || bisection=="smearmaxnoobj" || bisection=="smearsumrelnoobj" || bisection == "minlpsmearsumnoobj" ||  bisection == "minlpsmearsumrelnoobj" || bisection=="smearmaxrelnoobj" || bisection=="lsmearmgnoobj" )
	  bs1=  new OptimLargestFirst(ext_sys.goal_var(),false,prec);

	  
	if (bisection=="roundrobin")
	  bs = new RoundRobin (prec);
	else if (bisection== "largestfirst")
	  bs= new OptimLargestFirst(ext_sys.goal_var(),true,prec);
	else if (bisection== "largestfirstnoobj")
	  bs= new OptimLargestFirst(ext_sys.goal_var(),false,prec);
	else if (bisection== "minlplargestfirst")
	  bs= new MinlpLargestFirst(ext_sys,ext_sys.goal_var(),true,prec);
	else if (bisection== "minlplargestfirstnoobj")
	  bs= new MinlpLargestFirst(ext_sys,ext_sys.goal_var(),false,prec);


	else if (bisection=="smearsum" || bisection== "smearsumnoobj")
	  bs = new SmearSum(ext_sys,prec,*bs1);
	else if (bisection=="smearmax" || bisection == "smearmaxnoobj")
	  bs = new SmearMax(ext_sys,prec,*bs1);
	else if (bisection=="smearsumrel" || bisection=="smearsumrelnoobj")
	  bs = new SmearSumRelative(ext_sys,prec,*bs1);
	else if (bisection=="minlpsmearsumrel" || bisection=="minlpsmearsumrelnoobj")
          bs = new MinlpSmearSumRelative(ext_sys,prec,*bs1);
	else if (bisection=="minlpsmearsum" || bisection=="minlpsmearsumnoobj")
          bs = new MinlpSmearSum(ext_sys,prec,*bs1);

	
	else if (bisection=="smearmaxrel" || bisection=="smearmaxrelnoobj")
	  bs = new SmearMaxRelative(ext_sys,prec,*bs1);
	else if  (bisection=="lsmear" || bisection=="lsmearnoobj"){
	  bs = new LSmear(ext_sys,prec,*bs1,LSMEAR);
	  }
	else if (bisection=="lsmearmg"|| bisection=="lsmearmgnoobj"){
	  bs = new LSmear(ext_sys,prec,*bs1);
	  }
	else {cout << bisection << " is not an implemented  bisection mode "  << endl; return -1;}

	// The contractor
	CtcInteger integ (ext_sys.nb_var,*(ext_sys.get_integer_variables()));

	// the first contractor called
	CtcHC4 hc4(ext_sys.ctrs,0.01,true);
	CtcCompo hc4integ (integ, hc4, integ);
	// hc4 inside acid and 3bcid : incremental propagation beginning with the shaved variable
	CtcHC4 hc44cid(ext_sys.ctrs,0.1,true);
	//	CtcHC4 hc44cid(ext_sys.ctrs,0.01,true);
	// hc4 inside xnewton loop 
	CtcHC4 hc44xn (ext_sys.ctrs,0.01,false);

	// The 3BCID contractor on all variables (component of the contractor when filtering == "3bcidhc4") 
	Ctc3BCid c3bcidhc4(hc44cid);
	// hc4 followed by 3bcidhc4 : the actual contractor used when filtering == "3bcidhc4" 
	CtcCompo hc43bcidhc4 (integ,hc4, integ, c3bcidhc4, integ);

	// The ACID contractor (component of the contractor  when filtering == "acidhc4")
	CtcAcid acidhc4(ext_sys,hc44cid,true);
	// hc4 followed by acidhc4 : the actual contractor used when filtering == "acidhc4" 
	CtcCompo hc4acidhc4 (integ, hc4, integ, acidhc4, integ);

      

	Ctc* ctc;
	if (filtering == "hc4")
	  ctc= &hc4integ;
	else if
	  (filtering =="acidhc4")   
	  ctc= &hc4acidhc4;
	else if 
	  (filtering =="3bcidhc4")
	  ctc= &hc43bcidhc4;
	else {cout << filtering <<  " is not an implemented  contraction  mode "  << endl; return -1;}

	Linearizer* lr;
	Linearizer* lr1;



	if (linearrelaxation=="art")
	  lr= new LinearizerAffine2(ext_sys);
	else if  (linearrelaxation=="compo")
	  lr= new LinearizerCompo( *(new LinearizerXTaylor( ext_sys)),
				   *(new LinearizerAffine2(ext_sys)));

	else if (linearrelaxation=="xn")
	  lr= new LinearizerXTaylor (ext_sys);
	else if (linearrelaxation=="xnart"){
	  lr=new LinearizerXTaylor (ext_sys);
	  lr1=new LinearizerAffine2(ext_sys);
	}
	  
	//	else {cout << linearrelaxation  <<  " is not an implemented  linear relaxation mode "  << endl; return -1;}
	// fixpoint linear relaxation , hc4  with default fix point ratio 0.2
	//	CtcFixPoint* cxn;
	Ctc* cxn;
	CtcPolytopeHull* cxn_poly;
	CtcPolytopeHull* cxn_poly1;
	CtcCompo* cxn_compo;
	if (linearrelaxation=="compo" || linearrelaxation=="art"|| linearrelaxation=="xn")
          {
		cxn_poly = new CtcPolytopeHull(*lr);
		cxn_compo =new CtcCompo(integ,*cxn_poly,integ, hc44xn,integ);
		cxn = new CtcFixPoint (*cxn_compo, default_relax_ratio);
		//cxn =new CtcCompo(*cxn_poly, hc44xn);
	  }
	else if  (linearrelaxation=="xnart")
	  {
	    cxn_poly = new CtcPolytopeHull(*lr);
	    cxn_poly1 = new CtcPolytopeHull(*lr1);
	    cxn_compo =new CtcCompo(integ, *cxn_poly1, *cxn_poly, hc44xn, integ);
	    cxn = new CtcFixPoint (*cxn_compo, default_relax_ratio);
	    
	    cxn =new CtcCompo(*cxn_poly1, *cxn_poly, hc44xn);

	  }

	//  the actual contractor  ctc + linear relaxation 
	Ctc* ctcxn;
	if (linearrelaxation=="compo" || linearrelaxation=="art"|| linearrelaxation=="xn" || linearrelaxation=="xnart") 
          ctcxn= new CtcCompo  (*ctc, *cxn, integ); 
	
	else
	  ctcxn = ctc;
	/*
	Ctc* ctckkt = new CtcKhunTucker(norm_sys, true);
	ctcxn = new CtcCompo (*ctcxn , *ctckkt);
	*/

	// the optimizer : the same precision goalprec is used as relative and absolute precision
	Optimizer o(sys->nb_var,*ctcxn,*bs,*loupfinder,*buffer,ext_sys.goal_var(),prec,goalprec,goalprec);
       	cout << " sys.box " << sys->box << endl;

	// the trace 
	o.trace=1;

	//integer objective
	o.integerobj=integerobjective;
	// the allowed time for search
	o.timeout=timelimit;
	cout.precision(16);
	// the search itself 
	//	o.optimize(sys->box,initloup);
	//	std::ofstream Out("err.txt");
	//	std::streambuf* OldBuf = std::cerr.rdbuf(Out.rdbuf());
	o.optimize(sys->box);
	//	std::cerr.rdbuf(OldBuf);

	// printing the results     
	o.report();
        cout << o.get_time() << "  " << o.get_nb_cells()+1 << endl;

	//	if (filtering == "acidhc4"  )
	//cout    << " nbcidvar " <<  acidhc4.nbvar_stat() << endl;

	delete bs;

	if  (bisection=="lsmear" || bisection=="smearsum" || bisection=="smearmax" || bisection=="smearsumrel" || bisection=="smearmaxrel" || bisection=="lsmearmg" ||bisection =="minlpsmearsumrel" ||bisection =="minlpsmearsum" || bisection=="lsmearnoobj" || bisection=="smearsumnoobj" || bisection=="smearmaxnoobj" || bisection=="smearsumrelnoobj" || bisection =="minlpsmearsumnoobj" || bisection == "minlpsmearsumrelnoobj" || bisection=="smearmaxrelnoobj" || bisection=="lsmearmgnoobj" )

	  delete bs1;
	// bs1 deleted by SmearFunction destructor  (TO CHANGE)
	
	delete loupfinder;

	delete buffer;
	if (linearrelaxation=="compo" || linearrelaxation=="art"|| linearrelaxation=="xn" || linearrelaxation=="xnart") {
		delete lr;
	    delete ctcxn;
	    delete cxn;
	    delete cxn_poly;
	    delete cxn_compo;
	}
	if (linearrelaxation=="xnart"){
	  delete lr1;
	  delete cxn_poly1;
	}
	delete sys;
	return 0;
	
	}


	catch(ibex::SyntaxError& e) {
	  cout << e << endl;
	}
}
