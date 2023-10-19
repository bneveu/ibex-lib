//============================================================================
//                                  I B E X                                   
// File        : optimizer06.cpp
// Author      : Bertrand Neveu
// Copyright   : Ecole des Mines de Nantes (France)
// License     : See the LICENSE file
// Created     : Aug 2023
// Last Update : Sept 27, 2023
//============================================================================


#include "ibex.h"

#include "ibex_LinearizerAffine2.h"
#include "ibex_AmplInterface.h"



const double default_relax_ratio =0.2;
const double initbox_limit = 1.e8;  // TODO . parameter ??

using namespace std;
using namespace ibex;



int main(int argc, char** argv){
	// ------------------------------------------------
	// Parameterized Optimizer (with a system loaded from a file, and choice of contractor, linearization , bisector, and search strategy)
        // Load a problem to optimize (in format ampl .nl or minibex (.mbx or .bch ))
	// --------------------------
	try {
	  
	if (argc<13) {
		cerr << "usage: optimizer04int filename filtering linear_relaxation bisection upperbounding strategy [beamsize] integerobj prec goal_prec tolerance timelimit randomseed"  << endl;
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
	}
	    
	    
	#else
	sys = new System(argv[1]);
	#endif

	if (!(sys->goal)) {cout << " No goal " << endl; return -1;}
	for (int i=0; i< sys->box.size(); i++){
	  if (sys->box[i].lb() < -initbox_limit) 
	    sys->box[i]= Interval(-initbox_limit, sys->box[i].ub()) ;
	  if (sys->box[i].ub() >initbox_limit)
 	    sys->box[i] = Interval(sys->box[i].lb(), initbox_limit);
	}

	string filtering = argv[2];
	string linearrelaxation= argv[3];
	string bisection= argv[4];
	string loupfindermethod=argv[5];
	string strategy= argv[6];
	int nbinput=6;
	int beamsize;
	if (strategy=="bs" || strategy== "beamsearch") {beamsize=atoi(argv[7]); nbinput++;}

	int integerobjective= atoi(argv[nbinput+1]);
	double prec= atof(argv[nbinput+2]);
	double goalprec= atof (argv[nbinput+3]);
	double tolerance= atof (argv[nbinput+4]);
	double timelimit= atof(argv[nbinput+5]);

	int randomseed = atoi(argv[nbinput+6]);
	
	//	double initloup=atof(argv[nbinput+7]);
	RNG::srand(randomseed);
	//        cout << "fin lecture parametres " << endl;

	
	//	if (sys->minlp)	cout << " number of integer variables " << (sys->get_integer_variables())->size() << endl;
	if (sys->minlp)	cout << " integer variables " << *(sys->get_integer_variables()) << endl;
	// the extended system
	ExtendedSystem ext_sys(*sys,tolerance);
        NormalizedSystem norm_sys(*sys,tolerance);

	//	ext_sys.tolerance=tolerance;
	//	norm_sys.tolerance=tolerance;
	//	sys->tolerance=tolerance;
	//	cout << " sys " << endl;
	//	cout << *sys << endl;
	
	//	cout << norm_sys << endl;

	//	cout << "ext_sys" << ext_sys << endl;
	
	//	LoupFinderDefault loupfinder (norm_sys,true);
	LoupFinder* loupfinder;
	if (loupfindermethod=="ipoptxninhc4")
	  loupfinder = new LoupFinderDefaultIpopt (*sys,norm_sys,ext_sys,true,integerobjective);
	else if (loupfindermethod=="ipoptxn")
	  loupfinder = new LoupFinderDefaultIpopt (*sys,norm_sys,ext_sys,false,integerobjective);
	else if (loupfindermethod=="xninhc4")
	  loupfinder = new LoupFinderDefault (norm_sys,true,integerobjective);
	else if (loupfindermethod=="xn")
	  loupfinder = new LoupFinderDefault (norm_sys,false,integerobjective);
	else if (loupfindermethod=="prob")
	  loupfinder = new LoupFinderProbing (norm_sys);
	else if (loupfindermethod=="inhc4")
	  loupfinder = new LoupFinderInHC4 (norm_sys);
	else
	  {cout << loupfindermethod <<  " is not an implemented  feasible point finding method "  << endl; return -1;}

	CellBufferOptim* buffer;
	CellHeap futurebuffer (ext_sys);
       	CellHeap currentbuffer (ext_sys);
	if (strategy=="bfs")
	  buffer = new CellHeap   (ext_sys);
	else if (strategy=="dh")
	  buffer = new CellDoubleHeap  (ext_sys);
       	else if (strategy=="bs")
	  buffer = new CellBeamSearch  (currentbuffer, futurebuffer, ext_sys, beamsize);
	else
	  {cout << strategy <<  " is not an implemented  node selection strategy "  << endl; return -1;}
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

	if  (bisection=="lsmear" || bisection=="smearsum" || bisection=="smearmax" || bisection=="smearsumrel" || bisection=="smearmaxrel" || bisection == "minlpsmearsumrel" ||  bisection == "minlpsmearsum" || bisection=="lsmearmg" || bisection=="lsmearss" || bisection=="lsmearmgss" || bisection== "minlplsmear" || bisection== "minlplsmearmg")
	  bs1=  new OptimLargestFirst(ext_sys.goal_var(),true,prec);
        else if
	  (bisection=="lsmearnoobj" || bisection=="smearsumnoobj" || bisection=="smearmaxnoobj" || bisection=="smearsumrelnoobj" || bisection == "minlpsmearsumnoobj" ||  bisection == "minlpsmearsumrelnoobj" || bisection=="smearmaxrelnoobj" || bisection=="lsmearmgnoobj" || bisection == "minlplsmearmgnoobj" || bisection == "minlplsmearnoobj" )
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


	else if (bisection=="smearsum") 
	  bs = new SmearSum(ext_sys,prec,*bs1,true);
	else if (bisection== "smearsumnoobj")
          bs = new SmearSum(ext_sys,prec,*bs1,false);
	else if (bisection=="smearmax")
	  bs = new SmearMax(ext_sys,prec,*bs1,true);
	else if (bisection == "smearmaxnoobj")
	  bs = new SmearMax(ext_sys,prec,*bs1,false);
	else if (bisection=="smearsumrel")
	  bs = new SmearSumRelative(ext_sys,prec,*bs1,true);
	else if ( bisection=="smearsumrelnoobj")
	  bs = new SmearSumRelative(ext_sys,prec,*bs1,false);
	else if (bisection=="minlpsmearsumrel")
          bs = new MinlpSmearSumRelative(ext_sys,prec,*bs1,true);
	else if (bisection=="minlpsmearsumrelnoobj")
	  bs = new MinlpSmearSumRelative(ext_sys,prec,*bs1,false);
	else if (bisection=="minlpsmearsum")
          bs = new MinlpSmearSum(ext_sys,prec,*bs1,true);
	else if  (bisection=="minlpsmearsumnoobj")
	  bs = new MinlpSmearSum(ext_sys,prec,*bs1,false);
	
	
	else if (bisection=="smearmaxrel") 
	  bs = new SmearMaxRelative(ext_sys,prec,*bs1,true);
	else if (bisection=="smearmaxrelnoobj")
	  bs = new SmearMaxRelative(ext_sys,prec,*bs1,false);
	
	else if  (bisection=="lsmear" || bisection=="lsmearnoobj"){
	  bs = new LSmear(ext_sys,prec,*bs1,LSMEAR);
	  }
	else if (bisection=="lsmearmg"|| bisection=="lsmearmgnoobj"){
	  bs = new LSmear(ext_sys,prec,*bs1);
	  }


	
	else if (bisection=="minlplsmearmg")
	  bs = new MinlpLSmear(ext_sys,prec,*bs1,true );
        else if  (bisection=="minlplsmearmgnoobj")
	  bs = new MinlpLSmear(ext_sys,prec,*bs1,false );
	else if  (bisection=="minlplsmear")
	  bs = new MinlpLSmear(ext_sys,prec,*bs1,true ,MINLPLSMEAR);
	else if  (bisection=="minlplsmearnoobj")
	  bs = new MinlpLSmear(ext_sys,prec,*bs1,false ,MINLPLSMEAR);
	  
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
	else if (linearrelaxation=="no") {;}
	else {cout << linearrelaxation  <<  " is not an implemented  linear relaxation mode "  << endl; return -1;}

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
		//	cxn=cxn_compo;
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
	if (sys->nb_ctr==0 && sys->minlp==false){  // CtcKuhnTucker for unconstrained continuous problems
	  Ctc* ctckkt = new CtcKuhnTucker(norm_sys, true);
	  ctcxn = new CtcCompo (*ctcxn , *ctckkt, integ);
	}


	// the optimizer : the same precision goalprec is used as relative and absolute precision
	Optimizer o(sys->nb_var,*ctcxn,*bs,*loupfinder,*buffer,ext_sys.goal_var(),prec,goalprec,goalprec);

      
	// the trace 
	o.trace=1;
	//	if (o.trace) cout << " sys.box " << sys->box << endl;
	cout.precision(16);
	//integer objective
	loupfinder->integerobj=integerobjective;
	o.integerobj=integerobjective;
	//	o.integer_tolerance=goalprec;

	// ipopt preprocessing

	if (loupfindermethod=="ipoptxninhc4" || loupfindermethod=="ipoptxn")
	  ((LoupFinderDefaultIpopt*) loupfinder)->finder_ipopt.optimizer= &o;
	o.ipopt_preprocessing(*sys,norm_sys,ext_sys);
	if (loupfindermethod=="ipoptxn" ||loupfindermethod=="ipoptxninhc4" )
	  ((LoupFinderDefaultIpopt*)loupfinder)->finder_ipopt.recursive_call=true;
       
	// the allowed time for search
	o.timeout=timelimit;

	// the search itself 
	//	o.optimize(sys->box,initloup);
	//	std::ofstream Out("err.txt");
	//	std::streambuf* OldBuf = std::cerr.rdbuf(Out.rdbuf());
	if (o.trace) cout << " sys.box " << sys->box << endl;
	
	o.optimize(sys->box,o.get_loup());
	//	std::cerr.rdbuf(OldBuf);

	// printing the results     
	o.report();
	cout << "preprocessing ampl time " << o.preprocampltime << " ipopttime " << o.preprocipopttime << endl;
        cout << o.get_time() << "  " << o.get_nb_cells()+1 << endl;
        if (loupfindermethod == "ipoptxn" || loupfindermethod =="ipoptxninhc4"){
	  cout << " ipopttime " << ((LoupFinderDefaultIpopt*)loupfinder)->finder_ipopt.ipopttime << " ampltime " << ((LoupFinderDefaultIpopt*)loupfinder)->finder_ipopt.ampltime << endl;
	  cout << " correction nodes " << ((LoupFinderDefaultIpopt*)loupfinder)->finder_ipopt.correction_nodes << " correction time " << ((LoupFinderDefaultIpopt*)loupfinder)->finder_ipopt.correction_time << endl;
	}
	//	if (filtering == "acidhc4"  )
	//cout    << " nbcidvar " <<  acidhc4.nbvar_stat() << endl;

	delete bs;

	if  (bisection=="lsmear" || bisection=="smearsum" || bisection=="smearmax" || bisection=="smearsumrel" || bisection=="smearmaxrel" || bisection=="lsmearmg" ||bisection =="minlpsmearsumrel" ||bisection =="minlpsmearsum" || bisection=="lsmearnoobj" || bisection=="smearsumnoobj" || bisection=="smearmaxnoobj" || bisection=="smearsumrelnoobj" || bisection =="minlpsmearsumnoobj" || bisection == "minlpsmearsumrelnoobj" || bisection=="smearmaxrelnoobj" || bisection=="lsmearmgnoobj" )

	  delete bs1;
	
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
