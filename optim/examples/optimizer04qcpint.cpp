//============================================================================
//                                  I B E X                                   
// File        : optimizer04qcpeqint.cpp
// Author      : Bertrand Neveu
// Copyright   : Ecole des Mines de Nantes (France)
// License     : See the LICENSE file
// Created     : Jul 12, 2012
// Last Update : Jul 12, 2012
//============================================================================


#include "ibex.h"

#include "ibex_LinearizerAffine2.h"
#include "ibex_AmplInterface.h"

const double default_relax_ratio = 0.2;
//const double default_relax_ratio = 0.01;

using namespace std;
using namespace ibex;



int main(int argc, char** argv){
	// ------------------------------------------------
	// Parameterized Optimizer (with a system loaded from a file, and choice of contractor, linearization , bisector, and search strategy)
        // Load a problem to optimize (in format ampl .nl or minibex (.mbx or .bch ))
	// --------------------------
	try {
Timer timer;
	timer.start();
	if (argc<16) {
		cerr << "usage: optimizer04 filename integerfile filtering linear_relaxation bisection upperbounding strategy [beamsize] recontract rigor qibexwidth integerobj prec goal_prec tolerance timelimit randomseed"  << endl;
		exit(1);
	}

	System * sys;
	NormalizedSystem * norm_sys;
	#ifdef __IBEX_AMPL_INTERFACE_H__
	std::size_t found = string(argv[1]).find(".nl");
	if (found!=std::string::npos){
	  AmplInterface interface (argv[1]);
	  sys= new System(interface);
	}else
	  sys = new System(argv[1],0);
	#else
	sys = new System(argv[1],0);
	#endif

	
	
	for (int i=0; i< sys->box.size(); i++){
	  if (sys->box[i].lb() < -1.e8) 
	    sys->box[i]= Interval(-1.e8,sys->box[i].ub()) ;
	  if (sys->box[i].ub() >1.e8) 
	    sys->box[i] =  Interval(sys->box[i].lb(), 1.e8);
	}
	
	string integerfile = argv[2];
	ifstream ficint(integerfile);
	int nbint;
	ficint >> nbint;
	vector<int> l;
	int k;
	for (int i=0;i<nbint; i++){
	  ficint >> k;
	  l.push_back(k);}
        ficint.close();
	string filtering = argv[3];
	string linearrelaxation= argv[4];
	string bisection= argv[5];
	string loupfind=argv[6];
	string strategy= argv[7];
	int nbinput=7;
	int beamsize;
	if (strategy=="bs" || strategy== "beamsearch") {beamsize=atoi(argv[7]); nbinput++;}
	
	string recontraction= argv[nbinput+1];
	string rigormode = argv[nbinput+2];
	double qibexwidth= atof(argv[nbinput+3]);
	int integerobjective= atoi(argv[nbinput+4]);
	double prec= atof(argv[nbinput+5]);
	double goalprec= atof (argv[nbinput+6]);
	double tolerance=atof(argv[nbinput+7]);
	double timelimit= atof(argv[nbinput+8]);
	//	double eqeps= 1.e-6;
	double eqeps= tolerance;
	int randomseed = atoi(argv[nbinput+9]);
	//	double initloup=atof(argv[nbinput+10]);
	RNG::srand(randomseed);

	// the extended system 
	ExtendedSystem ext_sys(*sys,eqeps);
	// no need of normalized system if all constraints are LEQ
	int leq=0;
	for (int j=0; j < sys->nb_ctr; j++){
	  if (sys->ops[j] == LEQ)
	    leq++;
	}
	//	if (leq==sys->nb_ctr)  cout << "only leq " << endl;
	if (leq==sys->nb_ctr)
	  norm_sys=(NormalizedSystem*)sys;
	else
	  norm_sys= new NormalizedSystem (*sys,eqeps);

	ext_sys.tolerance=tolerance;
	norm_sys->tolerance=tolerance;
	sys->tolerance=tolerance;
	//	cout << "nor_sys" << norm_sys << endl;
	BitSet b (ext_sys.nb_var);
	for (int i=0; i<l.size();i++)
	  b.add(l[i]);
	//	sys->minlp=true;
	//	sys->set_integer_variables(b);
	ext_sys.set_integer_variables(b);
	ext_sys.minlp=true;
	//	cout << " integer variables " << *(ext_sys.get_integer_variables()) << endl;
	norm_sys->minlp=true;
	norm_sys->set_integer_variables(b);

	//	cout << "loupfind " << loupfind << endl;
	LoupFinder* loupfinder;
	if (loupfind=="xninhc4")
	  loupfinder= new LoupFinderDefault (*norm_sys,true);
	else if (loupfind=="xn")
	  loupfinder= new LoupFinderDefault (*norm_sys,false);
	  
	else if (loupfind=="inhc4")
	  loupfinder= new LoupFinderInHC4 (*norm_sys);
	else if (loupfind=="prob" || loupfind=="no")
	  loupfinder= new LoupFinderProbing (*norm_sys);
	  
	else 
	  {cout << " loupfinder not found " << endl;  return(-1) ;}

	CellBufferOptim* buffer;
	CellHeap futurebuffer (ext_sys);
	CellHeap currentbuffer (ext_sys);
	
       	if (strategy=="bfs")
	  buffer = new CellHeap   (ext_sys);
	else  if (strategy=="dh")
	  buffer = new CellDoubleHeap  (ext_sys);
	else if (strategy=="bs"){
	  buffer = new CellBeamSearch  (currentbuffer, futurebuffer, ext_sys, beamsize);
	}
	cout << "file " << argv[1] << endl;


	//	cout << " filtering " << filtering; 
	//      cout << " linearrelaxation " << linearrelaxation;
	//	cout << " bisection " << bisection ;
	//	cout << " strategy " << strategy ;
	//	cout << " randomseed " << randomseed << endl;


	// Build the bisection heuristic
	// --------------------------

	Bsc * bs;
	OptimLargestFirst * bs1;

	if  (bisection=="lsmear" || bisection=="smearsum" || bisection=="smearmax" || bisection=="smearsumrel" || bisection=="smearmaxrel" || bisection == "minlpsmearsumrel" ||  bisection == "minlpsmearsum" || bisection=="lsmearmg" || bisection=="lsmearss" || bisection=="lsmearmgss" || bisection=="qibexsmearsumrel" ||bisection=="qibexsmearsum" ||bisection=="qibexlargestfirst" )
	  bs1=  new OptimLargestFirst(ext_sys.goal_var(),true,prec);
	if  (bisection=="lsmearnoobj" || bisection=="lsmearmgnoobj" || bisection=="smearsumrelnoobj"|| bisection=="smearsumnoobj" || bisection == "minlpsmearsumnoobj" ||  bisection == "minlpsmearsumrelnoobj" || bisection=="qibexsmearsumrelnoobj" || bisection=="qibexsmearsumnoobj" || bisection == "qibexlargestfirstnoobj" )
	  bs1=  new OptimLargestFirst(ext_sys.goal_var(),false,prec);
	if (bisection=="roundrobin")
	  bs = new RoundRobin (prec,0.5);
	else if (bisection=="qibexroundrobin")
	  bs = new QibexRoundRobin (prec,0.5);
	else if (bisection== "largestfirst")
	  bs= new OptimLargestFirst(ext_sys.goal_var(),true,prec,0.5);
	else if (bisection== "largestfirstnoobj")
	  bs= new OptimLargestFirst(ext_sys.goal_var(),false,prec,0.5);
	else if (bisection== "qibexlargestfirst")
	  bs= new QibexLargestFirst(ext_sys.goal_var(),true,prec,0.5);
	else if (bisection== "qibexlargestfirstnoobj")
	  bs= new QibexLargestFirst(ext_sys.goal_var(),false,prec,0.5);
	else if (bisection=="smearsum" || bisection=="smearsumnoobj")
	  bs = new SmearSum(ext_sys,prec,*bs1);
	else if (bisection=="smearmax")
	  bs = new SmearMax(ext_sys,prec,*bs1);
	else if (bisection=="smearsumrel" || bisection=="smearsumrelnoobj")
	  bs = new SmearSumRelative(ext_sys,prec,*bs1);
	else if (bisection=="minlpsmearsumrel" || bisection=="minlpsmearsumrelnoobj")
          bs = new MinlpSmearSumRelative(ext_sys,prec,*bs1);
	else if (bisection=="minlpsmearsum" || bisection=="minlpsmearsumnoobj")
          bs = new MinlpSmearSum(ext_sys,prec,*bs1);

	
	else if (bisection=="smearmaxrel")
	  bs = new SmearMaxRelative(ext_sys,prec,*bs1);
	else if  (bisection=="lsmear" || bisection=="lsmearnoobj")
	  bs = new LSmear(ext_sys,prec,*bs1,LSMEAR);
	else if  (bisection=="qibexsmearsumrel" || bisection=="qibexsmearsumrelnoobj")
	  bs = new QibexSmearSumRelative(ext_sys,prec,*bs1);
	else if  (bisection=="qibexsmearsum" || bisection=="qibexsmearsumnoobj")
	  bs = new QibexSmearSum(ext_sys,prec,*bs1);
	else if (bisection=="lsmearmg"|| bisection=="lsmearmgnoobj" )
	  bs = new LSmear(ext_sys,prec,*bs1);
	else {cout << bisection << " is not an implemented  bisection mode "  << endl; return -1;}

	// The contractors
	CtcInteger integ (ext_sys.nb_var,*(ext_sys.get_integer_variables()));
	// the first contractor called
	CtcHC4 hc4(ext_sys.ctrs,0.01,true);
	CtcCompo hc4integ (integ, hc4, integ);
	// hc4 inside acid and 3bcid : incremental propagation beginning with the shaved variable
	CtcHC4 hc44cid(ext_sys.ctrs,0.1,true);
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
	else if
	   (filtering =="no")
	  ctc=&integ;
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
		cxn_compo =new CtcCompo(integ,*cxn_poly, integ, hc44xn, integ);
		cxn = new CtcFixPoint (*cxn_compo, default_relax_ratio);
		//	cxn =new CtcCompo(*cxn_poly, hc44xn);
	  }
	else if  (linearrelaxation=="xnart")
	  {
	    //	    cout << " xnart " << endl;
	    cxn_poly = new CtcPolytopeHull(*lr);
	    cxn_poly1 = new CtcPolytopeHull(*lr1);
	    cxn_compo =new CtcCompo(integ,*cxn_poly1, *cxn_poly, hc44xn, integ);
	    //	    cxn = new CtcFixPoint (*cxn_compo, default_relax_ratio);
	    
	    cxn =new CtcCompo(*cxn_poly1, *cxn_poly, hc44xn);
	  }

	//  the actual contractor  ctc + linear relaxation 
	Ctc* ctcxn;
	if (linearrelaxation=="compo" || linearrelaxation=="art"|| linearrelaxation=="xn" || linearrelaxation=="xnart") 
          ctcxn= new CtcCompo  (*ctc, *cxn,integ); 
	
	else
	  ctcxn = ctc;

	IntervalVector boxn (sys->nb_var);
	IntervalVector gradient(sys->nb_var);
	IntervalVector gradient1(sys->nb_var);
	for (int i=0; i < sys->nb_var; i++)
	  boxn[i]=Interval(0.0,0.0);
	gradient=  sys->goal->gradient(boxn);
	ofstream fic("diag_hessian.txt");
	for (int i=0; i < sys->nb_var; i++){
	  boxn[i]=Interval(1.0,1.0);
	  gradient1=sys->goal->gradient(boxn);
	  fic << ((gradient1[i]-gradient[i])/2).mid() << " " ;
	  boxn[i]=Interval(0.0,0.0);
	}
	fic.close();
	// the optimizer : the same precision goalprec is used as relative and absolute precision
	QibexOptimizer o(sys->nb_var,*ctcxn,*bs,*loupfinder,*buffer,ext_sys.goal_var(),qibexwidth,tolerance,prec,goalprec,goalprec);

	//	cout << " sys.box " << sys->box << endl;

	// the trace 
	o.trace=1;
	// rigor mode
	if (rigormode=="r")
	  o.rigor=true;
	else
	  o.rigor=false;

       	// contraction after relaxation
	if (recontraction=="r")
	  o.recontract=true;
	else
	  o.recontract=false;

	if (loupfind=="no") o.loupfinderp=false;

        //integer objective
	o.integerobj=integerobjective;

	// the allowed time for search
	o.timeout=timelimit;
	//	cout << " timelimit " << timelimit << endl;
	cout.precision(16);

	// the search itself 
	//	o.optimize(sys->box,initloup);
	//	std::ofstream Out("err.txt");
	//	std::streambuf* OldBuf = std::cerr.rdbuf(Out.rdbuf());
	timer.stop();
	double time = timer.get_time();
	cout << " presolve time " << time << endl;
	
	o.optimize(sys->box);
	//	std::cerr.rdbuf(OldBuf);

	// printing the results     
	o.report();
        cout << o.get_time() << "  " << o.get_nb_cells()+1 << endl;
	cout << "external solver time " << o.solvertime << " ampl time " << o.ampltime << endl;
	/*
	if (filtering == "acidhc4"  )
	cout    << " nbcidvar " <<  acidhc4.nbvar_stat() << endl;
	*/
	delete bs;

	if  (bisection=="lsmear" || bisection=="smearsum" || bisection=="smearmax" || bisection=="smearsumrel" || bisection=="smearmaxrel" || bisection=="lsmearmg" || bisection =="minlpsmearsumrel" || bisection =="minlpsmearsum" || bisection=="qibexsmearsumrel" || bisection=="lsmearnoobj" || bisection=="lsmearmgnoobj" || bisection=="smearsumrelnoobj"|| bisection=="smearsumnoobj"|| bisection =="minlpsmearsumnoobj"|| bisection == "minlpsmearsumrelnoobj"|| bisection=="qibexsmearsumrelnoobj" || bisection == "qibexsmearsumnoobj" || bisection== "qibexlargestfirstnoobj")
	  delete bs1;
	

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
	if (leq <sys->nb_ctr) delete norm_sys;
	
	}

	catch(ibex::SyntaxError& e) {
	  cout << e << endl;
	}
}
