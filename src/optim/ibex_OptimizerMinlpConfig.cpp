//============================================================================
//                                  I B E X                                   
// File        : ibex_OptimizerMinlpConfig.cpp
// Author      : Bertrand Neveu, Gilles Chabert
// Copyright   : IMT Atlantique (France)
// License     : See the LICENSE file
// Created     : Dec 11, 2014
// Last Update : Oct 15, 2019
//============================================================================

#include "ibex_OptimizerMinlpConfig.h"
#include "ibex_CtcHC4.h"
#include "ibex_CtcAcid.h"
#include "ibex_Ctc3BCid.h"
#include "ibex_CtcPolytopeHull.h"
#include "ibex_CtcCompo.h"
#include "ibex_CtcInteger.h"
#include "ibex_CtcNewton.h"
#include "ibex_CtcFixPoint.h"
#include "ibex_CtcKuhnTucker.h"
#include "ibex_Exception.h"
#include "ibex_OptimLargestFirst.h"
#include "ibex_RoundRobin.h"
#include "ibex_SmearFunction.h"
#include "ibex_MinlpSmearSumRelative.h"
#include "ibex_MinlpSmearSum.h"
#include "ibex_LSmear.h"
#include "ibex_Random.h"
#include "ibex_NormalizedSystem.h"
#include "ibex_LinearizerXTaylor.h"
#include "ibex_CellHeap.h"
#include "ibex_CellDoubleHeap.h"
#include "ibex_CellBeamSearch.h"
#include "ibex_LoupFinderDefault.h"
#include "ibex_SyntaxError.h"

#include <sstream>
#include <vector>

using namespace std;

namespace ibex {

  OptimizerMinlpConfig::OptimizerMinlpConfig() {;}
  
OptimizerMinlpConfig::OptimizerMinlpConfig(int argc, char** argv) {

	try {
		if (argc<9) {
			ibex_error("usage: optimizerMinlp filename integerfilename filtering linear_relaxation bisection strategy [beamsize] prec goal_prec timelimit randomseed");
		}

		load_sys(argv[1]);

	} catch(SyntaxError& e) {
		ibex_error(e.msg.c_str());
	}
        init(argc,argv);
}
  
  void OptimizerMinlpConfig::init(int arc, char** argv)
  {
  
	
//	if (simpl_level)
//		sys->set_simplification_level(simpl_level.Get());
        integerfilename = argv[2];

	ifstream fic(argv[2]);
	int nbint;
	fic >> nbint;
	int lint[nbint];
	for (int i=0;i<nbint; i++){
	  fic >> lint[i];
	  //	  cout << l[i] << " " ;
	}
	cout << endl;
        fic.close();

	
	filtering = argv[3];
	linearrelaxation = argv[4];
	bisection = argv[5];
	strategy = argv[6];
	

	int nbinput=6;

	beamsize = -1;
	if (strategy=="bs" || strategy== "beamsearch") { beamsize=atoi(argv[7]); nbinput++; }

	double prec= atof(argv[nbinput+1]);
	double goalprec= atof (argv[nbinput+2]);
	double timelimit= atof(argv[nbinput+3]);
	int trace=atoi(argv[nbinput+4]);
	int randomseed = atoi(argv[nbinput+5]);

	RNG::srand(randomseed);


	cout.precision(16);
	
	// the extended system
	double eqeps=1e-8;
	ext_sys = &rec(new ExtendedSystem(*sys,eqeps));
	norm_sys = &rec(new NormalizedSystem(*sys,eqeps));

       	// The integer variables 
      	BitSet b (ext_sys->nb_var);		

	//	cout << " nbint " << nbint << endl;
	for (int i=0; i< nbint;i++)
	    b.add(lint[i]);
	norm_sys->minlp=true;
	norm_sys->set_integer_variables(b);
	ext_sys->minlp=true;	
	ext_sys->set_integer_variables(b);

	
	/*
	cout << "file " << argv[1] << endl;
	cout << " filtering " << filtering;
        cout << " linearrelaxation " << linearrelaxation;
	cout << " bisection " << bisection ;
	cout << " strategy " << strategy ;
	cout << " randomseed " << randomseed << endl;
	 */

	set_eps_x(prec);
	set_rel_eps_f(goalprec);
	set_abs_eps_f(goalprec);

	// the optimizer : the same precision goalprec is used as relative and absolute precision
	//	Optimizer o(sys->nb_var,*ctcxn,*bs,loupfinder,*buffer,ext_sys.goal_var(),1.e-7,goalprec,goalprec);
	// the trace
	set_trace(trace);


	// the allowed time for search
	set_timeout(timelimit);
  }

unsigned int OptimizerMinlpConfig::nb_var() {
	return sys->nb_var;
}

void OptimizerMinlpConfig::load_sys(const char* filename) {
	std::size_t found = string(filename).find(".nl");

	if (found!=std::string::npos) {
		cerr << "\n\033[31mAMPL files can only be read with optimizerMinlp (ibex-opt-extra package).\n\n";
		exit(0);
	} else
		sys = &rec(new System(filename));
}

Linearizer* OptimizerMinlpConfig::get_linear_relax() {
	Linearizer* lr;

	if (linearrelaxation=="art")
		ibex_error("[OptimizerMinlpConfig]: ART mode available only in ibex-affine package\n");
	else if (linearrelaxation=="compo")
		ibex_error("[OptimizerMinlpConfig]: COMPO mode available only in ibex-affine package\n");
	else if (linearrelaxation=="xn")
		lr = &rec(new LinearizerXTaylor(*ext_sys));
	/*	else {
			stringstream ss;
			ss << "[optimizerMinlp] " << linearrelaxation << " is not an implemented relaxation mode ";
			ibex_error(ss.str().c_str());
		}
*/
	return lr;
}

Ctc& OptimizerMinlpConfig::get_ctc() {
  CtcInteger* integ= &rec(new CtcInteger(ext_sys->nb_var, *(ext_sys->get_integer_variables())));
	// the first contractor called
	CtcHC4* hc4 = &rec(new CtcHC4(ext_sys->ctrs,0.01,true));
	// hc4 inside acid and 3bcid : incremental propagation beginning with the shaved variable
	CtcHC4* hc44cid = &rec(new CtcHC4(ext_sys->ctrs,0.5,true));

	CtcCompo* hc44cidinteger = &rec(new CtcCompo (*integ,*hc44cid,*integ));
	CtcCompo* hc4integer = &rec (new CtcCompo (*integ,*hc4,*integ));
	// hc4 inside xnewton loop
	CtcHC4* hc44xn = &rec(new CtcHC4(ext_sys->ctrs,0.01,false));

	// The 3BCID contractor on all variables (component of the contractor when filtering == "3bcidhc4")
	Ctc3BCid* c3bcidhc4 = &rec(new Ctc3BCid(*hc44cid));
	// hc4 followed by 3bcidhc4 : the actual contractor used when filtering == "3bcidhc4"
	CtcCompo* hc43bcidhc4 = &rec(new CtcCompo(*integ, *hc4, *c3bcidhc4,*integ));

	// The ACID contractor (component of the contractor  when filtering == "acidhc4")
	CtcAcid* acidhc4 = &rec(new CtcAcid(*ext_sys, *hc44cidinteger, true));
	// hc4 followed by acidhc4 : the actual contractor used when filtering == "acidhc4"
	CtcCompo* hc4acidhc4 = &rec(new CtcCompo(*integ,*hc4, *acidhc4,*integ));

	Ctc* ctc;
	if (filtering == "hc4")
	  //		ctc = hc4;
	  ctc = hc4integer;
	else if (filtering =="acidhc4")
		ctc = hc4acidhc4;
	else if (filtering =="3bcidhc4")
		ctc = hc43bcidhc4;
	else {
		stringstream ss;
		ss << "[optimizerMinlp] " << filtering << " is not an implemented  contraction mode ";
		ibex_error(ss.str().c_str());
	}

	Linearizer* lr = get_linear_relax();
	Ctc* cxn;
	CtcPolytopeHull* cxn_poly;
	CtcCompo* cxn_compo;
	if (linearrelaxation=="compo" || linearrelaxation=="art"|| linearrelaxation=="xn") {
		cxn_poly = &rec(new CtcPolytopeHull(*lr));
		cxn_compo = &rec(new CtcCompo(*cxn_poly, *integ,*hc44xn,*integ));
	        if (sys->nb_ctr==0)
		  cxn = cxn_poly;
		else
		  cxn = &rec(new CtcFixPoint (*cxn_compo, relax_ratio));
	}
	//  the actual contractor  ctc + linear relaxation
	Ctc* ctcxn;
	if (linearrelaxation=="compo" || linearrelaxation=="art"|| linearrelaxation=="xn")
		ctcxn = &rec(new CtcCompo  (*ctc, *cxn));
	else
		ctcxn = ctc;
	Ctc* ctckkt;

	if (sys->nb_ctr == 0){
	  ctckkt = &rec(new CtcKuhnTucker(*norm_sys, true));
	  ctcxn = &rec(new CtcCompo  (*ctcxn, *ctckkt));
	  }

	return *ctcxn;
}

Bsc& OptimizerMinlpConfig::get_bsc() {
	Bsc* bs;

	double prec = get_eps_x();

	if (bisection=="roundrobin")
		bs = &rec(new RoundRobin (prec,0.45));
	else if (bisection== "largestfirst")
                bs = &rec(new OptimLargestFirst(ext_sys->goal_var(),true,prec,0.45));
	else if (bisection== "largestfirstnoobj")
                bs = &rec(new OptimLargestFirst(ext_sys->goal_var(),false,prec,0.45));
	else if (bisection=="smearsum")
		bs = &rec(new SmearSum(*ext_sys,prec,
				       rec(new OptimLargestFirst(ext_sys->goal_var(),true,prec,0.45))));
	else if (bisection=="smearmax")
		bs = &rec(new SmearMax(*ext_sys,prec,
				       rec(new OptimLargestFirst(ext_sys->goal_var(),true,prec,0.45))));
	else if (bisection=="smearsumrel")
                bs = &rec(new SmearSumRelative(*ext_sys,prec,
					       rec(new OptimLargestFirst(ext_sys->goal_var(),true,prec,0.45))));
	else if (bisection=="smearmaxrel")
		bs = &rec(new SmearMaxRelative(*ext_sys,prec,
					       rec(new OptimLargestFirst(ext_sys->goal_var(),true,prec,0.45))));
	else if  (bisection=="lsmear")
                bs = &rec (new LSmear(*ext_sys,prec,
				      rec(new OptimLargestFirst(ext_sys->goal_var(),true,prec,0.45)),
				      LSMEAR));
	else if (bisection=="lsmearmg")
	        bs = &rec (new LSmear(*ext_sys,prec,
				      rec(new OptimLargestFirst(ext_sys->goal_var(),true,prec,0.45))));
	else if (bisection=="minlpsmearsumrel") 
	  bs = &rec (new MinlpSmearSumRelative(*ext_sys,prec,
					       rec(new OptimLargestFirst(ext_sys->goal_var(),true,prec,0.45)),true));

        else if (bisection=="minlpsmearsumrelnoobj")
	  bs = &rec (new MinlpSmearSumRelative(*ext_sys,prec,
					       rec(new OptimLargestFirst(ext_sys->goal_var(),false,prec,0.45)),false));
	else if (bisection=="minlpsmearsum")
	  bs = &rec (new MinlpSmearSum(*ext_sys,prec,
			   rec(new OptimLargestFirst(ext_sys->goal_var(),true,prec,0.45)),true));
	
	  
	
        else if (bisection=="minlpsmearsumnoobj")
	  bs = &rec (new MinlpSmearSum(*ext_sys,prec,
			   rec(new OptimLargestFirst(ext_sys->goal_var(),true,prec,0.45))
					       ,false));
	else {
		stringstream ss;
		ss << "[optimizerMinlp] " << bisection << " is not an implemented  bisection mode ";
		ibex_error(ss.str().c_str());
	}

	return *bs;
}


LoupFinder& OptimizerMinlpConfig::get_loup_finder() {
	return rec(new LoupFinderDefault(*norm_sys, true));
	//LoupFinderDefault loupfinder (norm_sys,false);
}

CellBufferOptim& OptimizerMinlpConfig::get_cell_buffer() {
	CellBufferOptim* buffer;

	CellHeap* futurebuffer = &rec(new CellHeap(*ext_sys));
	CellHeap* currentbuffer = &rec(new CellHeap(*ext_sys));

	//cout << "strategy=[" << strategy << "]\n";
	if (strategy=="bfs")
		buffer = &rec(new CellHeap(*ext_sys));
	else if (strategy=="dh")
		buffer = &rec(new CellDoubleHeap(*ext_sys));
	else if (strategy=="bs")
		buffer = &rec(new CellBeamSearch(*currentbuffer, *futurebuffer, *ext_sys, beamsize));

	return *buffer;
}

int OptimizerMinlpConfig::goal_var() {
	return ext_sys->goal_var();
}

} // namespace ibex
