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

using namespace std;
using namespace ibex;

// write the ampl system for calling ipopt
void write_system_ampl(const System& sys){
	ofstream ficsys ("system.mod", ofstream::trunc);
	
	//	cout << *sys  <<endl;
	ficsys.precision(15);
	ficsys << sys  <<endl;
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


double correct_ipopt_sol (Vector&v , Optimizer& o, System& sys){
  IntervalVector box = sys.box;
  double eps=1.e-6;
  cout << " size " << v.size() << " boxsize " << sys.box.size() << endl;
  IntervalVector boxsol(v.size());
  for ( int i=0; i< v.size() ; i++){
    double epsi = eps;
    if (fabs(v[i])>1) epsi= eps*fabs(v[i]);
    boxsol[i]= box[i] & Interval(v[i]- epsi, v[i]+ epsi);
  }
  cout << "boxsol " << boxsol << endl;
  o.optimize(boxsol);
  o.report();
  o.set_uplo(NEG_INFINITY);
  return o.get_loup();
}

double check_ipopt(LoupFinder& loup_finder, Vector& v, bool integerobj, double integer_tolerance){
  double newub=POS_INFINITY;
    double loup=POS_INFINITY;
    if (loup_finder.integer_check(v)){
      newub=loup_finder.goal_ub(v);
      cout << "newub " << newub << endl;
      if ( loup_finder.is_inner(v)){
	cout << "is_inner " << endl;
	if (!integerobj
	    ||	    (std::ceil (newub) - newub) < integer_tolerance || (newub - std::floor(newub)) < integer_tolerance)
	  {	  
	  loup= newub;
	  cout << " loup= " << loup << endl;
	  }
      }
    }
    return loup;
}
 

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
     

bool ipopt_ampl_results(int n,string& status, double& obj, Vector & v){
    ifstream fic1 ("results.txt");
    if (fic1.good()){
      string a;
      int j;
      double otime=0;
      fic1 >> a;
      fic1 >> a;
      fic1 >> status; 
      cout << "status " << status << endl;
      fic1 >> a;
      fic1 >> a; 
      fic1 >> obj;
       cout << "obj " << obj << endl;
      
      fic1 >> a;
      fic1 >> a;
      fic1 >> otime;
      double solvertime =otime;

      fic1 >> a;
      fic1 >> a;
      fic1 >> otime;
      double ampltime =otime;

      cout << "ampltime " << ampltime << " ipopttime " << solvertime << endl;
      
      for (int i =0 ; i< n; i++){
	fic1 >> a;
	fic1 >> a;
	fic1 >> v[i];
	cout << "j" << a << "  v["<<i<<"] " <<  v[i] << endl;
      }
    
      fic1 >> a;
      fic1.close();
      return true;
    }
    return false;
  }
	  


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

	
	if (sys->minlp)	cout << " number of integer variables " << (sys->get_integer_variables())->size() << endl;
	if (sys->minlp)	cout << " integer variables " << *(sys->get_integer_variables()) << endl;
	// the extended system
	ExtendedSystem ext_sys(*sys,tolerance);
        NormalizedSystem norm_sys(*sys,tolerance);



	
	//	ext_sys.tolerance=tolerance;
	//	norm_sys.tolerance=tolerance;
	//	sys->tolerance=tolerance;
	
	/*	cout << "norm_sys" << endl;
	cout << norm_sys << endl;
	*/
	//	cout << "ext_sys" << ext_sys << endl;
	
	//	LoupFinderDefault loupfinder (norm_sys,true);
	LoupFinder* loupfinder;
	if (loupfindermethod=="ipoptxninhc4")
	  loupfinder = new LoupFinderDefaultIpopt (*sys,norm_sys,true);
	else if (loupfindermethod=="ipoptxn")
	  loupfinder = new LoupFinderDefaultIpopt (*sys,norm_sys,false);
	else if (loupfindermethod=="xninhc4")
	  loupfinder = new LoupFinderDefault (norm_sys,true);
	else if (loupfindermethod=="xn")
	  loupfinder = new LoupFinderDefault (norm_sys,false);
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
	if (o.trace) cout << " sys.box " << sys->box << endl;
	cout.precision(16);
	// ipopt preprocessing
        double obj=POS_INFINITY;
	Vector v(sys->nb_var);
	string status= "undefined";
	bool res;
	int res0;

	if (found!=std::string::npos){
	  string str = "time /libre/neveu/ibex-lib/ampl/ampl.linux-intel64/ipopt";
	  str.append(string(" "));
	  str.append(string(argv[1]));
	  str.append(" 'print_level 3 wantsol 2' > ipopt_res.txt");
	  //	  cout << str << endl;
	  res0= system(str.c_str());
          res=ipopt_direct_results(sys->nb_var,status,obj,v);
	}
	else{
	  write_system_ampl(*sys);
	  res0= system("/libre/neveu/ampl/ampl model_ipopt.run");
	  res=ipopt_ampl_results(sys->nb_var,status,obj,v);
	}
	
	cout << " resultats ipopt " << status;
	double initloup=POS_INFINITY;
	if (status=="solved"){
	  cout << "v " << v ;
	  cout << " obj " << obj;
	  initloup= check_ipopt(*loupfinder,v,integerobjective,goalprec);
	  if (initloup==POS_INFINITY)
	    initloup=correct_ipopt_sol(v,o,*sys);
	  cout << "loup_init " << initloup << endl;
	  o.set_loup(initloup);
	  IntervalVector lp(v);
	  if (initloup == POS_INFINITY) lp.set_empty();
	  o.set_loup_point(lp);
	}
	//integer objective
	o.integerobj=integerobjective;
	o.integer_tolerance=goalprec;
	// the allowed time for search
	o.timeout=timelimit;

	// the search itself 
	//	o.optimize(sys->box,initloup);
	//	std::ofstream Out("err.txt");
	//	std::streambuf* OldBuf = std::cerr.rdbuf(Out.rdbuf());
	if (o.trace) cout << " sys.box " << sys->box << endl;
	o.optimize(sys->box,initloup);
	//	std::cerr.rdbuf(OldBuf);

	// printing the results     
	o.report();
        cout << o.get_time() << "  " << o.get_nb_cells()+1 << endl;
        if (loupfindermethod == "ipoptxn" || loupfindermethod =="ipoptxninhc4")
	  cout << " ipopttime " << ((LoupFinderDefaultIpopt*)loupfinder)->finder_ipopt.ipopttime << " ampltime " << ((LoupFinderDefaultIpopt*)loupfinder)->finder_ipopt.ampltime << endl;
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
