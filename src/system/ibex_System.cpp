//============================================================================
//                                  I B E X                                   
// File        : ibex_System.cpp
// Author      : Gilles Chabert
// Copyright   : IMT Atlantique (France)
// License     : See the LICENSE file
// Created     : Jun 12, 2012
// Last Update : Apr 08, 2019
//============================================================================

#include "ibex_System.h"
#include "ibex_SystemFactory.h"
#include "ibex_SyntaxError.h"
#include "ibex_UnknownFileException.h"
#include "ibex_ExprCopy.h"
#include "ibex_Expr2Minibex.h"
#include "ibex_P_Struct.h"
#include "ibex_Domain.h"
#include "ibex_Exception.h"
#include "ibex_String.h"
#include "ibex_BitSet.h"

#include <stdio.h>
#include <sstream>
#include <fstream>


#ifndef _WIN32 // MinGW does not support mutex
#include <mutex>
namespace {
std::mutex mtx;
}
#define LOCK mtx.lock()
#define UNLOCK mtx.unlock()
#else
#define LOCK
#define UNLOCK
#endif

extern int ibexparse();
extern void ibexparse_string(const char* syntax);
//extern int ibex_lineno;
extern void ibexrestart(FILE *);
extern FILE* ibexin;

using namespace std;

namespace ibex {

namespace {

class SystemCopy : public SystemFactory {

public:

	SystemCopy(const System& sys, const System::copy_mode& mode) {

		// do not initialize variables with sys.f_ctrs.args
		// since f may be uninitialized (unconstrained problem)
		add_var(sys.args,sys.box);

		if (mode==System::COPY && sys.goal!=NULL)
			add_goal(*sys.goal);

		for (int i=0; i<sys.nb_ctr; i++) {
			if (mode==System::COPY ||
					(sys.ctrs[i].op==EQ && mode==System::EQ_ONLY) ||
					(sys.ctrs[i].op!=EQ && mode==System::INEQ_ONLY))

				// TODO: we lose DAG structure here...

				add_ctr(sys.ctrs[i]);
		}
	}
};

} // end anonymous namespace

  System::System() : id(next_id()), nb_var(0), nb_ctr(0), goal(NULL), ops(NULL), box(1), integer_variables(NULL) /* tmp */ {

}

System::System(const char* filename, int simpl_level) : id(next_id()), nb_var(0), nb_ctr(0), goal(NULL), ops(NULL), box(1) /* tmp */ {
	FILE *fd;
	if ((fd = fopen(filename, "r")) == NULL) throw UnknownFileException(filename);
	load(fd, simpl_level);

}

System::System(int n, const char* syntax, int simpl_level) : id(next_id()), nb_var(n), /* NOT TMP (required by parser) */
		                                    nb_ctr(0), goal(NULL), ops(NULL), box(1) /* tmp */ {
	LOCK;
	try {
		parser::pstruct = new parser::P_StructChoco(*this, simpl_level);
		ibexparse_string(syntax);
		delete parser::pstruct;
		parser::pstruct = NULL;
	} catch(SyntaxError& e) {
		delete parser::pstruct;
		parser::pstruct = NULL;
		UNLOCK;
		throw e;
	}
	UNLOCK;
	//TODO / find integer variables
	integer_variables= new BitSet(nb_var);
}

System::System(const System& sys, copy_mode mode) : id(next_id()), nb_var(0), nb_ctr(0), goal(NULL), ops(NULL), box(1) {

	switch(mode) {
	case COPY :      init(SystemCopy(sys,COPY)); break;
	case INEQ_ONLY:  init(SystemCopy(sys,INEQ_ONLY)); break;
	case EQ_ONLY:    init(SystemCopy(sys,EQ_ONLY)); break;
	}

}

std::string System::minibex(bool human) const {
	stringstream s;
	s << "variables\n";

	Array<Domain> domains(args.size());
	for (int i=0; i<args.size(); i++) {
		domains.set_ref(i,*new Domain(args[i].dim));
	}

	
	ibex::load(domains, box);

	int index=0;
	
	for (int i=0; i<args.size(); i++) {
		const ExprSymbol& x=args[i];
	        s << x.name;
	
		if (x.dim.nb_rows()>1 || x.dim.nb_cols()>1) s << '[' << x.dim.nb_rows() << ']';
		if (x.dim.nb_cols()>1) s << '[' << x.dim.nb_cols() << ']';
		if ((*integer_variables)[index])
		  s << " integer ";
		s << " in ";
		ExprPrinter().print(s,domains[i],human);
		s << ";\n";
		index+=x.dim.size();
	}

	s << '\n';

	if (goal) {
		s << goal->minibex(human) << "\n\n";
	}

	if (nb_ctr>0) {
		s << f_ctrs.minibex(human) << "\n\n";
	}

	if (goal) {
		s << "minimize\n " << goal->name << "(";
		for (int i=0; i<args.size(); i++) {
			if (i>0) s << ',';
			s << args[i];
		}
		s << ");\n\n";
	}

	if (nb_ctr>0) {
		s << "constraints\n";
		char* main_node = next_generated_var_name();
		s <<  main_node << "=" << f_ctrs.name << '(';
		for (int i=0; i<args.size(); i++) {
			if (i>0) s << ',';
			s << args[i];
		}
		s << ");\n";
		for (int i=0; i<f_ctrs.image_dim(); i++) {
			 s << main_node << '[' << i << ']'<< ops[i] << "0;\n";
		}

		free(main_node);
	}
	s << "end";
	s.flush();

	for (int i=0; i<args.size(); i++) {
		delete &domains[i];
	}

	return s.str();
}

std::ostream& operator<<(std::ostream& os, const System& sys) {

	os << "variables: " << endl << "  ";
	int index=0;
	for (int i=0; i<sys.args.size(); i++) {
		const ExprSymbol& x = sys.args[i];
		os << x;
		if((*(sys.get_integer_variables()))[index])
		  os << " integer" ;
		   
		if (x.dim.nb_rows()>1) os << '[' << x.dim.nb_rows() << ']';
		if (x.dim.nb_cols()>1) {
			if (x.dim.nb_rows()==1) os << "[1]";
			os << '[' << x.dim.nb_cols() << ']';
		}
		if (i<sys.args.size()-1) os << ", ";
		index+=x.dim.size();
	}
	os << endl;

	os << "box: " << endl << "  ";
	os << sys.box << endl;


	os << "goal: " << endl;
	if (sys.goal!=NULL)
	    os << "  " << sys.goal->expr() << endl;
	else
		os << "  (none)" << endl;
	if (sys.nb_ctr>0) {
		os << "constraints:" << endl;
		for (int i=0; i<sys.ctrs.size(); i++)
			os << "  " << sys.ctrs[i] << endl;
	}
	return os;
}

Domain& System::constant(const std::string& name) {
	const char* _name=name.c_str();
	try {
		Domain* domain = mutable_constants[_name];
		return *domain;
	} catch (SymbolMap<Domain*>::SymbolNotFound&) {
		ibex_error("No mutable constant named \""+name+"\")");
	}
}

void System::load(FILE* fd, int simpl_level) {

	LOCK;

	ibexin = fd;

	try {
		parser::pstruct = new parser::P_StructSystem(*this, simpl_level);
		ibexparse();
		delete parser::pstruct;
		parser::pstruct = NULL;
	}

	catch(SyntaxError& e) {
		delete parser::pstruct;
		parser::pstruct = NULL;
		fclose(fd);
		ibexrestart(ibexin);
		UNLOCK;
		throw e;
	}

	fclose(fd);

	UNLOCK;
}

System::~System() {

	if (goal) delete goal;

	for (int i=0; i<ctrs.size(); i++)
		delete &ctrs[i];

	if (nb_ctr==0) {
		// delete the symbols of the global function "f" that
		// has never been created
		for (int i=0; i<args.size(); i++) delete &args[i];
	}

	if (ops) delete[] ops;

	for (IBEXMAP(Domain*)::iterator it=mutable_constants.begin(); it!=mutable_constants.end(); ++it) {
		delete it->second;
	}
	if (integer_variables) delete integer_variables;
}



/*
 * Note:
 *  A constraint g(x)<=0 is considered as potentially active
 *  if [g]([x])<=0, to conform with the definition of activity
 *  in the context of optimization.
 *  Otherwise (that is, if the constraint was considered as inactive
 *  in this case), the KKT contractor would be incorrect (losing
 *  solutions).
 *  An alternative would be to create two different concepts:
 *  activity:      g(x)<=0 is inactive if [g]([x])<0
 *  effectiveness: g(x)<=0 is ineffective if [g]([x])<=0 (for
 *                 contractors like HC4).
 *  However, the case where a constraint is not effective but
 *  active almost never happens.
 */
   bool  System::is_inactive( Interval& gx, CmpOp op) const {
	bool inactive;
	switch (op) {
	case LT:  inactive=gx.ub()<0; break;
	case LEQ: inactive=gx.ub()<0; break;
	case EQ:  inactive=(gx==Interval::zero()); break;
	case GEQ: inactive=gx.lb()>0; break;
	case GT:  inactive=gx.lb()>0; break;
	}
	return inactive;
}
  bool System::is_ineffective( Interval& gx, CmpOp op)const {
	bool ineffective;
	switch (op) {
	 
	case LT:  ineffective=gx.ub()<0; break;
	case LEQ: ineffective=gx.ub()<= tolerance; break;
	case EQ:  ineffective=(gx==Interval::zero()); break;
	case GEQ: ineffective=gx.lb()>=-tolerance; break;
	case GT:  ineffective=gx.lb()>0; break;
	}
	return ineffective;
}


//TODO: improvement: don't create a bitset in return when not necessary?

BitSet System::active_ctrs(const IntervalVector& box) const {

	if (nb_ctr==0) return BitSet::empty(1);

	BitSet active(BitSet::all(f_ctrs.image_dim()));

	IntervalVector res = f_ctrs.eval_vector(box);

	for (int c=0; c<f_ctrs.image_dim(); c++) {
		if (is_inactive(res[c],ops[c])) active.remove(c);
	}
	return active;
}

  BitSet System::effective_ctrs(const IntervalVector& box) const {

	if (nb_ctr==0) return BitSet::empty(1);

	BitSet effective(BitSet::all(f_ctrs.image_dim()));

	IntervalVector res = f_ctrs.eval_vector(box);
	//	cout << " res " << res << endl;
	for (int c=0; c<f_ctrs.image_dim(); c++) {
		if (is_ineffective(res[c],ops[c])) effective.remove(c);
	}
	return effective;
}


  
IntervalVector System::active_ctrs_eval(const IntervalVector& box) const {

	BitSet b=active_ctrs(box);

	assert(!b.empty());

	IntervalVector ev(b.size());
	int c;
	for (int i=0; i<b.size(); i++) {
		c=(i==0? b.min() : b.next(c));
		ev[i] = f_ctrs[c].eval(box); // not: a call to f_ctrs.eval_vector would benefit from the DAG...
	}
	return ev;
}

IntervalMatrix System::active_ctrs_jacobian(const IntervalVector& box) const {

	BitSet b=active_ctrs(box);

	assert(!b.empty());

	IntervalMatrix J(b.size(),nb_var);

	J=f_ctrs.jacobian(box,b);

	return J;
}

  bool System::is_inner(const IntervalVector& box) const {
  //	return active_ctrs(box).empty();
    //    std::cout << " effective_ctrs  " <<  effective_ctrs(box) << std::endl;
  	return effective_ctrs(box).empty();

  }

  bool System::is_inner(const Vector & pt) const {
    IntervalVector box(pt);
    return is_inner(box);
  }
  
  void System::set_integer_variables(const BitSet & integer_vars){
    if (integer_variables) delete integer_variables;
    integer_variables=new BitSet(integer_vars);
  }

  void System::set_integer_variables(const vector<int> & integer_vars){
    if (integer_variables) delete integer_variables;
    integer_variables= new BitSet (nb_var);
    if (integer_vars.size() >0){
      minlp=true;      
      for (int i=0; i<integer_vars.size();i++)
	integer_variables->add(integer_vars[i]);
    }
  }

  BitSet* System::get_integer_variables() const {
    return integer_variables;
  }

  bool System::is_integer(int i) const{
    if  (minlp)
      return (*integer_variables)[i];
    else
      return false;
  }
  

  vector<string> System::var_names() const {

	vector<string> var_names;
	int v=0;
	for (int s=0; s<args.size(); s++) {
		const ExprSymbol& x=args[s];
		switch (x.dim.type()) {
		case Dim::SCALAR:
			var_names.push_back(x.name);
			break;
		case Dim::ROW_VECTOR:
		case Dim::COL_VECTOR:
			for (int i=0; i<x.dim.vec_size(); i++)
				var_names.push_back(string(x.name)+'('+to_string(i+1)+')');
			break;
		default: // MATRIX
			for (int i=0; i<x.dim.nb_rows(); i++)
				for (int j=0; j<x.dim.nb_cols(); j++)
					var_names.push_back(string(x.name)+'('+to_string(i+1)+','+to_string(j+1)+')');
			break;
		}
		v+=x.dim.size();
	}
	assert(v==nb_var);
	return var_names;
  }

  

  vector<int> System::find_integer_variables (char* filename) {
    vector<int> integers;
    std::size_t found = string(filename).find(".nl");
    if (found!=std::string::npos){
      ifstream fic (filename);
      string a;
      int nbvar;
      for (int i=0; i<7 ;i++){
	fic >> a;
    }

      fic >> nbvar;
      //  cout << " nbvar " << nbvar << endl;
      for (int i=0; i<10 ;i++){
	fic >> a;
	//    cout << "  " << a ;
      }
      for (int i=0; i<6 ;i++){
	fic >> a;
	//    cout << "  " << a ;
      }
      for (int i=0; i<7 ;i++){
	fic >> a;
	//    cout << "  " << a ;
      }
      int nlvarc;
      int nlvaro;
      int nlvarb;
      fic >> nlvarc;
      fic >> nlvaro;
      fic >> nlvarb;
      //  cout << " nlvarc " << nlvarc << " nlvaro " <<  nlvaro << " nlvarb " << nlvarb << endl;
      for (int i=0; i<7 ;i++){
	fic >> a;
	//    cout << "  " << a ;
      }
      for (int i=0; i<11 ;i++){
	fic >> a;
	//    cout << "  " << a ;
      }
      int lbin ;
      int lint ;
      int nlbd;
      int nlcd;
      int nlod;
      fic >> lbin;
      fic >> lint;
      fic >> nlbd;
      fic >> nlcd;
      fic >> nlod;
      //  cout << "lbin " << lbin << " lint " << lint << " nlbd " << nlbd << " nlcd " << nlcd << " nlod " << nlod << endl;
      // the integer variables are found among the variables by using the variable ordering of the nl files described in
      //Hooking your solver to ampl by David M Gay
  
  
      fic.close ();

      for (int i=nlvarb-nlbd;i<nlvarb;i++)
	{//cout << i << " " << endl;
	  integers.push_back(i);
	}
      for (int i=nlvarc-nlcd;i< nlvarc; i++)
	{//cout << i << " " << endl;
	  integers.push_back(i);
	}
      for (int i=nlvaro-nlod;i< nlvaro; i++)
	{//cout << i << " " << endl;
	  integers.push_back(i);
	}
      for (int  i= nbvar-lbin-lint; i< nbvar; i++)
	{//cout << i << " " << endl;
	  integers.push_back(i);
	}
    }  
  return integers;
  
  }

          

} // end namespace ibex
