//============================================================================
//                                  I B E X                                   
// File        : valdationvp.cpp
// Author      : Bertrand Neveu
// Copyright   : Ecole des Mines de Nantes (France)
// License     : See the LICENSE file
// Created     : Dec 12, 2022
// Last Update : Dec  12, 2022
//============================================================================

/*
 validation of eigenvalues using characteristic polynomial
input 3 arguments
file1 : vp.ampl
file2 : poly.ampl
int : n: number of eigenvalues in vp.ampl

output: 
the smallest eigenvalue is written in file vp_cert.ampl"
 */

#include "ibex.h"
using namespace std;
using namespace ibex;
int main(int argc, char** argv){

  cout << argv[1] << endl;
  ifstream fic (argv[1]);
  int n = atoi(argv[3]);
  string a;
  int  b;
  double vp[n];
  
  fic >> a;

  fic >> a;

  fic >> a;

  fic >> a;

  for (int i=0; i< n; i++){
    fic >> b;
    fic >> vp[i];
  }
  fic >> a;
  fic >> a;
  fic.close();
  for (int i=0; i< n; i++){
    //    cout << vp[i] << endl;
  }
  int k=2;
  double epsilon=0;
  Interval valprop[n];
  for (int i=0; i<n; i++)
    valprop[i]=Interval(vp[i]-epsilon,vp[i]+epsilon);
  for (int i=0; i< n; i++){
    //    cout << valprop[i] << endl;
  }
  ifstream fic2 (argv[2]);
  double poly[n+1];
  for (int i=0; i< n+1; i++){
    fic2 >> poly[i];
  }
  fic2.close();
  Interval polycar[n+1];
  for (int i=0; i< n+1; i++){
    polycar[i]=Interval(poly[i],poly[i]);
  }
  for (int i=0; i< n+1; i++){
    //    cout << polycar[i] << endl;
  }	    
  double minvalprop[n];
  Interval result[n];
  for (int i=0;i<n; i++)
    result[i]=Interval(0,0);
  for (int j=0 ; j< n; j++){
    for (int i=0; i< n; i++){
      result[j]+= polycar[i] * pow (valprop[j], n-i);
    }
    result[j] +=polycar[n];
    //    cout << j << " epsilon 0 " << result[j] << endl;
    if (! result[j].contains(0)){
      epsilon=1.e-20;
      while (! result[j].contains(0) && epsilon < 1000){
	epsilon = k*epsilon;
	result[j]=Interval(0,0);
	for (int i=0; i< n; i++){
	  result[j]+= polycar[i] * pow (Interval(valprop[j].lb()-epsilon,valprop[j].ub()+epsilon), n-i);
	}
	result[j] +=polycar[n];
      }
    }
    minvalprop[j]=valprop[j].lb()-epsilon ;
    
    //    cout << j << " "  << epsilon << " "  << result[j] << "  " << minvalprop[j] << endl;
      
  }
  double minvalprop0= minvalprop[0];
  int index=0;
  for (int j=1; j<n; j++)
    if (minvalprop[j]<  minvalprop0){
      minvalprop0=minvalprop[j]; index=j;
    }
  ofstream fic1 ("vp_cert.ampl", ofstream::trunc);
  fic1 << "data;" << endl;
  fic1 << "param vp_cert:= " << minvalprop0 << ";" << endl;;
  fic1 << "end;" << endl;
  fic1.close();
    
  cout.precision(12);
  
  cout << " index " << index << " minvalprop " << minvalprop0 << " epsilon " <<  valprop[index].lb() - minvalprop0 <<  endl;
    
}


      
  
  
  
