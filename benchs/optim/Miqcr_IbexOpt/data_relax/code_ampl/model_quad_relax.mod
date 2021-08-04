param n;
param n_y;
param m;
param p;
param cons;
param lb {i in 1..n-n_y}default 0; 
param ub {i in 1..n-n_y}default 0; 
param Q {i in 1..n,j in 1..n} default 0;
param c {i in 1..n} default 0;
param A {r in 1..m,i in 1..n} default 0;
param b {r in 1..m} default 0;
param D {s in 1..p,i in 1..n} default 0;
param e {s in 1..p} default 0;
param assoc_prod{s in 1..n_y, t in 1..2} default 0;
var x {i in 1..n};


minimize obj :
sum {i in 1..n,j in 1..n}  Q[i,j] * x[i] * x[j]
+ sum {i in 1..n}c[i] * x[i]
+cons;


subject to eq_lin{s in 1..m}: sum {i in 1..n} A[s,i]*x[i] = b[s];
subject to ineq_lin{s in 1..p}: sum {i in 1..n} D[s,i]*x[i] <= e[s];

subject to cor1{s in n-n_y+1..n}:
 x[s] <=  ub[assoc_prod[s-n+n_y,2]]*x[assoc_prod[s-n+n_y,1]] + lb[assoc_prod[s-n+n_y,1]] * x[assoc_prod[s-n+n_y,2]] - ub[assoc_prod[s-n+n_y,2]]*lb[assoc_prod[s-n+n_y,1]]; 
subject to cor2{s in n-n_y+1..n : assoc_prod[s-n+n_y,1] <> assoc_prod[s-n+n_y,2]}:
x[s] <=  ub[assoc_prod[s-n+n_y,1]]*x[assoc_prod[s-n+n_y,2]] + lb[assoc_prod[s-n+n_y,2]] * x[assoc_prod[s-n+n_y,1]] - ub[assoc_prod[s-n+n_y,1]]*lb[assoc_prod[s-n+n_y,2]]; 
subject to cor3{s in n-n_y+1..n}:
 x[s] >=  ub[assoc_prod[s-n+n_y,1]]*x[assoc_prod[s-n+n_y,2]] + ub[assoc_prod[s-n+n_y,2]] * x[assoc_prod[s-n+n_y,1]] - ub[assoc_prod[s-n+n_y,1]]*ub[assoc_prod[s-n+n_y,2]]; 
subject to cor4{s in n-n_y+1..n}:
 x[s] >=  lb[assoc_prod[s-n+n_y,1]]*x[assoc_prod[s-n+n_y,2]] + lb[assoc_prod[s-n+n_y,2]] * x[assoc_prod[s-n+n_y,1]] - lb[assoc_prod[s-n+n_y,1]]*lb[assoc_prod[s-n+n_y,2]]; 



subject to cor5{s in n-n_y+1..n : assoc_prod[s-n+n_y,1] == assoc_prod[s-n+n_y,2]}:
  x[s] >=  x[assoc_prod[s-n+n_y,2]];



subject to upper_bound{i in 1..n-n_y}: x[i] <= ub[i]; 
subject to lower_bound{i in 1..n-n_y}: x[i] >=lb[i]; 


