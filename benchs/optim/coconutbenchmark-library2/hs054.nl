g3 0 1 0	# problem hs054
 6 1 1 0 1	# vars, constraints, objectives, ranges, eqns
 0 1	# nonlinear constraints, objectives
 0 0	# network constraints: nonlinear, linear
 0 6 0	# nonlinear vars in constraints, objectives, both
 0 0 0 1	# linear network variables; functions; arith, flags
 0 0 0 0 0	# discrete variables: binary, integer, nonlinear (b,c,o)
 2 6	# nonzeros in Jacobian, gradients
 0 0	# max name lengths: constraints, variables
 0 0 0 0 0	# common exprs: b,c,o,c1,o1
C0	#constr1
n0
O0 0	#obj
o54	#sumlist
7
o2	#*
n1.0416666666666667
o5	#^
v0	#y[1]
n2
o2	#*
n1.0416666666666667
o2	#*
o2	#*
n0.4
v0	#y[1]
v1	#y[2]
o2	#*
n1.0416666666666667
o5	#^
v1	#y[2]
n2
o5	#^
v2	#y[3]
n2
o5	#^
v3	#y[4]
n2
o5	#^
v4	#y[5]
n2
o5	#^
v5	#y[6]
n2
x6	# initial guess
0 -0.5
1 0.5
2 0.2857142857142857
3 -0.16
4 0.04
5 -0.1
r	#1 ranges (rhs's)
4 0.45
b	#6 bounds (on variables)
0 -1.25 1.25
0 -11 9
0 -0.2857142857142857 1.1428571428571428
0 -0.2 0.2
0 -20.019999999999996 19.98
0 -0.2 0.2
k5	#intermediate Jacobian column lengths
1
2
2
2
2
J0 2
0 1
1 0.5
G0 6
0 0
1 0
2 0
3 0
4 0
5 0
