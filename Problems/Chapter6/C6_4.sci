//**********************************************************************
// File Name: C6_4.sci
//**********************************************************************

exec('utils.sci',-1);

//**********************************************************************
// Mainline code.
//**********************************************************************

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (a), generate first-order models.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Generate input sequence, x(n) = (n+1)u(n).
x = 1:60;

[g,e1Covm] = covm(x,1);
a1Covm = gtoa(g);

[g,e1Fcov] = fcov(x,1);
a1Fcov = gtoa(g);

[g,e1Bcov] = bcov(x,1);
a1Bcov = gtoa(g);

[g,e1Burg] = burg(x,1);
a1Burg = gtoa(g);

[g,e1Mcov] = mcov(x,1);
a1Mcov = gtoa(g);

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (b), generate second-order models.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
[g,e2Covm] = covm(x,2);
a2Covm = gtoa(g);

[g,e2Fcov] = fcov(x,2);
a2Fcov = gtoa(g);

[g,e2Bcov] = bcov(x,2);
a2Bcov = gtoa(g);

[g,e2Burg] = burg(x,2);
a2Burg = gtoa(g);

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (d), determine Burg error as order
// increases beyond 2.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
for p = 1:10
  [g,e] = burg(x,p);
  printf("Order: %d  Error: %f\n",p,e(p));
end

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (e), repeat Part (a) and Part (b)
// for different input sequence.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/



[g,e2Mcov] = mcov(x,2);
a2Mcov = gtoa(g);


