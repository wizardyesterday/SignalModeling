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

[a1Covm,e1Covm] = covm(x,1);
g1Covm = atog(a1Covm);

[g1Fcov,e1Fcov] = fcov(x,1);
a1Fcov = gtoa(g1Fcov);

[g1Bcov,e1Bcov] = bcov(x,1);
a1Bcov = gtoa(g1Bcov);

[g1Burg,e1Burg] = burg(x,1);
a1Burg = gtoa(g1Burg);

[a1Mcov,e1Mcov] = mcov(x,1);
g1Mcov = atog(a1Mcov);

printf("First Order Reflection Coefficients, x1(n)\n");
printf("Covm Fcov Bcov Burg Mcov\n");
disp([g1Covm g1Fcov g1Bcov g1Burg g1Mcov]);
printf("\n");

printf("First Order Errors, x1(n)\n");
printf("Covm Fcov Bcov Burg Mcov\n");
disp([e1Covm e1Fcov e1Bcov e1Burg e1Mcov]);
printf("\n");

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (b), generate second-order models.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
[a2Covm,e2Covm] = covm(x,2);
g2Covm = atog(a2Covm);

[g2Fcov,e2Fcov] = fcov(x,2);
a2Fcov = gtoa(g2Fcov);

[g2Bcov,e2Bcov] = bcov(x,2);
a2Bcov = gtoa(g2Bcov);

[g2Burg,e2Burg] = burg(x,2);
a2Burg = gtoa(g2Burg);

[a2Mcov,e2Mcov] = mcov(x,2);
//g2Mcov = atog(a2Mcov);

printf("Second Order Reflection Coefficients, x1(n)\n");
printf("Covm Fcov Bcov Burg Mcov\n");
disp([g2Covm g2Fcov g2Bcov g2Burg]);
printf("\n");

printf("Second Order Errors, x1(n)\n");
printf("Covm Fcov Bcov Burg Mcov\n");
disp([e2Covm e2Fcov(2) e2Bcov(2) e2Burg(2) e2Mcov]);
printf("\n");

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (d), determine Burg error as order
// increases beyond 2.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
for p = 1:10
  [g,e] = burg(x,p);
//  printf("Order: %d  Error: %f\n",p,e(p));
end

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (e), repeat Part (a) and Part (b)
// for different input sequence.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Generate noise sequence.
noisegen(1,60,0.1);
v = feval(1:60,Noise);

// Generate time vector.
n = 1:60;
x = n .* (0.9).^n;

// Generate corrupted sequence.
y = x + v;

//----------------------------------------------
// Generate noise-free models.
//----------------------------------------------
[a1eCovm,e1eCovm] = covm(x,1);
g1eCovm = atog(a1eCovm);

[g1eFcov,e1eFcov] = fcov(x,1);
a1eFcov = gtoa(g1eFcov);

[g1eBcov,e1eBcov] = bcov(x,1);
a1eBcov = gtoa(g1eBcov);

[g1eBurg,e1eBurg] = burg(x,1);
a1eBurg = gtoa(g1eBurg);

[a1eMcov,e1eMcov] = mcov(x,1);
g1eMcov = atog(a1eMcov);

[a2eCovm,e2eCovm] = covm(x,2);
g2eCovm = atog(a2eCovm);

[g2eFcov,e2eFcov] = fcov(x,2);
a2eFcov = gtoa(g2eFcov);

[g2eBcov,e2eBcov] = bcov(x,2);
a2eBcov = gtoa(g2eBcov);

[g2eBurg,e2eBurg] = burg(x,2);
a2eBurg = gtoa(g2eBurg);

[a2eMcov,e2eMcov] = mcov(x,2);
g2eMcov = atog(a2eBurg);

printf("First Order Reflection Coefficients, x2(n)\n");
printf("Covm Fcov Bcov Burg Mcov\n");
disp([g1eCovm g1eFcov g1eBcov g1eBurg g1eMcov]);
printf("\n");

printf("First Order Errors, x2(n)\n");
printf("Covm Fcov Bcov Burg Mcov\n");
disp([e1eCovm e1eFcov e1eBcov e1eBurg e1eMcov]);
printf("\n");

printf("Second Order Reflection Coefficients, x2(n)\n");
printf("Covm Fcov Bcov Burg Mcov\n");
disp([g2eCovm g2eFcov g2eBcov g2eBurg g2eMcov]);
printf("\n");

printf("Second Order Errors, x2(n)\n");
printf("Covm Fcov Bcov Burg Mcov\n");
disp([e2eCovm e2eFcov(2) e2eBcov(2) e2eBurg(2) e2eMcov]);
printf("\n");

//----------------------------------------------
// Generate noise-perturbed models.
//----------------------------------------------
[a1yCovm,e1yCovm] = covm(y,1);
g1yCovm = atog(a1yCovm);

[g1yFcov,e1yFcov] = fcov(y,1);
a1yFcov = gtoa(g1yFcov);

[g1yBcov,e1yBcov] = bcov(y,1);
a1yBcov = gtoa(g1yBcov);

[g1yBurg,e1yBurg] = burg(y,1);
a1yBurg = gtoa(g1yBurg);

[a1yMcov,e1yMcov] = mcov(y,1);
g1yMcov = atog(a1yMcov);

[a2yCovm,e2yCovm] = covm(y,2);
g2yCovm = atog(a2yCovm);

[g2yFcov,e2yFcov] = fcov(y,2);
a2yFcov = gtoa(g2yFcov);

[g2yBcov,e2yBcov] = bcov(y,2);
a2yBcov = gtoa(g2yBcov);

[g2yBurg,e2yBurg] = burg(y,2);
a2yBurg = gtoa(g2yBurg);

[a2yMcov,e2yMcov] = mcov(y,2);
g2yMcov = atog(a2yMcov);

printf("First Order Reflection Coefficients, y(n)\n");
printf("Covm Fcov Bcov Burg Mcov\n");
disp([g1yCovm g1yFcov g1yBcov g1yBurg g1yMcov]);
printf("\n");

printf("First Order Errors, y(n)\n");
printf("Covm Fcov Bcov Burg Mcov\n");
disp([e1yCovm e1yFcov e1yBcov e1yBurg e1yMcov]);
printf("\n");

printf("Second Order Reflection Coefficients, y(n)\n");
printf("Covm Fcov Bcov Burg Mcov\n");
disp([g2yCovm g2yFcov g2yBcov g2yBurg g2yMcov]);
printf("\n");

printf("Second Order Errors, y(n)\n");
printf("Covm Fcov Bcov Burg Mcov\n");
disp([e2yCovm e2yFcov(2) e2yBcov(2) e2yBurg(2) e2yMcov]);
printf("\n");
