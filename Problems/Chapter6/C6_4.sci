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
x1 = 1:60;

// Generate a unit sample sequence.
delta = zeros(1,60);
delta(1) = 1;

[a1Covm,e1Covm] = covm(x1,1);
g1Covm = atog(a1Covm);
xhat1Covm = filterBlock(delta,1,a1Covm(2:$));

[g1Fcov,e1Fcov] = fcov(x1,1);
a1Fcov = gtoa(g1Fcov);
xhat1Fcov = filterBlock(delta,1,a1Fcov(2:$));

[g1Bcov,e1Bcov] = bcov(x1,1);
a1Bcov = gtoa(g1Bcov);
xhat1Bcov = filterBlock(delta,1,a1Bcov(2:$));

[g1Burg,e1Burg] = burg(x1,1);
a1Burg = gtoa(g1Burg);
xhat1Burg = filterBlock(delta,1,a1Burg(2:$));

[a1Mcov,e1Mcov] = mcov(x1,1);
g1Mcov = atog(a1Mcov);
xhat1Mcov = filterBlock(delta,1,a1Mcov(2:$));

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
[a2Covm,e2Covm] = covm(x1,2);
g2Covm = atog(a2Covm);
xhat2Covm = filterBlock(delta,1,a2Covm(2:$));

[g2Fcov,e2Fcov] = fcov(x1,2);
a2Fcov = gtoa(g2Fcov);
xhat2Fcov = filterBlock(delta,1,a2Fcov(2:$));

[g2Bcov,e2Bcov] = bcov(x1,2);
a2Bcov = gtoa(g2Bcov);
xhat2Bcov = filterBlock(delta,1,a2Bcov(2:$));

[g2Burg,e2Burg] = burg(x1,2);
a2Burg = gtoa(g2Burg);
xhat2Burg = filterBlock(delta,1,a2Burg(2:$));

[a2Mcov,e2Mcov] = mcov(x1,2);
//g2Mcov = atog(a2Mcov);
xhat2Mcov = filterBlock(delta,1,a2Mcov(2:$));

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
  [g,e] = burg(x1,p);
//  printf("Order: %d  Error: %f\n",p,e(p));
end

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (e), repeat parts (a) and (b) for
// different input sequence.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Generate noise sequence.
noisegen(1,60,0.1);
v = feval(1:60,Noise);

// Generate time vector.
n = 1:60;
x2 = n .* (0.9).^n;

// Generate corrupted sequence.
y = x2 + v;

//----------------------------------------------
// Generate noise-free models.
//----------------------------------------------
[a1eCovm,e1eCovm] = covm(x2,1);
g1eCovm = atog(a1eCovm);

[g1eFcov,e1eFcov] = fcov(x2,1);
a1eFcov = gtoa(g1eFcov);

[g1eBcov,e1eBcov] = bcov(x2,1);
a1eBcov = gtoa(g1eBcov);

[g1eBurg,e1eBurg] = burg(x2,1);
a1eBurg = gtoa(g1eBurg);

[a1eMcov,e1eMcov] = mcov(x2,1);
g1eMcov = atog(a1eMcov);

[a2eCovm,e2eCovm] = covm(x2,2);
g2eCovm = atog(a2eCovm);

[g2eFcov,e2eFcov] = fcov(x2,2);
a2eFcov = gtoa(g2eFcov);

[g2eBcov,e2eBcov] = bcov(x2,2);
a2eBcov = gtoa(g2eBcov);

[g2eBurg,e2eBurg] = burg(x2,2);
a2eBurg = gtoa(g2eBurg);

[a2eMcov,e2eMcov] = mcov(x2,2);
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

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Plot results from parts (a) and (b).
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
subplot(5,2,1);
title('xhat(n), First-Order Covm');
plot(xhat1Covm);
xgrid();

subplot(5,2,3);
title('xhat(n), First-Order Fcov');
plot(xhat1Fcov);
xgrid();

subplot(5,2,5);
title('xhat(n), First-Order Bcov');
plot(xhat1Bcov);
xgrid();

subplot(5,2,7);
title('xhat(n), First-Order Burg');
plot(xhat1Burg);
xgrid();

subplot(5,2,9);
title('xhat(n), First-Order Mcov');
plot(xhat1Mcov);
xgrid();

subplot(5,2,2);
title('xhat(n), Second-Order Covm');
plot(xhat2Covm);
xgrid();

subplot(5,2,4);
title('xhat(n), Second-Order Fcov');
plot(xhat2Fcov);
xgrid();

subplot(5,2,6);
title('xhat(n), Second-Order Bcov');
plot(xhat2Bcov);
xgrid();

subplot(5,2,8);
title('xhat(n), Second-Order Burg');
plot(xhat2Burg);
xgrid();

subplot(5,2,10);
title('xhat(n), Second-Order Mcov');
plot(xhat2Mcov);
xgrid();

