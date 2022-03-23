//**********************************************************************
// File Name: C7_1.sci
//**********************************************************************

// Bring in all of the utility functions.
exec('utils.sci',-1);

//**********************************************************************
// Mainline code.
//**********************************************************************

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (b): Evaluate the error for filters of p = 1,2,...20, and
// compare these errors to those obtained from the causial
// Wiener filter.
// Note from Example 7.3.2 of the text that the mean-square error
// for the causal Wiener filter is 0.3.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
alpha = 0.8

// Generate model order sequence.
p = 1:20;

// Generate lag sequence.
k = 1:100;

// Compute autocorrelation sequence.
rd = alpha ^ k;

// Allocate space for mean-square errors.
eb = zeros(1,20);

// Generate error sequence.
for i = 1:20
  [w,eb(p(i))] = FirWienerFilter(rd,1,p(i));
end

// Display the results.
printf("p, eb(p)\n");
disp([p' eb']);

// Plot the results.
subplot(211);
title("Mean-square error Versus p (rd(k) = alpha^k) for Fir Wiener Filter");

plot(p,eb);

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (c): With p = 10, plot the mean-square error versus alpha
// for alpha = 0.1,0.2,...,0.9. 
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Generate parameters.
alpha = 0.1:0.1:0.9;

// Allocate space for mean-square errors.
ec = zeros(1,9);

// We want to evaluate the errors for tenth-order Wiener filters.
p = 10;

for i = 1:9
  rd = alpha(i)^k;
  [w,ec(i)] = FirWienerFilter(rd,1,p);
end

// Plot the results.
subplot(212);
title("Mean-square error Versus Alpha (rd(k) = alpha^k) for Fir Wiener Filter");
plot(alpha,ec);





