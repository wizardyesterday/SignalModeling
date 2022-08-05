//**********************************************************************
// File Name: C9_5.sci
//**********************************************************************

// Bring in all of the utility functions.
exec('utils.sci',-1);

//**********************************************************************
// Mainline code.
//**********************************************************************
// Set number of realizations.
N = 1;

// Set number of samples.
numberOfSamples = 1000;

// Generate unit variance white Gaussian noise.
X = generateGaussianProcess(N,numberOfSamples,1);

// Set MA(2) filter coefficients.
g = [1 1.8 0.8]';

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Generate reference signal, d(n).  This configures the
// system for a system idendification application.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
d = filterBlock(X(:,1),g,0);

// Run the normalized LMS algorithm.
[W,E] = nlms(X(:,1),d,0.1,3);

// Run the optimized normalized LMS algorithm.
[W1,E1] = nlms_optimized(X(:,1),d,0.1,3);

// Output the values.
printf("W (NLMS): %f  W1 (Optimized NLMS): %f\n",W($,1:$)',W1($,1:$)');

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Plot results.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
scf(1);

subplot(211);
title('NLMS Coefficient Trajectory');
plot(W);

subplot(212);
title('Optimized NLMS Coefficient Trajectory');
plot(W1);


