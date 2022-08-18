//**********************************************************************
// File Name: C9_7.sci
//**********************************************************************

// Bring in all of the utility functions.
exec('utils.sci',-1);

//**********************************************************************
// Mainline code.
//**********************************************************************
// Set number of samples.
numberOfSamples = 500;

// Generate unit variance white Gaussian noise. random process.
X = generateGaussianProcess(1,numberOfSamples,1);

// Coefficients for a1.
a1 = [1 1.5 -0.8]';

// Compute autocorrelation.
r = ator(a1,1);

// Construct autocorrelation matrix.
R = toeplitz(r);

// Compute eigenvalues.
lamda = spec(R);

// Compute maximum LMS step size.
muMax = 2 / max(lamda);

// Choise required step size.
mu1 = muMax / 800;

// Set constant coefficients.
c1_1 = 0.2;
c2_1 = .6;

// Generate AR(2) process.
X1 = filterBlock(X,1,-a1(2:$));

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (a): 
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Generate x(n-1).
xnm1 = filterBlock(X1,[0 1],0);

// Run variable step size LMS filter.
[W_c1,E_c1] = lms_variableStepSize(xnm1,X1,c1_1,1,2);

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (a): 
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Run variable step size LMS filter.
[W_c1c2,E_c1c2] = lms_variableStepSize(xnm1,X1,c1_1,c2_1,2);

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Plot results.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
scf(1);

subplot(211);
title('LMS Coefficient Trajectory, c1 = 50, c2 = 1');
plot(W_c1);

subplot(212);
title('LMS Coefficient Trajectory, c1 = 50, c2 = 0.1');
plot(W_c1c2);
