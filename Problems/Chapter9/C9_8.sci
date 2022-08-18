//**********************************************************************
// File Name: C9_8.sci
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
c1_1 = 0.01;
c1_2 = 0.2;
c1_3 = 50;
c2_1 = .6;
c2_2 = 10;
c2_3 = 50;

// Generate AR(2) process.
X1 = filterBlock(X,1,-a1(2:$));

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (a): Evaluate the effectiveness of the variable step LMS
// algorithm with mu(n) = 1 / (c + n).  Choose a reasonable
// value for c, and determine the limitations of this algorithm.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Generate x(n-1).
xnm1 = filterBlock(X1,[0 1],0);

// Run variable step size LMS filter.
[W_c1_1,E_c1_1] = lms_variableStepSize(xnm1,X1,c1_1,1,2);
[W_c1_2,E_c1_2] = lms_variableStepSize(xnm1,X1,c1_2,1,2);
[W_c1_3,E_c1_3] = lms_variableStepSize(xnm1,X1,c1_3,1,2);

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (b): Evaluate the effectiveness of the variable step LMS
// algorithm with mu(n) = 1 / (c1 + c2*n).  Note that this allows
// for better control of how fast the step size decreases to
// zero.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Run variable step size LMS filter.
[W_c1_1c2_1,E_c1_1c2_1] = lms_variableStepSize(xnm1,X1,c1_1,c2_1,2);
[W_c1_1c2_2,E_c1_1c2_2] = lms_variableStepSize(xnm1,X1,c1_1,c2_2,2);
[W_c1_1c2_3,E_c1_1c2_3] = lms_variableStepSize(xnm1,X1,c1_1,c2_3,2);

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Plot results.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
scf(1);

subplot(321);
title('LMS Coefficient Trajectory, c1 = 0.01, c2 = 1');
plot(W_c1_1);

subplot(323);
title('LMS Coefficient Trajectory, c1 = 0.2, c2 = 1');
plot(W_c1_1);

subplot(325);
title('LMS Coefficient Trajectory, c1 = 50, c2 = 1');
plot(W_c1_1);

subplot(322);
title('LMS Coefficient Trajectory, c1 = 0.2, c2 = 0.6');
plot(W_c1_1c2_1);

subplot(324);
title('LMS Coefficient Trajectory, c1 = 0.2, c2 = 10');
plot(W_c1_1c2_2);

subplot(326);
title('LMS Coefficient Trajectory, c1 = 0.2, c2 = 50');
plot(W_c1_1c2_3);

