//**********************************************************************
// File Name: C9_2.sci
//**********************************************************************

// Bring in all of the utility functions.
exec('utils.sci',-1);

//**********************************************************************
// Mainline code.
//**********************************************************************
// Set number of realizations.
N = 1;

// Set variance.
sigmavSq = 1;

//------------------------------------------
// Generate Gaussian random process.  We
// want 100 realizations of length 1000.
// The number of realizations may be scaled
// back due to memory limitations.
//------------------------------------------
X = generateGaussianProcess(N,1000,sigmavSq);

// Set MA(2) filter coefficients.
g = [1 1.8 0.8]';

//+++++++++++++++++++++++++++++++++++++++++
// Generate autocorrelation sequence.
//+++++++++++++++++++++++++++++++++++++++++
rg = convol(g,g($:-1:1));
k = find(g == max(g));
rg = rg(k:$);

// Construct autocorrelation matrix.
Rg = toeplitz(rg);

// Compute eigenvalues.
lamda = spec(Rg);
lamdaMax = max(lamda);

// Set step size
muMax = 2/lamdaMax;
mu1 = 0.1 * muMax;
mu2 = 0.01 * muMax;
mu3 = 0.2 * muMax;

// Compute LMS misadjustments.
M1 = computeLmsMisadjustment(mu1,lamda);
M2 = computeLmsMisadjustment(mu2,lamda);
M3 = computeLmsMisadjustment(mu3,lamda);

// Generate reference signal, d(n).
d = filterBlock(X(:,1),g,0);

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (a): Paper problem.  The range of values for mu have
// already been computed.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (b): Implement an adaptive filter of order, p = 4, using
// the LMS algorithm.  Use a step size of mu = 0.1muMax.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Run LMS adaptive filter.
[W_lms_mu1,E_lms_mu1] = lms(X(:,1),d,mu1,4);

// Output result.
printf("W_1lms_mu1\n");
disp(W_lms_mu1($,1:4));

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (c): Repeat part (b) using the normalized LMS algorithm
// with beta = 0.1
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Run normalized LMS adaptive filter.
[W_nlms_mu1,E_nlms_mu1] = nlms(X(:,1),d,0.1,4);

// Output result.
printf("W_nlms_mu1\n");
disp(W_nlms_mu1($,1:4));

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Plot results.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

