//**********************************************************************
// File Name: C9_7.sci
//**********************************************************************

// Bring in all of the utility functions.
exec('utils.sci',-1);

//**********************************************************************
// Mainline code.
//**********************************************************************
// Set number of realizations.
N = 30;

// Set number of samples.
numberOfSamples = 500;

// Set variance.
sigmavSq = 1;

//------------------------------------------
// Generate Gaussian random process.  We
// want 30 realizations of length 500.
// The number of realizations may be scaled
// back due to memory limitations.
//------------------------------------------
X = generateGaussianProcess(N,numberOfSamples,sigmavSq);

// Coefficients for a1.
a1 = [1 0.8 -0.9]';

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

//----------------------------------------
// Generate AR(2) processes.
//----------------------------------------
for j = 1:N
  X1(:,j) = filterBlock(X(:,j),1,-a1(2:$));
end

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (a): Generate a sequence, x(n), of length N = 500, using
// the difference equation, x(n) = 0.8x(n-1) - 0.9x(n-2) + v(n),
// where v(n) is unit variance white Gaussian noise.  Determine
// a value for the LMS step size so that w_n converges in the
// mean to w = [0.8 -0.9]' within about 200 iterations.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (b): Implement the LMS adaptive predictor, and plot
// w_n(k) versus n for k = 1,2.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (c): Estimate the mean-square error, E(infinity), and
// compare it to the theoretical value.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

// Generate multiple realizations.
for j = 1:N
  // Generate x(n-1).
  xnm1 = filterBlock(X1(:,j),[0 1],0);

  // run LMS adaptive filter.
  [W_a1_mu1_p2,E_a1_mu1_p2(:,j)] = lms(xnm1,X1(:,j),mu1,2);
end

// Compute squared error.
ESq_a1_mu1_p2 = E_a1_mu1_p2 .* E_a1_mu1_p2;

// Compute ensemble averages.  This is the learning curve.
for j = 1:numberOfSamples
 ESqAvg_a1_mu1_p2(j) = mean(ESq_a1_mu1_p2(j,:)); 
end

// Estimate the steady-state mean-square error.
EInf_a1_mu1_p2 = mean(ESqAvg_a1_mu1_p2(numberOfSamples-99:numberOfSamples));

// Compute the theoretical steady-state mean-square error
EInfTheor_a1_mu1_p2 = sigmavSq / (1 - mu1 * trace(R) / 2);

mprintf("\nLMS E(infinity): %f\n",EInf_a1_mu1_p2);
mprintf("\nLMS Theoretical E(infinity): %f\n",EInfTheor_a1_mu1_p2);

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Plot results.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
//scf(1);

