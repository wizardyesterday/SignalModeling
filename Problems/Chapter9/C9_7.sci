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
// Part (d): Paper exercise.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (e): Find the vector, rdx = E{d(n)x(n)}, for a one-step
// linear predictor, and implement the p-vector algorithm for
// x(n).  Plot w_n(k) versus n, and compare your results to the
// LMS adaptive linear predictor.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Set filter order.
p = 2;

// Generate unit variance white Gaussian noise.
v = generateGaussianProcess(1,numberOfSamples,1);

// Set desired signal, the AR(2) process.
d = filterBlock(v,1,-a1(2:$));

// Generate x(n-1).
xnm1 = filterBlock(d,[0 1],0);

for n = 1:numberOfSamples
  //--------------------------------------------------
  // Construct the next vector, whose entries are,
  // rdx(k) = E{d(n)xnim1(n-k)}; k = 0, ..., p-1.
  //--------------------------------------------------
  for k = 0:p-1 
    rdx(k + 1) = crosscorrelate(d(n:$),xnm1(n:$),numberOfSamples,k);
  end

  //--------------------------------------------------
  // Add rdx_ = [rdx(0) rdx(1) ... rdx(p-1)]' to the
  // next row of Rdx.  Note that rdx_ is a vector of
  // cross-correlations, thus, Rdx is a matrix whose
  // rows are these vectors.
  //--------------------------------------------------
  Rdx(n,:) = rdx';
end

// Run p-vector adaptive filter.
W = lms_pvector(xnm1,Rdx,mu1/3,p);

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (f): Investigate the sensitivity of the p-vector
// algorithm to errors in the vector rdx by making small changes
// in the values of rdx used in part (e).
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Perturb the values of rdx.
for j = 1:numberOfSamples
  for k = 1:p
    RdxPerturbed(j,k) = Rdx(j,k) + 0.5 * (-1)^k;
  end
end

// Run p-vector adaptive filter on the perturbed reference.
WPerturbed = lms_pvector(xnm1,RdxPerturbed,mu1/3,p);

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Plot results.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
scf(1);

subplot(411);
title('LMS Coefficients, p=2, a_ = [1 0.8 -0.9], mu=muMax / 800');
plot(W_a1_mu1_p2);

subplot(412);
title('LMS Learning Curve, p=2, a_ = [1 0.8 -0.9], mu=muMax / 800');
plot(ESqAvg_a1_mu1_p2);

subplot(413);
title('P-vector Coefficients, p=2, a_ = [1 0.8 -0.9], mu=muMax / 2400');
plot(W);

subplot(414);
title('P-vector Coefficients, Noisy, p=2, a_ = [1 0.8 -0.9], mu=muMax / 2400');
plot(WPerturbed);


