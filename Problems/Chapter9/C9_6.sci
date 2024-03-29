//**********************************************************************
// File Name: C9_6.sci
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
a1 = [1 1.2728 -0.81]';

// Compute autocorrelation.
r_a1 = ator(a1,1);

// Construct autocorrelation matrix.
R_a1 = toeplitz(r_a1);

// Compute eigenvalues.
lamda_a1 = spec(R_a1);

// Compute maximum LMS step size.
muMax = 2 / max(lamda_a1);

// Set typical values of the step size.
mu1 = muMax / 10;
mu2 = muMax / 25;

// Set exponential weighting factors.
lamda1 = 1
lamda2 = 0.99;
lamda3 = 0.95;
lamda4 = 0.92;
lamda5 = 0.90;

//----------------------------------------
// Generate AR(2) processes.
//----------------------------------------
for j = 1:N
  X1(:,j) = filterBlock(X(:,j),1,-a1(2:$));
end

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (a): Implement an RLS adaptive predictor with lamda = 1
// (growing window RLS), and plot w_n(k) for k = 1,2.  Compare
// the convergence of the coefficients, w_n(k), to those that
// are obtained using the LMS algorithm for several different
// values of the step size, mu.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (b): Make a plot of the learning curve for RLS, and
// compare to the LMS learning curve (See Example 9.2.2 on how
// to plot learning curves).
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

// Generate error sequences.
for j = 1:N
  // Generate x(n-1).
  xnm1 = filterBlock(X1(:,j),[0 1],0);

  // Run second-order RLS adaptive filter.
  [W_a1_lamda1_p2,E_a1_lamda1_p2(:,j)] = rls(xnm1,X1(:,j),2,lamda1);

  // run LMS adaptive filters.
  [W_a1_mu1_p2,E_a1_mu1_p2(:,j)] = lms(xnm1,X1(:,j),mu1,2);
  [W_a1_mu2_p2,E_a1_mu2_p2(:,j)] = lms(xnm1,X1(:,j),mu2,2);
end

// Compute squared errors.
ESq_a1_lamda1_p2 = E_a1_lamda1_p2 .* E_a1_lamda1_p2;
ESq_a1_mu1_p2 = E_a1_mu1_p2 .* E_a1_mu1_p2;
ESq_a1_mu2_p2 = E_a1_mu2_p2 .* E_a1_mu2_p2;

// Compute ensemble averages.  This is the learning curve.
for j = 1:numberOfSamples
 ESqAvg_a1_lamda1_p2(j) = mean(ESq_a1_lamda1_p2(j,:)); 
 ESqAvg_a1_mu1_p2(j) = mean(ESq_a1_mu1_p2(j,:)); 
 ESqAvg_a1_mu2_p2(j) = mean(ESq_a1_mu2_p2(j,:)); 
end


//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (c): Repeat part (b) for exponential weighting factors of
// lamda = 0.99, 0.95, 0.92, 0.90, and discuss the trade-offs
// involved in the choice of lamda.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Generate error sequences.
for j = 1:N
  // Generate x(n-1).
  xnm1 = filterBlock(X1(:,j),[0 1],0);

  // Run second-order RLS adaptive filters.
  [W_a1_lamda2_p2,E_a1_lamda2_p2(:,j)] = rls(xnm1,X1(:,j),2,lamda2);
  [W_a1_lamda3_p2,E_a1_lamda3_p2(:,j)] = rls(xnm1,X1(:,j),2,lamda3);
  [W_a1_lamda4_p2,E_a1_lamda4_p2(:,j)] = rls(xnm1,X1(:,j),2,lamda4);
  [W_a1_lamda5_p2,E_a1_lamda5_p2(:,j)] = rls(xnm1,X1(:,j),2,lamda5);
end

// Compute squared errors.
ESq_a1_lamda2_p2 = E_a1_lamda2_p2 .* E_a1_lamda2_p2;
ESq_a1_lamda3_p2 = E_a1_lamda3_p2 .* E_a1_lamda3_p2;
ESq_a1_lamda4_p2 = E_a1_lamda4_p2 .* E_a1_lamda4_p2;
ESq_a1_lamda5_p2 = E_a1_lamda5_p2 .* E_a1_lamda5_p2;

// Compute ensemble averages.  This is the learning curve.
for j = 1:numberOfSamples
 ESqAvg_a1_lamda2_p2(j) = mean(ESq_a1_lamda2_p2(j,:)); 
 ESqAvg_a1_lamda3_p2(j) = mean(ESq_a1_lamda3_p2(j,:)); 
 ESqAvg_a1_lamda4_p2(j) = mean(ESq_a1_lamda4_p2(j,:)); 
 ESqAvg_a1_lamda5_p2(j) = mean(ESq_a1_lamda5_p2(j,:)); 
end

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (d): The function, rls_slidingWindow() has been added to
// the AdaptiveFilters.sci file.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (e): Let x(n) = a_n(1)x(n-1) - 0.81x(n-2) + v(n),
// where v(n) is unit variance white noise, and a_n(1) is a
// time-varying coefficient given by,
// an_1) = 1.2728; 0 <= n < 50; 100 <= n <= 200, 0; 50 <= n < 100.
// Compare the effectiveness of the LMS, growing window RLS,
// exponentially weighted RLS, and sliding window RLS algorithms
// for adaptive linear prediction.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Generate unit variance white noise.
v = generateGaussianProcess(1,200,1);

a_n = [1.2728*ones(1,50) zeros(1,50) 1.2728*ones(1,100)];

// Instantiate the pipelines.
aState = zeros(1,2);
bState = zeros(1,1);

// Generate AR(2) process using time-varying coefficients.
for j = 1:200
  // Set the coefficients.
  a = [a_n(j) -0.81];

  // Run the sample through the filter.
  [X2(j),bState,aState] = iirFilter(v(j),1,-a,bState,aState);
end

// Generate x(n-1).
xnm1 = filterBlock(X2,[0 1],0);

// Run the second-order RLS adaptive filters.
[W_a2_lamda1_p2,E_dummy] = rls(xnm1,X2,2,lamda1);
[W_a2_lamda2_p2,E_dummy] = rls(xnm1,X2,2,lamda2);
[W_a2_lamda3_p2,E_dummy] = rls(xnm1,X2,2,lamda3);
[W_a2_lamda4_p2,E_dummy] = rls(xnm1,X2,2,lamda4);
[W_a2_lamda5_p2,E_dummy] = rls(xnm1,X2,2,lamda5);

// Run the LMS filter.
[W_a2_mu_p2,E_dummy] = lms(xnm1,X2,muMax/150,2);

// Run the second-order sliding window RLS filter with a window of 30 samples.
Wrls_sw = rls_slidingWindow(xnm1,X2,2,30);

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Plot results.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
scf(1);

// Part (a).
subplot(321);
title('RLS Coefficient Trajectory, Lamda = 1, a1, p = 2');
plot(W_a1_lamda1_p2);

subplot(323);
title('LMS Coefficient Trajectory, mu = muMax / 10, a1, p = 2');
plot(W_a1_mu1_p2);

subplot(325);
title('LMS Coefficient Trajectory, mu = muMax / 25, a1, p = 2');
plot(W_a1_mu2_p2);

// Part (b).
subplot(322)
title('RLS Learning Curve, Lamda = 1, a1, p = 2');
plot(ESqAvg_a1_lamda1_p2);

subplot(324)
title('LMS Learning Curve, mu = muMax / 10, a1, p = 2');
plot(ESqAvg_a1_mu1_p2);

subplot(326)
title('LMS Learning Curve, mu = muMax / 25, a1, p = 2');
plot(ESqAvg_a1_mu2_p2);

scf(2);

// Part (c).
subplot(421);
title('RLS Coefficient Trajectory, Lamda = 0.99, a1, p = 2');
plot(W_a1_lamda2_p2);

subplot(423);
title('RLS Coefficient Trajectory, Lamda = 0.95, a1, p = 2');
plot(W_a1_lamda3_p2);

subplot(425);
title('RLS Coefficient Trajectory, Lamda = 0.92, a1, p = 2');
plot(W_a1_lamda4_p2);

subplot(427);
title('RLS Coefficient Trajectory, Lamda = 0.90, a1, p = 2');
plot(W_a1_lamda5_p2);

subplot(422);
title('RLS Learning Curve, Lamda = 0.99, a1, p = 2');
plot(ESqAvg_a1_lamda2_p2);

subplot(424);
title('RLS Learning Curve, Lamda = 0.95, a1, p = 2');
plot(ESqAvg_a1_lamda3_p2);

subplot(426);
title('RLS Learning Curve, Lamda = 0.92, a1, p = 2');
plot(ESqAvg_a1_lamda4_p2);

subplot(428);
title('RLS Learning Curve, Lamda = 0.90, a1, p = 2');
plot(ESqAvg_a1_lamda5_p2);

scf(3);

// Part (e).
subplot(421);
title('RLS Coefficient Trajectory, Lamda = 1, a2, p = 2');
plot(W_a2_lamda1_p2);

subplot(422);
title('RLS Coefficient Trajectory, Lamda = 0.99, a2, p = 2');
plot(W_a2_lamda2_p2);

subplot(423);
title('RLS Coefficient Trajectory, Lamda = 0.95, a2, p = 2');
plot(W_a2_lamda3_p2);

subplot(424);
title('RLS Coefficient Trajectory, Lamda = 0.92, a2, p = 2');
plot(W_a2_lamda4_p2);

subplot(425);
title('RLS Coefficient Trajectory, Lamda = 0.90, a2, p = 2');
plot(W_a2_lamda5_p2);

subplot(426);
title('LMS Coefficient Trajectory, mu = muMax/150, a2, p = 2');
plot(W_a2_mu_p2);

subplot(427);
title('Sliding Window RLS Coefficient Trajectory, L = 30, a2, p = 2');
plot(Wrls_sw);






