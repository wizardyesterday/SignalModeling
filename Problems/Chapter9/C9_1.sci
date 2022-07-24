//**********************************************************************
// File Name: C9_1.sci
//**********************************************************************

// Bring in all of the utility functions.
exec('utils.sci',-1);

//**********************************************************************
// Mainline code.
//**********************************************************************
// Set variances.
sigmavSq = 0.25;

// Initialize white noise generator.
noisegen(1,500,sqrt(sigmavSq));

// Set AR(2) filter coefficients.
a1 = [1 .1 .8]';
a2 = [1 .1 -.8]';

// Set step sizes.
mu1 = 0.05
mu2 = 0.01;

//----------------------------------------
// Generate AR(2) processes.
//----------------------------------------
for j = 1:2
  // Generate the noise.
  v = feval([1:500],Noise);

  X1(:,j) = filterBlock(v,1,a1(2:$));
  X2(:,j) = filterBlock(v,1,a2(2:$));
end

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (a): Paper problem.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (b): Plot the learning curve of the adaptive filter for
// 100 different realizations of x(n).  Estimate the misadjustment
// by time-averaging over the final iterations of the ensemble-
// average learning curves.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Initialize running sums.
W1avg_1 = 0;
W1avg_2 = 0;

// Generate error sequences.
for j = 1:2
  // Generate reference signal looking at this as an FIR filter.
  d = filterBlock(X1(:,j),a1(2:$),0);

  [W1_1,E1_1(:,j)] = lms(X1(:,j),d,mu1,2);
  [W1_2,E1_2(:,j)] = lms(X1(:,j),d,mu2,2);

  // Update running sums.
  W1avg_1 = W1avg_1 + W1_1($,1:2);
  W1avg_2 = W1avg_2 + W1_2($,1:2);
end

// Compute squared error.
E1Sq_1 = E1_1 .* E1_1;
E1Sq_2 = E1_2 .* E1_2;

// Compute learning curves.
for j = 401:500
  E1avg_1(j) = sum(E1Sq_1(j,:));
  E1avg_2(j) = sum(E1Sq_2(j,:));
end

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (c): Estimate the steady-state values of the adaptive
// filter coefficients for step sizes of mu = 0.05 and 0.01.
// Note that running sums were computed in Part (b).
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Compute average frequency estimates.
W1avg_1 = W1avg_1 / length(W1avg_1);
W1avg_2 = W1avg_2 / length(W1avg_2);

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (c): Repeat the experiments in parts (b) and (c) with
// a(1) = 0.1, a(2) = -0.8, and sigmaV^2 = 0.25.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Initialize running sums.
W2avg_1 = 0;
W2avg_2 = 0;

// Generate error sequences.
for j = 1:2
  // Generate reference signal looking at this as an FIR filter.
  d = filterBlock(X2(:,j),a2(2:$),0);

  [W2_1,E2_1(:,j)] = lms(X2(:,j),d,mu1,2);
  [W2_2,E2_2(:,j)] = lms(X2(:,j),d,mu2,2);

  // Update running sums.
  W2avg_1 = W2avg_1 + W2_1($,1:2);
  W2avg_2 = W2avg_2 + W2_2($,1:2);
end

// Compute squared error.
E2Sq_1 = E2_1 .* E2_1;
E2Sq_2 = E2_2 .* E2_2;

// Compute learning curves.
for j = 401:500
  E2avg_1(j) = sum(E2Sq_1(j,:));
  E2avg_2(j) = sum(E2Sq_2(j,:));
end

// Compute average frequency estimates.
W2avg_1 = W2avg_1 / length(W2avg_1);
W2avg_2 = W2avg_2 / length(W2avg_2);

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Plot results.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

//scf(1);

