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

// Initialize white noise generator with unit variance.
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
// Generate error sequences.
for j = 1:2
  // Generate reference signal looking at this as an FIR filter.
  d1 = filterBlock(X1(:,j),a1(2:$),0);

  [W1_1,E1_1(:,j)] = lms(X1(:,j),d1,mu1,2);
  [W1_2,E1_2(:,j)] = lms(X1(:,j),d1,mu2,2);
end

// Compute learning curves.
for j = 1:500
  E1avg_1(j) = sum(E1_1(j,:));
  E1avg_2(j) = sum(E1_2(j,:));
end

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (c): 
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Plot results.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

//scf(1);

