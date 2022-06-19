//**********************************************************************
// File Name: C8_8.sci
//**********************************************************************

// Bring in all of the utility functions.
exec('utils.sci',-1);

//**********************************************************************
// Mainline code.
//**********************************************************************
// Initialize white noise generator with unit variance.
noisegen(1,1024,1);

// Generate time vector.
n_100 = 0:99;
n_100 = n_100(:);

// Construct filter.
b2 = [1 1.5 1];
a2 = [1 0 0.9];

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (a): Generate 100 samples of an ARMA(2,2) process, x(n),
// by filtering unit variance white Gaussian noise with a filter
// that has a system function of the form,
// H(z) = {1 + 1.5z^(-1) + z^(-1)} / {1 + 0.9z^(-2).
// Using the method of iterative prefiltering, find the second-
// order ARMA model for x(n).
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

// Generate white noise sequence.
w = feval([1:100],Noise);

// Generate ARMA(2,2) process.
x = filterBlock(w,b2,a2(2:$));

// Estimate model parameters using 7 iterations.
[a,b,err] = ipf(x,2,2,7,[1 1]);

// Generate estimated ARMA(2,2) process.
xhat = filterBlock(w,b,a(2:$));

// Construct power spectra.
Px = constructPowerSpectrum(x);
Pxhat = constructPowerSpectrum(xhat);

// Convert to decibels.
Px = 20*log10(Px);
Pxhat = 20*log10(Pxhat);

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (b): Repeat your experiment in part (a) for 20 different
// realizations of x(n).  How accurate are the estimated model
// coefficients to the true model, on the average?  How accurate
// are the estimates of the power spectrum that are generated
// from these models?  What happens if the process is over-
// modeled by letting p = q = 3?  What about p = q = 5?
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Initialize accumulated sums.
bAvg = 0;
aAvg = 0;

// Preallocate matrices.
Px_20 = zeros(1024,1);
Pxhat_20 = zeros(1024,1);

for j = 1:20
  // Generate white noise sequence.
  w_20 = feval([1:100],Noise);

  // Generate ARMA(2,2) process.
  x_20 = filterBlock(w,b2,a2(2:$));

  // Estimate model parameters using 7 iterations.
  [a_20,b_20,err] = ipf(x_20,2,2,7,[1 1]);

  // Update accumulated sums.
  bAvg = bAvg + b_20;
  aAvg = aAvg + a_20;

  // Generate estimated ARMA(2,2) process.
  xhat_20 = filterBlock(w_20,b_20,a_20(2:$));

  // Construct power spectra.
  Px_20(:,j) = constructPowerSpectrum(x_20);
  Pxhat_20(:,j) = constructPowerSpectrum(xhat_20);
end

// Compute mean values.
bAvg = bAvg / 20;
aAvg = aAvg / 20;

// Convert to decibels.
Px_20 = 20*log10(Px_20);
Pxhat_20 = 20*log10(Pxhat_20);

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Plot results.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

scf(1);

// Part (a).
subplot(231);
title('Actual Spectrum');
plot(Px);

subplot(234);
title('Estimated Lpf Spectrum');
plot(Pxhat);

// Part (b).
subplot(232);
title('Actual Spectrum, Overlay: 20');
plot(Px_20);

subplot(235);
title('Estimated Lpf Spectrum, Overlay: 20');
plot(Pxhat_20);

subplot(233);
title('Average Actual Spectrum, Ensemble: 20');
plot(sum(Px_20,'c') / 20);

subplot(236);
title('Average Estimated Lpf Spectrum, Ensemble: 20');
plot(sum(Pxhat_20,'c') / 20);
