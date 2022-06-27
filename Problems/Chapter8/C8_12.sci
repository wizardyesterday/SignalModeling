//**********************************************************************
// File Name: C8_12.sci
//**********************************************************************

// Bring in all of the utility functions.
exec('utils.sci',-1);

//**********************************************************************
//
//  Name:  computePisarenkoPowers
//
//  The purpose of this function is to esimate the powers of complex
//  exponentials within a random process whose frequency estimates
//  are computed using the Pisarenko harmonic decomposition method.
//  This function performs all of the parameter computations that are
//  needed by the computeEigenPowers() function.
//
//  Calling Sequence: P = computePisarenkoPowers(x,p)
//
//  Inputs:
//
//    x - The input sequence.
//
//    p - The number of complex exponentials.
//
//  Outputs:
//
//    P - The power of each sinusoid.
//
//**********************************************************************
function P = computePisarenkoPowers(x,p)

  // Enforce a column vector.
  x = x(:);

  // Compute the autocorrelation matrix of order, p+1.
  R = Covar(x,p+1);

  // Compute the eigenvectors of R and the diagonal matrix of eigenvalues.
  [v,d] = spec(R);

  // Determine the variance of the noise.
  sigma2 = min(abs(diag(d)));

  // Find the eigenvector that spans the noise subspace.
  index = find(abs(diag(d)) == sigma2);

  vmin = v(:,index);

  // Remove vmin from the eigenvector matrix.
  v(:,index) = [];

  // Remove lamdaMin from the eigenvalue matrix.
  d(:,index) = [];
  d(index,:) = [];

  // Remove the insigificant values.
  d = clean(d);

  // Compute roots of vmin.
  ro = roots(vmin);

  // Compute frequency estimates.
  omega = atan(imag(ro),real(ro));

  // Construct vector of eigenvalues.
  lamda = diag(d);

  // Compute power estimates.
  P = computeEigenPowers(v,lamda,sigma2,omega);

endfunction

//**********************************************************************
// Mainline code.
//**********************************************************************
// Initialize unit variance noise source.
noisegen(1,100,1);

// Set up parameters
A1 = 4;
A2 = 3;
A3 = 1;
w1 = 0.4 * %pi;
w2 = 0.45 * %pi;
w3 = 0.8 * %pi;

// Generate time vector.
n = 0:99;

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (a): Assuming that x(n) contains three complex
// exponentials, use the PIsarenko harmonic decomposition to
// estimate the frequencies.  Repeat for 20 realizations of x(n).
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
for j = 1:20
  // Generate white noise process.
  w = feval([1:100],Noise);

  // Generate random phases.
  phi1 = generateBipolarNoise(%pi,1);
  phi2 = generateBipolarNoise(%pi,1);
  phi3 = generateBipolarNoise(%pi,1);

  // Generate realization.
  x = A1*exp(%i*w1*n + phi1) + + A2*exp(%i*w2*n + phi2) ...
      + A3*exp(%i*w3*n + phi3) + w;

  // Compute spectral estimate.
  [v(:,j),sigma] = phd(x,3);

  // Compute roots.
  ro = roots(v(:,j));

  // Estimate the frequencies of the complex exponentials.
  fr_phd(:,j) = atan(imag(ro),real(ro)) / %pi;

  // Compute the powers of the sinusoids.
  P(:,j) = computePisarenkoPowers(x,3);
end

for j = 1:3
  // Compute the variance of the frequency estimates.
  vr_phd(j) = variance(fr_phd(j,:));

  // Compute the mean values of the frequency estimates.
  mean_phd(j) = mean(fr_phd(j,:));
end

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (b): Estimate the powers of the signals using the
// frequency estimates derived in part (a).
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
for j = 1:3
  pMean(j,:) = mean(P(j,:));
end

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (c): Repeat part(a) using the MUSIC algorithm, the
// engenvector method, and the minimum norm algorithm on 20
// different realizations of x(n).  Compare the accuracy of the
// estimates that are produced with each method.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
for j = 1:20
  // Generate white noise process.
  w = feval([1:100],Noise);

  // Generate random phases.
  phi1 = generateBipolarNoise(%pi,1);
  phi2 = generateBipolarNoise(%pi,1);
  phi3 = generateBipolarNoise(%pi,1);

  // Generate realization.
  x = A1*exp(%i*w1*n + phi1) + + A2*exp(%i*w2*n + phi2) ...
      + A3*exp(%i*w3*n + phi3) + w;

  // Compute spectral estimates using 10 noise eigenvectors.
  Px_music(:,j) = music(x,3,10);
  Px_ev(:,j) = ev(x,3,10);
  Px_min_norm(:,j) = min_norm(x,3,10);
end

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Plot results.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

// Part (c).
subplot(311);
title('MUSIC, p:3, M:10');
plot(Px_music);

subplot(312);
title('Ev, p:3, M:10');
plot(Px_ev);

subplot(313);
title('Min_norm, p:3, M:10');
plot(Px_min_norm);



