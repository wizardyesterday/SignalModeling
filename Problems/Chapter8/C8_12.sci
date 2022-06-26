//**********************************************************************
// File Name: C8_12.sci
//**********************************************************************

// Bring in all of the utility functions.
exec('utils.sci',-1);

//**********************************************************************
//
//  Name:  computeDtft
//
//  The purpose of this function is to compute the discrete Fourier
//  transform of a short sequence.
//
//  Calling Sequence: V = computeDtft(v,p,w)
//
//  Inputs:
//
//    V - The input sequence. This is typically an eigenvector with
//    a length of p+1.
//
//    p - The order of the DFT.
//
//    w - The radian frequency for which the DTFT is to be evaluated.
//
//  Outputs:
//
//    V - The discrete time Fourier transform of v.
//
//**********************************************************************
function V = computeDtft(v,p,w)

  // Clear sum.
  V = 0;

  // Compute the DTFT of v.
  for k = 0:p
    V = V + v(k+1) * exp(-%i*k*w);
  end

endfunction

//**********************************************************************
//
//  Name:  computePowers
//
//  The purpose of this function is to estimate the spectrum of a
//  process using the modified covariance method. Additionally, this
//  function computeds the bandwidth of the bandpass filters that are
//  "designed" by the modified minimum  ariance method.  Note that
//  this function provides a good mapping of the equations described
//  in the textbook.
//
//  Calling Sequence: P = computePowers(x,p)
//
//  Inputs:
//
//    x - The input sequence.
//
//    p - The order of the modified minimum variance estimate.
//
//  Outputs:
//
//    P - The power of each sinusoid.
//
//**********************************************************************
function P = computePowers(x,p)

  // Enforce a column vector.
  x = x(:);

  // Compute the autocorrelation matrix of order, p+1.
  R = Covar(x,p+1);

  // Compute the eigenvectors of R and the diagonal matrix of eigenvalues.
  [v,d] = spec(R);

  // Determine the variance of the noise.
  sigma = min(abs(diag(d)));

  // Find the eigenvector that spans the noise subspace.
  index = find(abs(diag(d)) == sigma);

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

  for i = 1:p
    for k = 1:p
      // Compute V_i(omega_k).
      V(i,k) = computeDtft(v(:,i),p,omega(k));

      // Compute |V_i(omega_k)|^2.
      V(i,k) = V(i,k) .* conj(V(i,k));
    end
  end

  // Construct right hand side vector of eigenvalues.
  lam = diag(d);

  // Subtract the noise.
  lam = lam - sigma;

  // Compute power estimates.
  P = V \ lam;

  // Remove infinitesimal imaginary component.
  P = real(P);

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
  P(:,j) = computePowers(x,3);
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



