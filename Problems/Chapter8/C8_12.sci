//**********************************************************************
// File Name: C8_12.sci
//**********************************************************************

// Bring in all of the utility functions.
exec('utils.sci',-1);

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
//  Calling Sequence: [B,Px] = modminvar(x,p)
//
//  Inputs:
//
//    x - The input sequence.
//
//    p - The order of the modified minimum variance estimate.
//
//  Outputs:
//
//    B - A vector of bandwidths of the bandpass filters.
//
//    Px - The minimum variance estimate of the power spectrum of
//    x(n) using a decibel scale.
//
//**********************************************************************
function P = computePowers(x,p)

  // Enforce a column vector.
  x = x(:);

  // First, find the pseudospectrum.
  [vmin,sigma] = phd(x,p)


  // Generate frequency index.
  n = 0:p-1;

  for j = 1:p
    for k = 0:p-1
      // Construct signal vector.
      e(k+1) = exp(%i*j*k*%pi/p);
    end

  end

  P = 1;

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
end

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (b): The power estimation function resides in this file.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (c): Repeat part(a) using the MUSIC algorithm, the
// engenvector method, and the minimum norm algorithm on 20
// different realizations of x(n).  Compare the accuracy of the
// estimates that are produced with each method.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/


//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Plot results.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/



