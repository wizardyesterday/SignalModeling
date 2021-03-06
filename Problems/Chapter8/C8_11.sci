//**********************************************************************
// File Name: C8_11.sci
//**********************************************************************

// Bring in all of the utility functions.
exec('utils.sci',-1);

//**********************************************************************
//
//  Name:  modminvar
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
function [B,Px] = modminvar(x,p)

  // Enforce a column vector.
  x = x(:);

  // Compute the autocorrelation matrix of order, p.
  R = Covar(x,p);

  // Compute inverse of autocorrelation matrix.
  Rinv = inv(R);

  // Square it.
  Rinv2 = Rinv * Rinv;

  // Generate frequency index.
  n = 0:p-1;

  for j = 1:p
    for k = 0:p-1
      // Construct signal vector.
      e(k+1) = exp(%i*j*k*%pi/p);
    end

    // Compute filter.
    g = (Rinv * e)  / (e' * Rinv * e);

    // compute bandwidth.
    B(j) = g' * g;

    // Compute power spectrum estimate.
    Px(j) = (10 * log10(e' * Rinv * e)) - (10*  log10(e' * Rinv2 * e));
  end

  // Ensure that we have real values.
  B = real(B);
  Px = real(Px);

endfunction

//**********************************************************************
// Mainline code.
//**********************************************************************

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (a): Pencil and paper problem.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (b): The written function is modminvar().
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (c): Pencil and paper problem.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (d): Let x(n) be a harmonic process consisting of a sum
// of two sinusoids in white noise.  Let the sinusoid frequencies
// be w1 = 0.25PI and w2=0.25PI, and assume that the signal-to-
// noise ratio for each sinusoid is 10dB.  Find the
// autocorrelation sequence of this process, and compute the
// 40th-order minimum variance estimate and a 40th-order MEM
// estimate.  Note that the problem stated to use 10th-order
// models, but the results were unsatisfactory for analysis of
// a record length of 256.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Initialize unit variance noise source.
noisegen(1,256,1);

// Construct noise.
v = feval([1:256],Noise);

// Generate time index.
n = 0:255;

// Generate noisy signal.
x = sqrt(20)*cos(0.2*%pi*n) + sqrt(20)*cos(0.25*%pi*n) + v;

// Compute modified minimum variance spectrum.
[B,Px_modminvar_d] = modminvar(x,40);

// Compute minimum variance spectrum.
Px_minvar_d = minvar(x,40);

// Compute maximum entropy spectrum.
Px_mem_d = mem(x,40);

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (e): Make a plot of the bandwidths of the bandpass
// filters versus omega that are used in the modified MV
// spectrum estimate generated in (d), and compare them to the
// fixed bandwidth, deltaOmega = 2*PI/p, used in the MV spectrum
// estimate.  The bandwidths were already computed in part (d)
// using the modminvar() function.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (f): Create an ensemble of 25 different harmonic
// processes of the form given in (d), and generate an overlay
// plot of the modified MV spectrum estimates.  Repeat for the MV
// spectrum estimate, and comment on any differences that you
// observe.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

for j = 1:25
  // Construct noise.
  v = feval([1:256],Noise);

  // Generate noisy signal.
  x = sqrt(20)*cos(0.2*%pi*n) + sqrt(20)*cos(0.25*%pi*n) + v;

  // Compute modified minimum variance spectrum.
  [Bdummy,Px_modminvar_f(:,j)] = modminvar(x,40);

  // Compute minimum variance spectrum.
  Px_minvar_f(:,j) = minvar(x,40);

  // Compute maximum entropy spectrum.
  Px_mem_f(:,j) = mem(x,40);
end

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Plot results.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

scf(1);

// Part (d).
subplot(421);
title('Modified MV Spectrum');
plot(Px_modminvar_d);

subplot(423);
title('MV Spectrum');
plot(Px_minvar_d);

subplot(425);
title('MEM Sectrum');
plot(Px_mem_d);

// Part (e).
subplot(427);
title('Bandwidth of Bandpass Filters');
plot(B);

// Part (f).
subplot(422);
title('Modified MV Spectrum Overlay');
plot(Px_modminvar_f);

subplot(424);
title('MV Spectrum Overlay');
plot(Px_minvar_f);

subplot(426);
title('MEM Sectrum Overlay');
plot(Px_mem_f);
