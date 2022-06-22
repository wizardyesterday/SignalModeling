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
//  process using the modified minimum variance method.
//
//  Calling Sequence: Px = modminvarminvar(x,p)
//
//  Inputs:
//
//    x - The input sequence.
//
//    p - The order of the minimum variance estimate.  For short
//    sequences, p is typically about length(x) / 3.
//
//  Outputs:
//
//    Px - The minimum variance estimate of the power spectrum of
//    x(n) using a decibel scale.
//
//**********************************************************************
function Px = modminvar(x,p)

  // Enforce a column vector.
  x = x(:);

  // Compute the autocorrelation matrix of order, p.
  R = Covar(x,p);

  // Compute the eigenvectors of R and the diagonal matrix of eigenvalues.
  [v,d] = spec(R);

  // Construct a column vector of eigenvalue reciprocals.
  U = diag(inv(abs(d)+ %eps));

  // Compute the 1024-pointspectrum for each column.
  V = matrixFft(v,1024);
  V = V .* conj(V);

  // Construct numerator.
  num = V * U;

  // Construct denominator.
  den = num .* num;

  // Estimate the power spectrum in decibels.
  Px = (10 * log10(num)) - (10 * log10(den));

endfunction

//**********************************************************************
// Mainline code.
//**********************************************************************
// Set model order.
p = 10;

// Construct time vector.
k = 0:p;

// Construct autocorrelation sequence.
rx = cos(0.35*%pi*k);

// Add delta(k) to rx(0).
rx(1) = rx(1) + 1;

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (a): The written function is modminvar().
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Construct noise free autocorrelation sequence.
rx = cos(0.35*%pi*k);
rx = rx(:);

// Add delta(k) to rx(0).
rx(1) = rx(1) + 1;

// Compute minimum variance spectrum.
Px_a = minvar_auto(rx);

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (b): Generate an overlay plot of 25 minimum variance
// spectrum estimates with sigmaW^2 = 0.1, and compare these
// estimates to that generated in part (a).
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Preallocate matrix.
Py_b = zeros(1024,1);

// Set variance.
sigmaW2 = 0.1;

// Compute limits of uniform distribution.
b = sqrt(3*sigmaW2);

for j = 1:25
  w = generateBipolarNoise(b,p+1);

  // Generate noisy autocorrelation sequence.
  ry = rx + w;

  // Compute minimum variance spectrum.
  Py_b(:,j) = minvar_auto(ry);
end

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (c1): Compute the maximum entropy spectrum with p = 10
// and sigmaW^2 = 0.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Initialize unit variance noise source.
noisegen(1,256,1);

// Construct noise.
v = feval([1:11],Noise);

// Generate time index.
n = 0:10;

// Generate noisy signal.
x = cos(0.35*%pi*n); + v;

// Compute maximum entropy spectrum.
Px_c1 = mem2(x,10,0);

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (b): Generate an overlay plot of 25 maximum entropy
// spectrum estimates with sigmaW^2 = 0.1, and compare these
// estimates to that generated in part (a).
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Preallocate matrix.
Py_c2 = zeros(1024,1);

// Set variance.
sigmaW2 = 0.1;

for j = 1:25
  // Construct noise.
  v = feval([1:11],Noise);

  // Generate noisy signal.
  x = cos(0.35*%pi*n); + v;

  // Compute maximum entropyspectrum.
  Py_c2(:,j) = mem2(x,10,0.1);
end

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Plot results.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

scf(1);

// Part (a).
subplot(221);
title('Noise Free  Minvar Spectrum');
plot(Px_a);

// Part (b).
subplot(223);
title('Minvar Spectrum, Variance: 0.1');
plot(Py_b);

// Part (c).
subplot(222);
title('Noise Free Mem Spectrum');
plot(Px_c1);

subplot(224);
title('Mem Spectrum, Variance: 0.1');
plot(Py_c2);
