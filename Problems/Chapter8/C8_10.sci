//**********************************************************************
// File Name: C8_10.sci
//**********************************************************************

// Bring in all of the utility functions.
exec('utils.sci',-1);

//**********************************************************************
//
//  Name:  minvar_auto
//
//  The purpose of this function is to estimate the spectrum of a
//  process using the minimum variance method.
//
//  Calling Sequence: Px = minvar_auto(r)
//
//  Inputs:
//
//    r - The input autocorrelation sequence.
//
//  Outputs:
//
//    Px - The minimum variance estimate of the power spectrum of
//    x(n) using a decibel scale.
//
//**********************************************************************
function Px = minvar_auto(r)

  // Compute the autocorrelation matrix.
  R = autocorrelationMatrix(r)

  // Compute the eigenvectors of R and the diagonal matrix of eigenvalues.
  [v,d] = spec(R);

  // Construct a column vector of eigenvalue reciprocals.
  U = diag(inv(abs(d)+ %eps));

  // Compute the 1024-pointspectrum for each column.
  V = matrixFft(v,1024);
  V = V .* conj(V);

  // Estimate the power spectrum in decibels.
  Px = (10 * log10(p)) - (10 * log10(V * U));

endfunction

//**********************************************************************
//
//  Name: acm2
//
//  Purpose: The purpose of this function is to compute an all pole
//  model of a signal using the autocorrelation method.
//  The signal is modeled as the unit sample response of a system
//  represented by, H(z) = b(0) / A(z), such that the coefficients of
//  A(z) are contained in the vector, [1 a(1), a(2),... a(p)].  The
//  difference between this function and the acm() function is that
//  we're working with the autocorrelation matrix rather than the data
//  matrix.  This allows to do clever trickery to the autocorrelation
//  sequence.
//
//  Calling Sequence: [a = acm2(x,p,sigma2)
//
//  Inputs:
//
//    x - The input vector to be processed.
//
//    p - The number of poles in the model.
//
//    sigma2 - The variance of the white noise that has been added
//    to the vector, x.  This results in a noisy autocorrelation
//    sequence.
//
//  Outputs:
//
//    a - The denominator coefficients in the model.
//
//    err - The model error.
//
//**********************************************************************
function [a,err] = acm2(x,p,sigma2)

  // Ensure that we have a column vector.
  x   = x(:);

  N   = length(x);

  if p < length(x)
    // Construct bipolar noise source.
    v = generateBipolarNoise(sqrt(3*sigma2),p+1);
    v = v(:);

    // Construct autocorrelation sequence.
    rx = convol(x,flipud(x));
    k = find(rx == max(rx));
    rx = rx(k:$);
    rx = rx(:);

    // Make the autocorrelation sequence noisy.
    rx = rx + v;
    
    // Construct autocorrelation matrix.
    Rx = autocorrelationMatrix(rx(1:p));

    // Compute denominator coefficients, a(1) ... a(p).
    a   = -Rx \ rx(2:p+1);

    // Insert a(0).
    a = [1; a];

    err = rx(1) + a(1:p)' * conj(rx(1:p));
  else
    error('Model order too large');
  end

endfunction

//**********************************************************************
//
//  Name:  mem2
//
//  The purpose of this function is to estimate the spectrum of a
//  process using the maximum entropy method.  The autocorrelation
//  method is used to find a pth-order all-pole model for x(n), and
//  the spectral estimate is formed from the following equation:
//
//    Px = b^2(0) / |A(omega)|^2
//
//  This routine is a modified varion of the mem() function so that
//  a variance can be passed to the modified version of the acm()
//  function.
//
//  Calling Sequence: Px = mem2(x,p,sigma2)
//
//  Inputs:
//
//    x - The input sequence.
//
//    p - The order of the all-pole model.
//
//    sigma2 - The variance of the white noise that has been added
//    to the vector, x.
//
//  Outputs:
//
//    Px - The maximum entropy estimate of the power spectrum of
//    x(n) using a decibel scale.
//
//**********************************************************************
function Px = mem2(x,p,sigma2)

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Construct an all-pole model using the autocorrelation 
  // matrix method.
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  [a,e] = acm2(x,p,sigma2);
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

  // Compute the FFT of the all-pole model.
  A = matrixFft(a,1024);

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Compute the estimated power spectral density.  Note that
  // e = b^(0), and that e is being normalized by the length
  // of the input vector.  Also, notice that the second log()
  // term is being multiplied by 20 to account that A
  // represents a voltage rather than a power (|A|^2),
  // Had we been working with |A|^2, the multiplier would have
  // been 10. 
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  Px = 10 * log10(e / length(x)) - 20 * log10(abs(A));
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

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
// Part (a): Compute the mimimum variance spectrum with p = 10
// and sigmaW^2 = 0.
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
