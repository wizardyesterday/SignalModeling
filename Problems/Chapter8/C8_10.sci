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
// Part (a): Compute the nimimum variance spectrum with p = 10
// and sigmaW^2 = 0.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Construct noise free autocorrelation sequence.
rx = cos(0.35*%pi*k);

// Add delta(k) to rx(0).
rx(1) = rx(1) + 1;

Px_a = minvar_auto(rx);

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Plot results.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

//scf(1);

