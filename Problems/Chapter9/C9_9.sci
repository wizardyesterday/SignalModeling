//**********************************************************************
// File Name: C9_9.sci
//**********************************************************************

// Bring in all of the utility functions.
exec('utils.sci',-1);

//**********************************************************************
//
//  Name: nlms_noiseCanceller
//
//  The purpose of this function is to perform noise cancelling using
//  the normalized LMS algorithm.  The reference signal is formed by
//  delaying the input signal by a specified number of samples such
//  that the noise portion of the delayed input signal is uncorrelated
//  with the nondelayed input signal.  This way, the adaptive filter
//  of the noise canceller will prodide an estimate of the desired
//  signal as its output.
//
//  Calling Sequence: [W,dyat] = nlms_noiseCanceller(x,n0,Beta,nord,w0)
//
//  Inputs:
//
//    x - The input data to the adaptive filter.
//
//    n0 - The number of samples to delay the input data, x, so
//    that the reference signal, d(n) = x(n - n0), can be formed.
//
//    Beta - The adaptive filtering update (normalized step-size)
//    parameter.
//
//    nord - The number of filter coefficients.
//
//    w0 - An optional row vector that serves as the initial guess
//    for FIR filter coefficients.  If w0 is omitted, then w0 = 0 is
//    assumed.
//    
//  Outputs:
//
//    W - A matrix of filter coefficients as they evolve over time.
//    Each row of this matrix contains the coefficients at the
//    iteration that is associated with the row.  For example, row 1
//    contains the coefficients for iteration 1, row 2 contains the
//    coefficients for iteration 2, and so on.
//
//    dhat - A vector of errors as they evolve over time.  Each entry
//    contains the error that is associated with each iteration.  For
//    example, the first entry contains the error for iteration 1,
//    the second entry contains the error for iteration 2, and so on.
//    Note that E(n) = d(n) - dhat(n).
//
//**********************************************************************
function [W,dhat] = nlms_noiseCanceller(x,n0,Beta,nord,w0)

  // Create filter coefficients for the delay line.
  b = zeros(1,n0+1);
  b($) = 1;

  // Construct the reference signal, d(n) = x(n-n0).
  d = filterBlock(x,b,0);

  // Construct the data matrix.
  X = convm(x,nord);

  // Retrieve the size of the data matrix.
  [M,N] = size(X);

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // argn(2) returns is the number of arguments passed to the 
  // function.
  // If 4 arguments were passed to the function, it is implied
  // that the last parameter was not passed.  In this
  // case, the initial condition for the filter coefficients
  // is set to a default value of all zeros. 
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  if argn(2) < 5
    w0 = zeros(1,N);
  end

  // Force a row vector without altering the values.
  w0 = w0(:).';

  // Compute the first output value.
  dhat(1) = w0 * X(1,:).';

  // Compute the first iteration for the error.
  E(1) = d(1) - dhat(1); 

  // Construct the normalizing denominator for the first iteration.
  DEN = X(1,:) * X(1,:)' + 0.0001;

  // Perform the first iteration for the filter vector.
  W(1,:) = w0 + (Beta / DEN) * E(1) * conj(X(1,:));

  if M > 1
    for k = 2:M - nord + 1
      // Compute the next output.
      dhat(k) = W(k-1,:) * X(k,:).';

      // Compute the next error.
      E(k) = d(k) - dhat(k);

      // Update the normalizing denominator.
      DEN = X(k,:) * X(k,:)' + 0.0001;

      // Update the filter vector.
      W(k,:) = W(k-1,:) + (Beta / DEN) * E(k) * conj(X(k,:));
    end
  end

endfunction

//**********************************************************************
// Mainline code.
//**********************************************************************
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (a): Paper exercise.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (b): Create NLMS noise canceller.  Given an input signal,
// x(n), the reference signal is x(n-k0), where k0 is the
// specified delay.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (c): Generate 1000 samples of v(n) and x(n), and use the
// NLMS noise canceller to estimate d(n) = sin(0.01PI*n), for
// filter orders p = 5, 10, 15, 20.  Use values for k0 that
// range from the minimum value determined in part (a) to
// k0 = 25.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Generate unit variance white Gaussian noise. random process.
g = generateGaussianProcess(1,1000,1);

// Coefficients for difference equation v(n) = g(n) + 0.5g(n-2).
b = [1 0 0.5];

// Generate MA(2) process.
v = filterBlock(g,b,0);

// Construct desired signal.
d = sin(0.01*%pi*[0:999]);

// Construct input signal.
x = d + v';

[W,dhat] = nlms_noiseCanceller(x,10,1.5,3)

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Plot results.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
