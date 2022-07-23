//**********************************************************************
// File Name: C9_1.sci
//**********************************************************************

// Bring in all of the utility functions.
exec('utils.sci',-1);

//**********************************************************************
//
//  Name: lms
//
//  The purpose of this function is to perform adaptive filtering using
//  the Widrow-Hoff LMS algorithm.
//
//  Calling Sequence: [W,E] = lms(x,d,mu,nord,w0)
//
//  Inputs:
//
//    x - The input data to the adaptive filter.
//
//    d - The desired output.
//
//    mu - The adaptive filtering update (step-size) parameter.
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
//    E - A vector of errors as they evolve over time.  Each entry
//    contains the error that is associated with each iteration.  For
//    example, the first entry contains the error for iteration 1,
//    the second entry contains the error for iteration 2, and so on.
//    Note that E(n) = d(n) - dhat(n).
//
//**********************************************************************
function [W,E] = lmsblah(x,d,mu,nord,w0)

  // Construct the data matrix.
  X = convm(x,nord);

  // Retrieve the size of the data matrix.
  [M,N] = size(X);

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // argn(2) returns is the number of arguments passed to the 
  // function.
  // If 4 arguments were passed to the function, it is implied
  // that the last two parameter was not passed.  In this
  // case, the initial condition for the filter coefficients
  // is set to a default value of all zeros. 
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  if argn(2) < 5
    w0 = zeros(1,N);
  end

  // Force a row vector without altering the values.
  w0   = w0(:).';

  // Perform the first iteration for the error vector.
  E(1) = d(1) - w0 * X(1,:).'; 

  // Perform the first iteration for the filter vector.
  W(1,:) = w0 + mu * E(1) * conj(X(1,:));

  if M > 1
    for k = 2:M - nord + 1
      // Update the error.
      E(k) = d(k) - W(k-1,:) * X(k,:).';

      // Update the filter vector.
      W(k,:) = W(k-1,:) + mu * E(k) * conj(X(k,:));
    end
  end

endfunction

//**********************************************************************
// Mainline code.
//**********************************************************************
// Set variances.
sigmavSq = 0.25;

// Initialize white noise generator with unit variance.
noisegen(1,500,sqrt(sigmavSq));

// Set AR(2) filter coefficients.
a1 = [1 -0.1 -0.8];
a2 = [1 -0.1 0.8];
a1 = [1 .1 .8]
a2 = [1 .1 -.8]

// Set step sizes.
mu1 = 0.05
mu2 = 0.01;

//----------------------------------------
// Generate AR(2) processes.
//----------------------------------------
for j = 1:1
  // Generate the noise.
  v = feval([1:500],Noise);

  X1(:,j) = filterBlock(v,1,a1(2:$));
  X2(:,j) = filterBlock(v,1,a2(2:$));

  d1(:,j) = X1(:,j) - v';
  d2(:,j) = X2(:,j) - v';
end
//----------------------------------------

[W1,E1] = lms(X1(:,1),d1,mu1,2)
[W2,E2] = lms(X2(:,1),d2,mu1,2)

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (a): Generate 100 samples of an ARMA(2,2) process, x(n),
// by filtering unit variance white Gaussian noise with a filter
// that has a system function of the form,
// H(z) = {1 + 1.5z^(-1) + z^(-1)} / {1 + 0.9z^(-2).
// Using the method of iterative prefiltering, find the second-
// order ARMA model for x(n).
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/


//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Plot results.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

//scf(1);

