//**********************************************************************
// File Name: SpectrumEstimation.sci
//**********************************************************************

//**********************************************************************
//
//  Name: zeroPadSequence
//
//  The purpose of this function is to perform adaptive filtering using
//  the Widrow-Hoff LMS algorithm.
//
//  Calling Sequence: [A,E] = lms(x,d,mu,nord,a0)
//
//  Inputs:
//
//    x - The input data to the adaptive filter.
//
//    d - The desired output.
//
//    mu - adaptive filtering update (step-size) parameter.
//
//    nord - The number of filter coefficients.
//
//    a0 - An optional row vector that serves as the initial guess
//    for FIR filter coefficients.  If a0 is omitted, then a0 = 0 is
//    assumed.
//    
//  Outputs:
//
//    xpadded - The zero-padded sequence.
//
//  A - A matrix of filter coefficients as they evolve over time.
//  Each row of this matrix contains the coefficients at the
//  iteration that is associated with the row.  For example, row 1
//  contains the coefficients for iteration 1, row 2 contains the
//  coefficients for iteration 2, and so on.
//
//  E - A vector of errors as they evolve over time.  Each entry
//  contains the error that is associated with each iteration.  For
//  example, the first entry contains the error for iteration 1,
//  the second entry contains the error for iteration 2, and so on.
//  Not that E(n) = d(n) - dhat(n).
//
//**********************************************************************
function [A,E] = lms(x,d,mu,nord,a0)

  // Construct the data matrix.
  X=convm(x,nord);

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
  if if argn(2) < 5
    a0 = zeros(1,N);
  end

  // Force a column vector without altering the values.
  a0   = a0(:).';

  E(1) = d(1) - a0*X(1,:).'; 

  // Perform the first iteration for the filter vector.
  A(1,:) = a0 + mu*E(1)*conj(X(1,:));

  if M > 1
    for k=2:M;
      // Update the error for each iteration.
      E(k) = d(k) - A(k-1,:)*X(k,:).';

      // Update the filter vector for each iteration.
      A(k,:) = A(k-1,:) + mu*E(k)*conj(X(k,:));
    end
  end


endfunction

