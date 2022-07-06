//**********************************************************************
// File Name: SpectrumEstimation.sci
//**********************************************************************

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
function [W,E] = lms(x,d,mu,nord,w0)

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

  // Force a column vector without altering the values.
  w0   = w0(:).';

  // Perform the first iteration for the error vector.
  E(1) = d(1) - w0 * X(1,:).'; 

  // Perform the first iteration for the filter vector.
  W(1,:) = w0 + mu * E(1) * conj(X(1,:));

  if M > 1
    for k = 2:M;
      // Update the error for each iteration.
      E(k) = d(k) - W(k-1,:) * X(k,:).';

      // Update the filter vector for each iteration.
      W(k,:) = W(k-1,:) + mu * E(k) * conj(X(k,:));
    end
  end

endfunction

//**********************************************************************
//
//  Name: lms
//
//  The purpose of this function is to perform adaptive filtering using
//  the normalized LMS algorithm.
//
//  Calling Sequence: [W,E] = lms(x,d,beta,nord,w0)
//
//  Inputs:
//
//    x - The input data to the adaptive filter.
//
//    d - The desired output.
//
//    beta - The adaptive filtering update (normalized step-size)
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
//    E - A vector of errors as they evolve over time.  Each entry
//    contains the error that is associated with each iteration.  For
//    example, the first entry contains the error for iteration 1,
//    the second entry contains the error for iteration 2, and so on.
//    Note that E(n) = d(n) - dhat(n).
//
//**********************************************************************

function [W,E] = nlms(x,d,beta,nord,w0)

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

  // Force a column vector without altering the values.
  w0   = w0(:).';

  // Perform the first iteration for the error vector.
  E(1) = d(1) - w0 * X(1,:).'; 

  // Construct the normalizing denominator for the first iteration.
  DEN = X(1,:) * X(1,:)' + 0.0001;

  // Perform the first iteration for the filter vector.
  W(1,:) = w0 + beta / DEN * E(1) * conj(X(1,:));

  if M>1
    for k = 2:M - nord+1;
      // Update the error for each iteration.
      E(k) = d(k) - W(k-1,:) * X(k,:).';

      // Update the normalizing denominator for each iteration.
      DEN = X(k,:) * X(k,:)' + 0.0001;

      // Update the filter vector for each iteration.
      W(k,:) = W(k-1,:) + beta / DEN * E(k) * conj(X(k,:));
    end
  end

endfunction

