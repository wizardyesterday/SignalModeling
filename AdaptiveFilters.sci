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
//
//  Name: nlms
//
//  The purpose of this function is to perform adaptive filtering using
//  the normalized LMS algorithm.
//
//  Calling Sequence: [W,E] = nlms(x,d,beta,nord,w0)
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

  // Force a row vector without altering the values.
  w0   = w0(:).';

  // Perform the first iteration for the error vector.
  E(1) = d(1) - w0 * X(1,:).'; 

  // Construct the normalizing denominator for the first iteration.
  DEN = X(1,:) * X(1,:)' + 0.0001;

  // Perform the first iteration for the filter vector.
  W(1,:) = w0 + beta / DEN * E(1) * conj(X(1,:));

  if M > 1
    for k = 2:M - nord + 1
      // Update the error.
      E(k) = d(k) - W(k-1,:) * X(k,:).';

      // Update the normalizing denominator.
      DEN = X(k,:) * X(k,:)' + 0.0001;

      // Update the filter vector.
      W(k,:) = W(k-1,:) + beta / DEN * E(k) * conj(X(k,:));
    end
  end

endfunction

//**********************************************************************
//
//  Name: rls
//
//  The purpose of this function is to perform adaptive filtering
//  using the exponentially weighted recursive least squares algorithm.
//  If the parameter, lamda has a value of 1, the growing window RLS
//  will result.
//
//  Calling Sequence: [W,E] = rls(x,d,nord,lambda)
//
//  Inputs:
//
//    x - The input data to the adaptive filter.
//
//    d - The desired output.
//
//    nord - The number of filter coefficients.
//
//    lamda - The exponential forgetting factor.  A value of 1
//    results in the growing window RLS algorithm being executed.
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
function [W,E] = rls(x,d,nord,lambda)

  // Set initial reciprocal value for the P matrix.
  delta = 0.001;

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
  if argn(2) < 4
    // Default to growing window RLS algorithm.
    lambda = 1.0;
  end

  // Initial value.
  P = eye(N,N) / delta;

  // Initialize first iteration of the filter vector.
  W(1,:) = zeros(1,N);

  for k = 2:M - nord + 1
    // Update the filtered information vector.
    z = P * X(k,:)';

    // Update gain vector.
    g = z / (lambda + X(k,:) * z);

    // Update the apriori error.
    alpha = d(k) - X(k,:) * W(k-1,:).';

      // Update the filter vector.
    W(k,:) = W(k-1,:) + alpha * g.';

    // Update inverse autocorrelation matrix.
    P = (P - g * z.') / lambda;
  end

endfunction

//**********************************************************************
//
//  Name: lms_signError
//
//  The purpose of this function is to perform adaptive filtering using
//  the Widrow-Hoff LMS algorithm. A sign-error variant has been
//  incorporated.
//
//  Calling Sequence: [W,E] = lms_signEerror(x,d,mu,nord,w0)
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
function [W,E] = lms_signError(x,d,mu,nord,w0)

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

      // Update the filter vector using the sign-error variant.
      W(k,:) = W(k-1,:) + mu * sign(E(k)) * X(k,:);
    end
  end

endfunction

//**********************************************************************
//
//  Name: lms_signData
//
//  The purpose of this function is to perform adaptive filtering using
//  the Widrow-Hoff LMS algorithm.  A sign-data variant has been
//  implemented.
//
//  Calling Sequence: [W,E] = lms_signData(x,d,mu,nord,w0)
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
function [W,E] = lms_signData(x,d,mu,nord,w0)

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

      // Update the filter vector using the sign-data variant.
      W(k,:) = W(k-1,:) + mu * E(k) * sign(X(k,:));
    end
  end

endfunction

//**********************************************************************
//
//  Name: lms_signSign
//
//  The purpose of this function is to perform adaptive filtering using
//  the Widrow-Hoff LMS algorithm.  A sign-sign variant has been
//  implemented.
//
//  Calling Sequence: [W,E] = lms_signSign(x,d,mu,nord,w0)
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
function [W,E] = lms_signSign(x,d,mu,nord,w0)

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
  W(1,:) = w0 + mu * sign(E(1)) * sign(X(1,:));

  if M > 1
    for k = 2:M - nord + 1
      // Update the error.
      E(k) = d(k) - W(k-1,:) * X(k,:).';

      // Update the filter vector using the sign-sign variant.
      W(k,:) = W(k-1,:) + mu * sign(E(k)) * sign(X(k,:));
    end
  end

endfunction

