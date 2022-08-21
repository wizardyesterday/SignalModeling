//**********************************************************************
// File Name: SpectrumEstimation.sci
//**********************************************************************

//**********************************************************************
//
//  Name: computeLmsMisadjustment
//
//  The purpose of this function is to compute the theoretical LMS
//  misadjustment given the eigenvalues of the autocorrelation matrix
//  for a random process and the LMS step size parameter.
//
//  Calling Sequence: M = computeLmsMisadjustment(mu,lamda)
//
//  Inputs:
//
//    mu - The adaptive filtering update (step-size) parameter.
//
//    lamda - A vector of eigenvalues of the autocorrelation matrix
//    for a random process.
//    
//  Outputs:
//
//    M - The LMS misadjustment.
//
//**********************************************************************
function M = computeLmsMisadjustment(mu,lamda)

  // Compute the summation portion of the denominator.
  SUM = sum(lamda ./ (2 - mu * lamda));

  // Complete computation of the denominator.
  den = 1 - mu * SUM;

  // Compute misadjustment.
  M = 1 / den;

endfunction

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
  // that the last parameter was not passed.  In this
  // case, the initial condition for the filter coefficients
  // is set to a default value of all zeros. 
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  if argn(2) < 5
    w0 = zeros(1,N);
  end

  // Force a row vector without altering the values.
  w0 = w0(:).';

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
//  Calling Sequence: [W,E] = nlms(x,d,Beta,nord,w0)
//
//  Inputs:
//
//    x - The input data to the adaptive filter.
//
//    d - The desired output.
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
//    E - A vector of errors as they evolve over time.  Each entry
//    contains the error that is associated with each iteration.  For
//    example, the first entry contains the error for iteration 1,
//    the second entry contains the error for iteration 2, and so on.
//    Note that E(n) = d(n) - dhat(n).
//
//**********************************************************************
function [W,E] = nlms(x,d,Beta,nord,w0)

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

  // Perform the first iteration for the error vector.
  E(1) = d(1) - w0 * X(1,:).'; 

  // Construct the normalizing denominator for the first iteration.
  DEN = X(1,:) * X(1,:)' + 0.0001;

  // Perform the first iteration for the filter vector.
  W(1,:) = w0 + Beta / DEN * E(1) * conj(X(1,:));

  if M > 1
    for k = 2:M - nord + 1
      // Update the error.
      E(k) = d(k) - W(k-1,:) * X(k,:).';

      // Update the normalizing denominator.
      DEN = X(k,:) * X(k,:)' + 0.0001;

      // Update the filter vector.
      W(k,:) = W(k-1,:) + Beta / DEN * E(k) * conj(X(k,:));
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
function [W,E] = rls(x,d,nord,Lambda)

  // Set initial reciprocal value for the P matrix.
  delta = 0.001;

  // Construct the data matrix.
  X = convm(x,nord);

  // Retrieve the size of the data matrix.
  [M,N] = size(X);

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // argn(2) returns is the number of arguments passed to the 
  // function.
  // If 3 arguments were passed to the function, it is implied
  // that the last parameter was not passed.  In this
  // case, the initial condition for the filter coefficients
  // is set to a default value of all zeros. 
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  if argn(2) < 4
    // Default to growing window RLS algorithm.
    Lambda = 1.0;
  end

  // Initial value.
  P = eye(N,N) / delta;

  // Initialize first iteration of the filter vector.
  W(1,:) = zeros(1,N);

  for k = 2:M - nord + 1
    // Update the filtered information vector.
    z = P * X(k,:)';

    // Update gain vector.
    g = z / (Lambda + X(k,:) * z);

    // Update the a priori error.
    alpha = d(k) - X(k,:) * W(k-1,:).';

      // Update the filter vector.
    W(k,:) = W(k-1,:) + alpha * g.';

    // Update inverse autocorrelation matrix.
    P = (P - g * z.') / Lambda;

    // Set the returned error value.
    E(k) = alpha;
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
  // that the last parameter was not passed.  In this
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
  // that the last parameter was not passed.  In this
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
  // that the last parameter was not passed.  In this
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

//**********************************************************************
//
//  Name: lms_iirFilteredSignal
//
//  The purpose of this function is to perform adaptive filtering using
//  an IIR filter.  The filtered signal approach to recursive filtering
//  is implemented.
//
//  Calling Sequence: [A,B,e] = lms_iirFilteredSignal(x,d,p,q,mu)
//
//  Inputs:
//
//    x - The input data to the adaptive filter.
//
//    d - The desired output.
//
//    p - The number of denominator coefficients.
//
//    q - The number of numerator coefficients.
//
//  Outputs:
//
//    A - A matrix of filter coefficients as they evolve over time.
//    Each row of this matrix contains the coefficients at the
//    iteration that is associated with the row.  For example, row 1
//    contains the coefficients for iteration 1, row 2 contains the
//    coefficients for iteration 2, and so on.
//
//    B - A matrix of filter coefficients as they evolve over time.
//    Each row of this matrix contains the coefficients at the
//    iteration that is associated with the row.  For example, row 1
//    contains the coefficients for iteration 1, row 2 contains the
//    coefficients for iteration 2, and so on.
//
//    e - A vector of errors as they evolve over time.  Each entry
//    contains the error that is associated with each iteration.  For
//    example, the first entry contains the error for iteration 1,
//    the second entry contains the error for iteration 2, and so on.
//    Note that E(n) = d(n) - dhat(n).
//
//**********************************************************************
function [A,B,e] = lms_iirFilteredSignal(x,d,p,q,mu)

  // Save length of input vector.
  N = length(x);

  // Allocate state memory for IIR filters.
  xState = zeros(1,q+1);
  yState = zeros(1,p);
  xfState = zeros(1,p);
  yfState = zeros(1,p);
  xfInState = 0;
  yfInState = 0;

  // Set initial conditions.
  a = zeros(1,p);
  b = zeros(1,q+1);

  for n = 1:N
    // Compute output value.
    [y(n),xState,yState] = iirFilter(x(n),b,a,xState,yState);

    // Compute error.
    e(n) = d(n) - y(n);

    //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
    // Update filtered signals.  Now, yf(n) is the result of passing
    // y(n) through an all-pole filter, 1/An(z).  Similarly xf(n) is
    // the result of passing x(n) through an all-pole filter, 1/An(z).
    //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
    [yf(n),yfInState,yfState] = iirFilter(y(n),1,a,yfInState,yfState); 
    [xf(n),xfInState,xfState] = iirFilter(x(n),1,a,xfInState,xfState);
    //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

    //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
    // Update the filter coefficients.  Note that E(n) * yf(n-k) and
    // E(n) * xf(n-k) are gradient approximations.
    //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
    for k = 1:p
      if (n-k+1) > 0
        // Update recursive taps.
        a(k) = a(k) + (mu * e(n) * yf(n-k+1));
      end
    end

    for k = 1:q+1
      if (n-k+1) > 0
        // Update nonrecursive taps.
        b(k) = b(k) + (mu * e(n) * xf(n-k+1));
      end

    // Copy coefficients to returned values.
    A(n,:) = a;
    B(n,:) = b;
   //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
    end
  end

endfunction

//**********************************************************************
//
//  Name: nlms_optimized
//
//  The purpose of this function is to perform adaptive filtering using
//  the normalized LMS algorithm.  An optimization has been made for the
//  computation of ||X||^2.  Instead of computing the norm-squared on
//  each iteration, we start with ||X(0)||^2 initially computed.  During
//  loop execution, we compute,
//
//    ||X(n+1)||^2 = ||X(n)||^2 + |x(n+1)|^2 - |x(n-p)|^2
//
//  Now, x(n-p) is the last element of the kth row of the convolution
//  matrix.  What remains to determine is x(n+1).  I used x(k) within
//  the loop, noting that k is used as the loop index rather than
//  n.  Experiments show that the coefficients converge as expected.
//  When k+1, rather than k, was used instability occasionally occurred.
//  It might be that Eq. (9.51) has a typographical error.  This really
//  needs to be investigated.
//
//  Calling Sequence: [W,E] = nlms_optimized(x,d,Beta,nord,w0)
//
//  Inputs:
//
//    x - The input data to the adaptive filter.
//
//    d - The desired output.
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
//    E - A vector of errors as they evolve over time.  Each entry
//    contains the error that is associated with each iteration.  For
//    example, the first entry contains the error for iteration 1,
//    the second entry contains the error for iteration 2, and so on.
//    Note that E(n) = d(n) - dhat(n).
//
//**********************************************************************
function [W,E] = nlms_optimized(x,d,Beta,nord,w0)

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
  w0   = w0(:).';

  // Perform the first iteration for the error vector.
  E(1) = d(1) - w0 * X(1,:).'; 

  // Construct the normalizing denominator for the first iteration.
  norm_X = X(1,:) * X(1,:)';
  DEN = norm_X + 0.0001;

  // Perform the first iteration for the filter vector.
  W(1,:) = w0 + Beta / DEN * E(1) * conj(X(1,:));

  if M > 1
    for k = 2:M - nord + 1
      // Update the error.
      E(k) = d(k) - W(k-1,:) * X(k,:).';

      //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
      // Use the optimization,
      // ||X(n+1)||^2 =
      // ||X(n)||^2 + |x(n+1)|^2 - |x(n-p)|^2.
      //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
      norm_X = norm_X + x(k)^2 - X(k,$)^2;

      // Update the normalizing denominator.
      DEN = norm_X + 0.0001;

      // Update the filter vector.
      W(k,:) = W(k-1,:) + Beta / DEN * E(k) * conj(X(k,:));
    end
  end

endfunction

//**********************************************************************
//
//  Name: lms_pvector
//
//  The purpose of this function is to perform adaptive filtering using
//  the LMS algorithm with a Griffiths p-vector variant.
//
//  Calling Sequence: W = lms_pvector(x,Rdx,mu,nord,w0)
//
//  Inputs:
//
//    x - The input data to the adaptive filter.
//
//    Rdx - The cross-correlation between d(n) and x(n).  This is
//    a matrix of nord columns per row such that a cross-correlation
//    was constructed nord elements at a time.
//    Note: this parameter may have to be reworked.
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
//**********************************************************************
function W = lms_pvector(x,Rdx,mu,nord,w0)

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
  w0   = w0(:).';

  // Perform the first iteration for the filter vector.
  W(1,:) = w0 + mu * Rdx(1,:) - mu * (w0 * X(1,:).') * conj(X(1,:));

  if M > 1
    for k = 2:M - nord + 1

      // Do this to simplify notation.
      correction = mu * Rdx(k,:) - mu * (W(k-1,:) * X(k,:).') * conj(X(k,:));

      // Update the coefficient vector.
      W(k,:) = W(k-1,:) + correction;
    end
  end

endfunction

//**********************************************************************
//
//  Name: lms_variableStepSize
//
//  The purpose of this function is to perform adaptive filtering using
//  the Widrow-Hoff LMS algorithm, but with a variable step size.  The
//  step size is defined as mu(n) = 1 / (c1 + c2*n).
//
//  Calling Sequence: [W,E] = lms_variableStepSize(x,d,c1,c2l,nord,w0)
//
//  Inputs:
//
//    x - The input data to the adaptive filter.
//
//    d - The desired output.
//
//    c1 - A constant that sets the initial convergence rate.
//
//    c2 - A constant that determins the convergence rate as time passes.
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
function [W,E] = lms_variableStepSize(x,d,c1,c2,nord,w0)

  // Construct the data matrix.
  X = convm(x,nord);

  // Retrieve the size of the data matrix.
  [M,N] = size(X);

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // argn(2) returns is the number of arguments passed to the 
  // function.
  // If 5 arguments were passed to the function, it is implied
  // that the last parameter was not passed.  In this
  // case, the initial condition for the filter coefficients
  // is set to a default value of all zeros. 
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  if argn(2) < 6
    w0 = zeros(1,N);
  end

  // Force a row vector without altering the values.
  w0 = w0(:).';

  // Perform the first iteration for the error vector.
  E(1) = d(1) - w0 * X(1,:).'; 

  // Set initial value of the step size.
  mu = 1 / c1;

  // Perform the first iteration for the filter vector.
  W(1,:) = w0 + mu * E(1) * conj(X(1,:));

  if M > 1
    for k = 2:M - nord + 1
      // Update the error.
      E(k) = d(k) - W(k-1,:) * X(k,:).';

      // Set the step size for the current iteration.
      mu = 1 / (c1 + c2*(k - 1));

      // Update the filter vector.
      W(k,:) = W(k-1,:) + mu * E(k) * conj(X(k,:));
    end
  end

endfunction

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
//    dhat - A vector of estimations of the desired signal, d(n).
//    Each entry is associated with iteration n of the adaptive
//    filter.
//
//**********************************************************************
function [W,dhat] = nlms_noiseCanceller(x,n0,Beta,nord,w0)

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Since a reference signal is not provided to this function,
  // it is constructed by delaying the input signal by n0
  // samples.  Advantge is taken of the fact that the noise
  // component of the signal is no longer correlated with the
  // input signal, x(n) delayed by n0 samples.
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

  // Create filter coefficients for the delay line.
  b = zeros(1,n0+1);
  b($) = 1;

  // Construct the reference signal, d(n) = x(n-n0).
  d = filterBlock(x,b,0);
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

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


//*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/
//*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/
//*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/
//
// The following code is still under development. Don't use it yet.
// Here's the status of things.  Formation of the current windows seems
// to be fine.  We slide a window one sample at a time (maybe that's
// correct?).  The previous window is represented by x_(n - L), where
// x_ represents a vector.
// That is, x_(n - L) = [x(n-L) x(n-L-1) ... x(n-l-p]'.  This vector
// needs to be formed, but it's not being done correctly.
//
//*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/
//*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/
//*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/

//**********************************************************************
//
//  Name: rls_slidingWindow
//
//  The purpose of this function is to perform adaptive filtering
//  using the sliding window  recursive least squares algorithm.
//
//  Calling Sequence: W = rls_slidingWindow(x,d,nord,L)
//
//  Inputs:
//
//    x - The input data to the adaptive filter.
//
//    d - The desired output.
//
//    nord - The number of filter coefficients.
//
//    L - The window size in samples.
//    
//  Outputs:
//
//    W - A matrix of filter coefficients as they evolve over time.
//    Each row of this matrix contains the coefficients at the
//    iteration that is associated with the row.  For example, row 1
//    contains the coefficients for iteration 1, row 2 contains the
//    coefficients for iteration 2, and so on.
//
//**********************************************************************
function W = rls_slidingWindow(x,d,nord,L)

  // Set initial reciprocal value for the P matrix.
  delta = 0.001;

  // Construct the data matrix.
  X = convm(x,nord);

  // Retrieve the size of the data matrix.
  [M,N] = size(X);

  // Initial value.
  P = eye(N,N) / delta;

  // Initialize first iteration of the filter vectors.
  W(1,:) = zeros(1,N);
  W_t(1,:) = zeros(1,N);

  // Allocate previous windows.
  xWindowPrev = zeros(1,L);
  dWindowPrev = zeros(1:L);

  // Construct the data matrix.
  X_prev = convm(x,nord);

  // Initialize iteration number.
  n = 2;
 
//  for k = 2:L:M - nord - L
  for k = 2:M - nord - L

    // Set current windows.
    xWindow = x(k:k+L+1);
    dWindow = d(k:k+L);

    // Construct convolution matrices.
    X = convm(xWindow,nord);
    X_prev = convm(xWindowPrev,nord);

    for j = 1:L
      //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
      //++++++++++++++++++++++++++++++++++++++++++++++++++++++
      // First step of the algorithm.  This is the updating
      // step.
      //++++++++++++++++++++++++++++++++++++++++++++++++++++++
      //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

      // Update the filtered information vector.
      z = P * X(j,:)';

      // Update gain vector.
      g = z / (1 + X(j,:) * z);

      // Update the a priori error.
      alpha = dWindow(j) - X(j,:) * W(n-1,:).';

      // Update the filter vector.
      W_t = W(n-1,:) + alpha * g.';

      // Update inverse autocorrelation matrix.
      P_t = P - g * z.';

      //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
      //++++++++++++++++++++++++++++++++++++++++++++++++++++++
      // Second step of the algorithm.  This is the
      // downdating step.
      //++++++++++++++++++++++++++++++++++++++++++++++++++++++
      //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

      // Update the filtered information vector.
      z_t = P_t * X_prev(j,:)';

      // Update gain vector.
      g_t = z_t / (1 - X_prev(j,:) * z_t);

      // Update the a priori error.
      alpha_t = dWindowPrev(j) - X_prev(j,:) * W_t.';

      // Update the filter vector.
      W(n,:) = W_t - alpha_t * g_t.';

      // Update inverse autocorrelation matrix.
      P = P_t + g * z_t.';

      // Increment to the next iteration of w_n.
      n = n + 1;
    end

    // Update the previous buffers.
    xWindowPrev = xWindow(2:$);
    dWindowPrev = dWindow;
  end

endfunction

//*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/
//*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/
//*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/
//
// Remember, don't use the rls_slidingWindow() function, yet.
//
//*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/
//*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/
//*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/

