//**********************************************************************
// File Name: SignalModeling.sci
//**********************************************************************

//**********************************************************************
//
//  Name: pade
//
//  Purpose: The purpose of this function is to solve for the
//  coefficients of a filter that models a signal given the number of
//  poles and zeros of the model.
//  The signal is modeled as the unit sample response of a system
//  represented by, H(z) = B(z) / A(z), such that the coefficients of
//  B(z) and A(z) are contained in the vectors, [b(0), b(1),.. b(q)],
//  and [1 a(1), a(2),... a(p)] respectively.
//
//  Calling Sequence: [a,b] = pade(x,p,q)
//
//  Inputs:
//
//    x - The input vector to be processed.
//
//    p - The number of poles in the model.
//
//    q - The number of zeros in the model.
//
//  Outputs:
//
//    a - The denominator coefficients in the model.
//
//    b - The numerator coefficients in the model.
//
//**********************************************************************
function [a,b] = pade(x,p,q)

  // Ensure that we have a column vector.
  x = x(:);

  if (p + q) < length(x)
    X = convm(x,p+1);
    
    // Construct data matrix.
    Xq = X(q+2:q+p+1,2:p+1);

    //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
    // ap= -inv(Xq'*Xq) * Xq' * x1
    //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
    a = [1;-Xq \ X(q+2:q+p+1,1)];

    //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
    // bq = X0 * ap
    //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
    b = X(1:q+1,1:p+1) * a;
  else
    error('Model order too large');
  end

endfunction

//**********************************************************************
//
//  Name: prony
//
//  Purpose: The purpose of this function is to solve for the
//  coefficients of a filter that models a signal given the number of
//  poles and zeros of the model.
//  The signal is modeled as the unit sample response of a system
//  represented by, H(z) = B(z) / A(z), such that the coefficients of
//  B(z) and A(z) are contained in the vectors, [b(0), b(1),.. b(q)],
//  and [1 a(1), a(2),... a(p)] respectively.
//
//  Calling Sequence: [a,b,err] = prony(x,p,q)
//
//  Inputs:
//
//    x - The input vector to be processed.
//
//    p - The number of poles in the model.
//
//    q - The number of zeros in the model.
//
//  Outputs:
//
//    a - The denominator coefficients in the model.
//
//    b - The numerator coefficients in the model.
//
//    err - The mean square error of the approximation to the signal.
//
//**********************************************************************
function [a,b,err] = prony(x,p,q)

  // Ensure that we have a column vector.
  x = x(:);

  N = length(x);

  if (p + q) < length(x) then
    // Create convolution matrix.
    X = convm(x,p+1);

    // Construct the data matrix.
    Xq = X(q+1:N+p-1,1:p);

    //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
    // ap = -inv(Xq'*Xq) * Xq' * xq+1
    //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
    a = [1;-Xq \ X(q+2:N+p,1)];

    //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
    // bq = X0 * ap
    //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
    b = X(1:q+1,1:p+1) * a;

    //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
    // err = Rx*ap
    //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
    err = x(q+2:N)' * X(q+2:N,1:p+1) * a;
    //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  else
    error('Model order too large');
  end

endfunction

//**********************************************************************
//
//  Name: shanks
//
//  Purpose: The purpose of this function is to solve for the
//  coefficients of a filter that models a signal given the number of
//  poles and zeros of the model.
//  The signal is modeled as the unit sample response of a system
//  represented by, H(z) = B(z) / A(z), such that the coefficients of
//  B(z) and A(z) are contained in the vectors, [b(0), b(1),.. b(q)],
//  and [1 a(1), a(2),... a(p)] respectively.
//
//  Calling Sequence: [a,b,err] = shanks(x,p,q)
//
//  Inputs:
//
//    x - The input vector to be processed.
//
//    p - The number of poles in the model.
//
//    q - The number of zeros in the model.
//
//  Outputs:
//
//    a - The denominator coefficients in the model.
//
//    b - The numerator coefficients in the model.
//
//    err - The mean square error of the approximation to the signal.
//
//**********************************************************************
function [a,b,err] = shanks(x,p,q)

  // Ensure that we have a column vector.
  x = x(:);

  N = length(x);

  if (p + q) < length(x)
    // Construct the denominator coefficients.
    a = prony(x,p,q);

    // Construct unit sample vector.
    u = [1; zeros(N-1,1)];

    //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
    // Generate unit sample response. Note that my IIR filter assumes
    // that a monic polynomial is used for A(z), hence, the leading
    // unity term is not needed to be passed to the filter function.
    //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
    g = filterBlock(u,[1],a(2:$));
    //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

    G = convm(g,q+1);

    //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
    // b = inv(G0'*G0)*G0'*x0
    //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
    b = G(1:N,:) \ x;
    
    //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
    // err = rx(0) - rxg * b
    //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
    err = x'*x -x'*G(1:N,1:q+1) * b;
  else
    error('Model order too large');
  end

endfunction

//**********************************************************************
//
//  Name: spike
//
//  Purpose: The purpose of this function is to compute a least squares
//  inverse FIR filter that approximates a delayed unit sample.
//
//  Calling Sequence: [h,err] = spike(g,n0,n)
//
//  Inputs:
//
//    g - The filter for which an inverse it to be computed.
//
//    n0 - The delay of the unit sample to be approximated.
//
//    n - The order of the inverse filter.
//
//  Outputs:
//
//    h - The inverse filter coefficients.
//
//    err - The mean square error of the approximation to the delayed
//    unit sample.
//
//**********************************************************************
function [h,err] = spike(g,n0,n)

  g = g(:);

  m = length(g);

  if (m + n - 1) > n0
    // Construct convolution matrix.
    G = convm(g,n);

    // Construct delayed unit sample.
    d = zeros(m+n-1,1);
    d(n0+1) = 1;

    // Compute coefficients of inverse filter.
    h = G \ d;

    // Compute error.
    err = 1 - G(n0+1,:) * h;
  end

endfunction

//**********************************************************************
//
//  Name: ipf
//
//  Purpose: The purpose of this function is to compute a pole-zero
//  model for a sequence of values that represent a signal using the
//  iterative prefiltering method.
//  The signal is modeled as the unit sample response of a system
//  represented by, H(z) = B(z) / A(z), such that the coefficients of
//  B(z) and A(z) are contained in the vectors, [b(0), b(1),.. b(q)],
//  and [1 a(1), a(2),... a(p)] respectively.
//  
//  Calling Sequence: [a,b,err] = ipf(x,p,q,n,a)
//
//  Inputs:
//
//    x - The input vector to be processed.
//
//    p - The number of poles in the model.
//
//    q - The number of zeros in the model.
//
//    n - The number of iterations to execute.
//
//    a - The initial estimate for the denominator coefficients.
//
//  Outputs:
//
//    a - The denominator coefficients in the model.
//
//    b - The numerator coefficients in the model.
//
//    err - The model error.
//
//**********************************************************************
function [a,b,err] = ipf(x,p,q,n,a)

  // Ensure that we have a column vector.
  x = x(:);

  N = length(x);

  if (p+q) < length(x)
    //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
    // argn() returns is the number of arguments passed to the
    // function.
    // If 5 arguments were not passed to the function, it is
    // implied that the last parameter was not passed.
    //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
    if argn() < 5
      // Fall through to Prony's method for the initial estimate.
      a   = prony(x,p,q);
    end
    //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

    // Construct unit sample function.
    delta = [1; zeros(N-1,1)];

    // Perform n iterations of the estimator.
    for i=1:n
      // f = a(n) * x(n).
      f = filterBlock(x,1,a);

      // g(n) = a(n) * delta(n).
      g = filterBlock(delta,1,a);

      // Construct convolution matrix of f.
      u = convm(f,p+1); 

      // Construct convolution matrix of g.
      v = convm(g,q+1);

      // Compute composite of numerator and denominator coefficients.
      ab = -[u(1:N,2:p+1) -v(1:N,:) ]\u(1:N,1);

      // Extract denominator coefficients.
      a   = [1; ab(1:p)];

      // Extract numerator cofficients.
      b = ab(p+1:p+q+1);

      // Compute error.
      err = norm( u(1:N,1) + [u(1:N,2:p+1) -v(1:N,:)]*ab);
    end
  else
    error('Model order too large');
  end

endfunction

//**********************************************************************
//
//  Name: acm
//
//  Purpose: The purpose of this function is to compute an all pole
//  model of a signal using the autocorrelation method.
//  The signal is modeled as the unit sample response of a system
//  represented by, H(z) = b(0) / A(z), such that the coefficients of
//  A(z) are contained in the vector, [1 a(1), a(2),... a(p)].
//
//  Calling Sequence: [a,err] = acm(x,p)
//
//  Inputs:
//
//    x - The input vector to be processed.
//
//    p - The number of poles in the model.
//
//  Outputs:
//
//    a - The denominator coefficients in the model.
//
//    err - The model error.
//
//**********************************************************************
function [a,err] = acm(x,p)

  // Ensure that we have a column vector.
  x   = x(:);

  N   = length(x);

  if p < length(x)
    // Construct convolution matrix,
    X   = convm(x,p+1);

    // Construct data matrix.
    Xq  = X(1:N+p-1,1:p);

    // Compute denominator coefficients.
    a   = [1;-Xq\X(2:N+p,1)];

    // Compute error.
    err = abs(X(1:N+p,1)'*X*a);
  else
    error('Model order too large');
  end

endfunction

//**********************************************************************
//
//  Name: covm
//
//  Purpose: The purpose of this function is to compute an all pole
//  model of a signal using the covariance method.
//  The signal is modeled as the unit sample response of a system
//  represented by, H(z) = b(0) / A(z), such that the coefficients of
//  A(z) are contained in the vector, [1 a(1), a(2),... a(p)].
//
//  Calling Sequence: [a,err] = covm(x,p)
//
//  Inputs:
//
//    x - The input vector to be processed.
//
//    p - The number of poles in the model.
//
//  Outputs:
//
//    a - The denominator coefficients in the model.
//
//    err - The model error.
//
//**********************************************************************
function [a,err] = covm(x,p)

  // Ensure that we have a column vector.
  x   = x(:);
  N   = length(x);

  if p < length(x)
    // Construct convolution matrix,
    X = convm(x,p+1);

    // Construct data matrix.
    Xq = X(p:N-1,1:p);

    // Compute denominator coefficients.
    a = [1;-Xq\X(p+1:N,1)];

    // Compute error.
    err = abs(X(p+1:N,1)'*X(p+1:N,:)*a);
  else
    error('Model order too large');
  end

endfunction

//**********************************************************************
//
//  Name: durbin
//
//  Purpose: The purpose of this function is to solve for the
//  moving average model parameters for an input sequence using
//  Durbin's method.
//  The signal is modeled as the unit sample response of a system
//  represented by, H(z) = B(z) , such that the coefficients of
//  B(z)  are contained in the vector, [b(0), b(1),.. b(q)].
//
//  Calling Sequence: b = durbin(x,p,q)
//
//  Inputs:
//
//    x - The input vector to be processed.
//
//    p - The order of the all-pole model used to approximate 1/B(z).
//    This parameter should be at least 4*q.
//
//    q - The order of the moving average model.
//
//  Outputs:
//
//
//    b - The parameters of the moving average model.
//
//**********************************************************************
function b = durbin(x,p,q)

  // Ensure that we have a column vector.
  x = x(:);

  if p < length(x)
    //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
    // Use the autocorrelation method to estimate the
    // denominator of an all-pole model for x(n).
    //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
    [a,epsilon] = acm(x,p);

    //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
    // Use the autocorrelation method to estimate the model
    // parameters of an all-pole model for a(n).  This will
    // be an estimate of b(n).
    //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
    [b,epsilon] = acm(length(x)*a/sqrt(epsilon),q);
    b = b*length(x)/sqrt(epsilon);
  else
    error('Model order too large');
  end

endfunction



