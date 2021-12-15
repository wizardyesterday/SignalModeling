//**********************************************************************
// File Name: Lattice.sci
//**********************************************************************

//**********************************************************************
//
//  Name: fcov
//
//  Purpose: The purpose of this function is to compute an all-pole
//  model for the input sequence, x(n), using the forward covariance
//  method.
//
//  Calling Sequence: [gamn,err] = fcov(x,p)
//
//  Inputs:
//
//    x - A vector of signal values that are to be modeled.
//
//    p - The order of the model.
//
//  Outputs:
//
//    gamm - The vector of reflection coefficients.
//
//   err - The vector of modeling errors.
//
//**********************************************************************
function [gamm,err] = fcov(x,p)

  // Force column vector.
  x = x(:);

  N = length(x);

  // Initialize forward prediction error.
  eplus = x(2:N);

  // Initialize backward prediction error.
  eminus = x(1:N-1);

  N = N - 1;

  for j=1:p;
    // Compute reflection coefficient.
    gamm(j) = -eminus'*eplus/(eminus'*eminus);

    // eplus{j}(n) = eplus{j-1}(n) + gamms{j}*eminus{j-1}(n-1).
    temp1 = eplus  + gamm(j)*eminus;

    // eminus{j}(n) = eminus{j-1}(n-1) + gamma*{j}*eplus{j-1}(n).
    temp2 = eminus + conj(gamm(j))*eplus;

    // Update error.
    err(j) = temp1'*temp1;

    // Update forward prediction error.
    eplus = temp1(2:N);

    // Update backward prediction error.
    eminus = temp2(1:N-1);

    N = N - 1;
  end

endfunction

//**********************************************************************
//
//  Name: bcov
//
//  Purpose: The purpose of this function is to compute an all-pole
//  model for the input sequence, x(n), using the backward covariance
//  method.
//
//  Calling Sequence: [gamn,err] = bcov(x,p)
//
//  Inputs:
//
//    x - A vector of signal values that are to be modeled.
//
//    p - The order of the model.
//
//  Outputs:
//
//    gamm - The vector of reflection coefficients.
//
//   err - The vector of modeling errors.
//
//**********************************************************************
function [gamm,err] = bcov(x,p)

  // Force column vector.
  x = x(:);

  N = length(x);

  // Initialize forward prediction error.
  eplus = x(2:N);

  // Initialize backward prediction error.
  eminus = x(1:N-1);

  N = N - 1;

  for j=1:p;
    // Compute reflection coefficient.
    gamm(j) = -eminus'*eplus/(eplus'*eplus);

    // eplus{j}(n) = eplus{j-1}(n) + gamms{j}*eminus{j-1}(n-1).
    temp1 = eplus  + gamm(j)*eminus;

    // eminus{j}(n) = eminus{j-1}(n-1) + gamma*{j}*eplus{j-1}(n).
    temp2 = eminus + conj(gamm(j))*eplus;

    // Update error.
    err(j) = temp2'*temp2;

    // Update forward prediction error.
    eplus = temp1(2:N);

    // Update backward prediction error.
    eminus = temp2(1:N-1);

    N = N - 1;
  end

endfunction

//**********************************************************************
//
//  Name: burg
//
//  Purpose: The purpose of this function is to compute an all-pole
//  model for the input sequence, x(n).
//
//  Calling Sequence: [gamn,err] = burg(x,p)
//
//  Inputs:
//
//    x - A vector of signal values that are to be modeled.
//
//    p - The order of the model.
//
//  Outputs:
//
//    gamm - The vector of reflection coefficients.
//
//   err - The vector of modeling errors.
//
//**********************************************************************
function [gamm,err] = burg(x,p)

  // Force column vector.
  x = x(:);

  N = length(x);

  // Initialize forward prediction error.
  eplus  = x(2:N);

  // Initialize backward prediction error.
  eminus = x(1:N-1);

  N = N - 1;

  for j = 1:p;
    // Compute reflection coefficient.
    gamm(j) = -2*eminus'*eplus/(eplus'*eplus + eminus'*eminus);

    // eplus{j}(n) = eplus{j-1}(n) + gamms{j}*eminus{j-1}(n-1).
    temp1 = eplus  + gamm(j)*eminus;

    // eminus{j}(n) = eminus{j-1}(n-1) + gamma*{j}*eplus{j-1}(n).
    temp2 = eminus + conj(gamm(j))*eplus;

    // Update error.
    err(j) = temp1'*temp1 + temp2'*temp2;

    // Update forward prediction error.
    eplus = temp1(2:N);

    // Update backward prediction error.
    eminus = temp2(1:N-1);

    N = N - 1;
  end

endfunction

//**********************************************************************
//
//  Name: mcov
//
//  Purpose: The purpose of this function is to compute an all-pole
//  model for the input sequence, x(n), using the modified covariance
//  method.  The model is of the form H(z) = b(0)/A(z).
//
//  Calling Sequence: [a,err] = mcov(x,p)
//
//  Inputs:
//
//    x - A vector of signal values that are to be modeled.
//
//    p - The order of the model.
//
//  Outputs:
//
//    a - The model parameters.
//
//    err - The modeling error.
//
//**********************************************************************
function [a,err] = mcov(x,p)

  // Force column vector.
  x   = x(:);

  N   = length(x);

 if p < length(x)
    // Construct the data matrix.
    X  = toeplitz(x(p+1:N),flipud(x(1:p+1)));

    // Construct the autocorrelation matrix.
    R  = X'*X;

    R1 = R(2:p+1,2:p+1);
    R2 = flipud(fliplr(R(1:p,1:p)));
    b1 = R(2:p+1,1);
    b2 = flipud(R(1:p,p+1));

    // Compute model coeeficients.
    a = [1 ; -(R1 + R2) \ (b1 + b2)];

    // Compute modeling error.
    err = R(1,:)*a + fliplr(R(p+1,:))*a;
  else
    error('Model order too large');
  end

endfunction

//**********************************************************************
//
//  Name: ctob
//
//  Purpose: The purpose of this function is to convert the feedforward
//  polynomial, of a pole-zero lattice filter into the numerator
//  polynomial of H(z) = B(z) / A(z).
//
//  Calling Sequence: b = ctob(c,g)
//
//  Inputs:
//
//    g - The reflection coefficients of the lattice filter.
//
//    c - The coefficients of the feedforward path of the lattice
//    filter.
//
//  Outputs:
//
//    filter.
//    b - The numerator coefficients of H(z).
//
//**********************************************************************
function b = ctob(c,g)

  // Force column vectors.
  c = c(:);
  g = g(:);

  q = length(g);

  // Preallocate storage.
  b = zeros(1,q+1)';
  A = zeros(q+1,q);

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Construct matrix, A, such that,
  // A = [a1 a2 a3 ... ap].
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  for i = 1:q
    a = gtoa([g(1:i)]);
    A(1:length(a),i) = a;
  end
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Add a0 so that A = [a0 a1 a2 a3 ... ap].
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  v = zeros(1:q+1)';
  v(1) = 1;
  A = [v A];
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
 
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Construct b.
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  for k = 1:q+1
    for j = k:q+1
      b(k) = b(k) + c(j)*conj(A(j - k + 1,j));
    end
  end

endfunction
