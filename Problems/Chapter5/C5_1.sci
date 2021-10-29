//**********************************************************************
// File Name: C5+1.sci
//**********************************************************************

// Bring in all of the utility functions.
exec('utils.sci',-1);

//**********************************************************************
//
//  Name: Ldu
//
//  Purpose: The purpose of this function is to perform a Cholesky
//  factorization of a positive definite matrix R.
//
//  Calling Sequence: [A,L,D,U] = Ldu(r)
//
//  Inputs:
//
//    r - The input autocorrelation sequence.
//
//  Outputs:
//
//    A - The model matrix.
//
//    L - The lower triangular matrix for the factorization.
//
//    D - The diagonal matrix that contains the model errors along
//    the main diagonal.
//
//    U - The upper triangular matrix for the factorization.
//
//**********************************************************************
function [A,L,D,U] = Ldu(r)

  // Compute the length or r
  p = length(r);

  // Preallocate the matrix.
  A = zeros(p,p);

  for i = 1:p
    // Compute the (i - 1)th order model.
    [a,e(i)] = rtoa(r(1:i));

    // Add the model to column i, and pad with zeros.
    A(1:i,i) = conj(flipud(a));
  end

  // Compute return values.
  L = inv(A');
  U = L';
  D = diag(e);

endfunction

//**********************************************************************
// Mainline code.
//**********************************************************************

// Select a simple autocorrelation sequence.
r = [1 .5 .2];

// Create the Toeplitz matrix from r(n).
R = toeplitz(r);

// Perform the factorization.
[A,L,D,U] = Ldu(r);

