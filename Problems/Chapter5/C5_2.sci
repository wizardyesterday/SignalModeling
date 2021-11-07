//**********************************************************************
// File Name: C5_2.sci
//**********************************************************************

// Bring in all of the utility functions.
exec('utils.sci',-1);

//**********************************************************************
//
//  Name: newAtoR
//
//  Purpose: The purpose of this function is to provide an alternate
//  way of retrieving an autocorrelation sequence from the model
//  coefficients and the modeling error.
//
//  NOTE: Let's specialize this to a model of order 4 for now.  There
//  seems to be no regular structure that I can see.... yet.
//
//  Calling Sequence: r = newAtoR(a,e)
//
//  Inputs:
//
//    a - The model coefficients.
//
//    e - The modeing error.
//
//  Outputs:
//
//    r - The input autocorrelation sequence.
//
//**********************************************************************
function [r,A] = newAtoR(a,e)

  // Compute the length or a.
  n = length(a);
  p = n - 1;

  // Do this until the function is generalized.
  if p <> 4
    error('Model order must be 4');
  end

  // Preallocate the matrix.
  A = zeros(p+1,p+1);

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Construct the diagonal.
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  A(1,1) = 1;

  for i = 2:p
    A(i,i) = 1 + a(i+1);
  end

  for i = p:p+1
    A(i,i) = 1;
  end
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Construct the first row, first column, and last row.
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  for i = 1:p+1
    A(1,i) = a(i);
    A(i,1) = a(i);
    A(p+1,p-i+2) = a(i);
  end
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Kludge the rest of the entries for a 5x5 matrix.
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  A(2,3) =  a(4);
  A(2,4) = a(5);
  A(3,2) = A(2,1) + a(4);
  A(4,2) = A(3,1) + a(5);
  A(4,3) = a(2);

  A(3,2) = a(2) + a(4);
  A(3,3) = 1 + a(5);
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

  r = A \ e;

endfunction

//**********************************************************************
// Mainline code.
//**********************************************************************

// Select a simple autocorrelation sequence.
r = [1 .5 .2 .1 .09]';

[a,e] = rtoa(r);

e = [e zeros(1,length(a)-1)]';

[r1,A] = newAtoR(a,e);

