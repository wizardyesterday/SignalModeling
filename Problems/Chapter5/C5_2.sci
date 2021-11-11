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
//  coefficients and the modeling error.  Matrix, Ap, is first
//  constructed, and once constructed, the autocorrelation sequence
//  is computed from Ap and the modeling error, e.
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
function r = newAtoR(a,e)

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Compute the length or a.
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  n = length(a);
  p = n - 1;
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

  // Preallocate the matrix, Ap.
  A = zeros(p+1,p+1);

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Construct the first row, first column, and last row.
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  for i = 1:p+1
    // Construct row 1.
    A(1,i) = a(i);

    // Construct column 1.
    A(i,1) = a(i);

    // Construct row p+1.
    A(p+1,p-i+2) = a(i);
  end
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Construct the main diagonal and lower trangular portion.
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Iterate through the rows.
  for i = 1:p
    // Iterate through the columns.
    for j = 2:p+1-i;
      // Compute index into a(n).
      aIndex = i + 2 * (j - 1);

      if (aIndex > p+1)
        A(i+j-1,j) = A(i,1);
      else
        A(i+j-1,j) = A(i,1) + a(aIndex);
     end
    end
  end
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Construct the upper trangular portion.  Note that the
  // loop limits are somewhat naive.  The test for overruns
  // by aIndex keeps data access within bounds.
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Iterate throught the columns.
  for i = 2:p
    // Iterate through the rows.
    for j = 2:p
      // Compute index into a(n).
      aIndex = i + 2 * (j - 1);

      // Ensure that accesses are within bounds.
      if aIndex <= length(a)
        A(j,j+i-1) = a(aIndex);
      end
    end
  end
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

  // Compute autocorrelation sequence.
  r = A \ e;

endfunction

//**********************************************************************
// Mainline code.
//**********************************************************************

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Select a simple autocorrelation sequence.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Model of order 1.
//r = [1 .5]';

// Model of order 2.
//r = [.1 .05 .001]';

// Model of order 3.
//r = [1 .7 .65 .3]';

// Model of order 8.
r = [1 .5 .2 .1 .09 .05 .01 0.009 0.001]';
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

// Compute model estimate from autocorrelation sequence.
[a,e] = rtoa(r);

// Convert into a vector.
e = [e zeros(1,length(a)-1)]';

// Recover autocorrelation sequence.
r1 = newAtoR(a,e);

