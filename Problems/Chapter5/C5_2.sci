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

  // Preallocate the matrix.
  A = zeros(p+1,p+1);

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
  // Construct the diagonal and lower trangular portion.
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Iterate through the rows.
  for i = 1:p
    // Iterate through the columns.
    for j = 2:p+1-i;
      aIndex= i + 2 * (j - 1);

      if (aIndex > p+1)
          A(i+j-1,j) = A(i,1);
//printf(">  A(%d,%d): %f\n",i+j-1,j,A(i+j-1,j));
      else
          A(i+j-1,j) = A(i,1) + a(aIndex);
//printf("<= A(%d,%d): %f\n",i+j-1,j,A(i+j-1,j));
      end
    end
  end
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Construct the upper trangular portion.
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Iterate throught the columns.
  for i = 2:p
    // Iterate through the rows.
    for j = 2:p-1
printf("(i,j): (%d,%d)\n",i,j);
      column = i + j - 1;
      if j <> i
printf("In j: (i,column): (%d,%d)\n",i,j);

        A(j,column) = A(i,1);
      end    
    end
  end
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/




//  r = A \ e;

endfunction

//**********************************************************************
// Mainline code.
//**********************************************************************

// Select a simple autocorrelation sequence.
r = [1 .5 .2 .1 .09]';

[a,e] = rtoa(r);

e = [e zeros(1,length(a)-1)]';

[r1,A] = newAtoR(a,e);

