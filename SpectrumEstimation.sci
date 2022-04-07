//**********************************************************************
// File Name: SpectrumEstimation.sci
//**********************************************************************
//**********************************************************************
//
//  Name: per
//
///  The purpose of this function is to compute the spectrum of a
//   process using the periodogram.

//  Calling Sequence: Px = per(x,n1,n2)
//
//  Inputs:
//
//    x - The input sequence.
//
//    n1 - The starting index such that processing starts x(n1)
//
//    n2 - The ending index, such that processing ends at x(n2).  If
//    n1 and n2 are not specified, the whole input sequence is
//    processed.
//
//  Outputs:
//
//    Px - The periodogram estimate of the power spectrum of x(n).
//
//**********************************************************************
function Px = per(x,n1,n2)

  // Force a column vector.

  x = x(:);

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // argn(2) returns is the number of arguments passed to the
  // function.
  // If 1 argument was passed to the function, it is implied
  // that the last two parameters were not passed.  In this
  // case, the whole input sequence is processed.
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  if argn(2) == 1
    n1 = 1;
    n2 = length(x);
  end
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Compute the periodogram.
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Compute Fourier transform.
  X = fft(x(n1:n2),-1);

  // Compute power spectral density.
  Px = (X .* conj(X)) / (n2 - n1 + 1) ^2;

  Px(1)=Px(2);
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

endfunction






