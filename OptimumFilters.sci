//**********************************************************************
// File Name: OptimumFilters.sci
//**********************************************************************

//**********************************************************************
//
//  Name: FirWienerFilter
//
//  Purpose: The purpose of this function is to design a pth-order FIR
//  Wiener filter.
//
//  Calling Sequence: [w,e] = FirWienerFilter(rd,sigmav2,p)
//
//  Inputs:
//
//    rd - A sequence of autocorrelation values for the desired
//    signal, d(n).
//
//    sigmav2 - The variance of the noise signal that perturbs d(n).
//
//    p - The order of the desired Wiener filter.
//
//  Outputs:
//
//    w - The Wiener filter coefficients.
//
//    e - The mean-square error associated with the filter.
//
//**********************************************************************
function [w,e] = FirWienerFilter(rd,sigmav2,p)

  //Force column a vector.
  rd = rd(:);

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Construct the autocorrelation sequence of
  // x(n) = d(n) + v(n).  The result will be
  // rx(k) = rd(k) + sigmav2 * delta(k),
  // where delta(k) is the delta function.
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  rx = rd(1:p);
  rx(1) = rx(1) + sigmav2;
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

  // Construct autocorrelation matrix.
  Rx = toeplitz(rx);

  // Compute Wiener filter coefficients noting that rdx(k) = rd(k)
  w = Rx \ rd(1:p);

  // Compute mean-square error.
  e = rd(1) - w' * rd(1:p);

endfunction

