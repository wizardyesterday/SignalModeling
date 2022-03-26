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

  // Force column a vector.
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

//**********************************************************************
//
//  Name: WienerNoiseCanceller
//
//  Purpose: The purpose of this function is to design a pth-order FIR
//  Wiener filter.
//
//  Calling Sequence: [w,e] = WienerNoiseCanceller(x,v1,v2,p)
//
//  Inputs:
//
//    x - An observed signal represented by, x(n) = d(n) + v1(n).
//
//    v1 - The noise, v1(n), that is combined with the observed signal.
//
//    v2 - A signal, v2(n), that is measured by a secondary sensor.
//
//    p - The order of the desired Wiener filter.
//
//  Outputs:
//
//    dhat - The estimate, dhat(n), of the desired signal.
//
//    e - The mean-square error associated with the filter.
//
//**********************************************************************
function [dhat,v1hat,e] = WienerNoiseCanceller(x,v1,v2,p)

  //Force column vectors.
  x = x(:);
  v1 = v1(:);
  v2 = v2(:);

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Compute autocorrelation sequences, and cross-
  // correlation sequences.
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Compute autocorrelation of v1(n).
  rv1 = convol(v1,v1($:-1:1));
  kmax = find(rv1 == max(rv1));
  rv1 = rv1(kmax:kmax+p-1);

  // Compute autocorrelation of v2(n).
  rv2 = convol(v2,v2($:-1:1));
  kmax = find(rv2 == max(rv2));
  rv2 = rv2(kmax:kmax+p-1);

  // Construct cross-correlation of x(n) and v2(n).
  rv1v2 = convol(x,v2($:-1:1));
  kmax = find(rv1v2 == max(rv1v2));
  rv1v2 = rv1v2(kmax:kmax+p-1);
  rv1v2 = rv1v2(:);
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

  // Construct autocorrelation matrix.
  Rv2 = toeplitz(rv2);

  // Compute Wiener filter coefficients.
  w = Rv2 \ rv1v2(1:p);

  // Compute autocorrelation of d(n);
  rd = convol(d,d($:-1:1));
  kmax = find(rd == max(rd));
  rd = rd(kmax:$);

  // Compute estimate of v1(n).
  v1hat = filterBlock(v2,w,0);

  // Compute estimate of d(n).
  dhat = x - v1hat;

  // Compute mean-square error.
  e = rv1(1) - w' * rv1v2(1:p);

endfunction

