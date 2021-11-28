//**********************************************************************
// File Name: C5_5.sci
//**********************************************************************

// Bring in all of the utility functions.
exec('utils.sci',-1);

//**********************************************************************
//
//  Name: atos
//
//  Purpose: The purpose of this function is to compute the singular
//  predictor polynomials given the linear prediction coefficients.
//
//  Calling Sequence: [sS,sA] = atos(a)
//
//  Inputs:
//
//    a - The linear prediction coefficients.
//
//  Outputs:
//
//    sS - The singular predictor polynomial that represents the
//    symmetric part of a(n).
//
//    sA - The singular predictor polynomial that represents the
//    antisymmetric part of a(n).
//
//**********************************************************************
function [sS,sA] = atos(a)

  // Force a column vector.
  a = a(:);

  // Construct symmetric singular predictor values.
  sS = [a; 0] + flipud([a; 0]);

  // Construct antisymmetric singular predictor values.
  sA = [a; 0] - flipud([a; 0]);

endfunction

//**********************************************************************
//
//  Name: lpctolsp
//
//  Purpose: The purpose of this function is to compute the
//  frequencies associated with the singular predictor polynomials
//  given the linear prediction coefficients.
//
//  Calling Sequence: [sS,sA] = lpctolsp(a)
//
//  Inputs:
//
//    a - The linear prediction coefficients.
//
//  Outputs:
//
//  freqS - The frequencies associated with the symmetric singular
//  predictor polynomial.
//
//  freqS - The frequencies associated with the antisymmetric singular
//  predictor polynomial.
//
//**********************************************************************
function [freqS,freqA] = lpctolsp(a)

  // Force a column vector.
  a = a(:);

  // Compute singular predictor polynomials.
  [sS,sA] = atos(a);

  // Compute roots of singular predictor polynomials.
  rS = roots(sS);
  rA = roots(sA);

  // Compute frequency vectors.
  freqS = atan(imag(rS),real(rS));
  freqA = atan(imag(rA),real(rA));

endfunction

//**********************************************************************
//
//  Name: atos
//
//  Purpose: The purpose of this function is to compute the linear
//  prediction coefficients given the line spectral pair frequencies
//  associated with the symmetric and antisymmetric singular predictor
//  polynomials.
//
//  Calling Sequence: a = lsptolpc(freqS,freqA)
//
//  Inputs:
//
//  freqS - The frequencies associated with the symmetric singular
//  predictor polynomial.
//
//  freqS - The frequencies associated with the antisymmetric singular
//  predictor polynomial.
//
//  Outputs:
//
//    a - The linear prediction coefficients.
//
//**********************************************************************
function a = lsptolpc(freqS,freqA)

  // Force a column vectors.
  freqS = freqS(:);
  freqA = freqA(:);

  // Compute roots of the singular predictor polynomials.
  rS = exp(%i*freqS);
  rA = exp(%i*freqA);

  // Construct the singular predictor polynomials from the roots.
  sS = poly(rS,'z');
  sA = poly(rA,'z');

  // Retrieve the coefficients of the polynomials.
  sS = coeff(sS);
  sA = coeff(sA);

  // Compensate for the poly() function produces assending order.
  sS = fliplr(sS);
  sA = fliplr(sA);

  // Compute the linear predictor coefficients.
  a = (sS + sA) / 2;

  // Remove the extraneous element.
  a($) = [];

  // Remove any undesired imaginary part due to roundoff error.
  a = real(a);

  // Force a column vector.
  a = a(:);

endfunction

//**********************************************************************
// Mainline code.
//**********************************************************************
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (a): Generate several different
// singular predictor polynomials, and
// compute the angle of their roots.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// System of order 1.
a1 = [1 -.2]';
[freq1S,freq1A] = lpctolsp(a1)

// System of order 2.
a2 = [1 .5 -.2]';
[freq2S,freq2A] = lpctolsp(a2)

// System of order 3.
a3 = [1 -1/3 -1/3 2/3]';
[freq3S,freq3A] = lpctolsp(a3)

// System of order 5.
a5 = [1 -1/3 -1/3 2/3 1/2]';
[freq5S,freq5A] = lpctolsp(a5)
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (b): See what happens when a zero of
// sS is close to a zero of sA.
// It seems like the zeros of A(z) become
// closer and closer together for the
// second order case.  The result is a
// strong resonant peak in the frequency
// response of 1/A(z).
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
a4_2 = [1 1.999 1]';
//a4_2 = [1 1 1]';
[freq4_2S,freq4_2A] = lpctolsp(a4_2)

// Compute frequency response.
[h4_2,fr] = frmag(1,a4_2,1000);
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (c): compute the linear prediction
// filter coefficients given the frequencies
// associated with the singular predictor
// coefficients.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
a1hat = lsptolpc(freq1S,freq1A);
a2hat = lsptolpc(freq2S,freq2A);
a3hat = lsptolpc(freq3S,freq3A);
a5hat = lsptolpc(freq5S,freq5A);
a4_2hat = lsptolpc(freq4_2S,freq4_2A);
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (d): Investigate quantization.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/






