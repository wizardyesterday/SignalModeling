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
//  prediction polynomials given the linear prediction coefficients.
//
//  Calling Sequence: [sS,sA] = atos(a)
//
//  Inputs:
//
//    a - The linear prediction coefficients.
//
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
//  Name: atos
//
//  Purpose: The purpose of this function is to compute the singular
//  prediction polynomials given the linear prediction coefficients.
//
//  Calling Sequence: [sS,sA] = atos(a)
//
//  Inputs:
//
//    a - The linear prediction coefficients.
//
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
function a = lsptolpc(freqS,freqA)

  // Force a column vectors.
  freqS = freqS(:);
  freqA = freqA(:);

  // Compute roots in Cartesian form.
  rS = exp(%i*freqS);
  rA = exp(%i*freqA);

  sS = poly(rS,'z');
  sA = poly(rA,'z');

  sS = coeff(sS);
  sA = coeff(sA);

  sS = fliplr(sS);
  sA = fliplr(sA);

  a = (sS + sA) / 2;
  a($) = [];
  a = real(a);
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
a1 = [1 -.2]';
[sS1,sA1] = atos(a1);
rS1 = roots(sS1);
rA1 = roots(sA1);
aS1 = atan(imag(rS1),real(rS1));
aA1 = atan(imag(rA1),real(rA1));

a1hat = lsptolpc(aS1,aA1);

a2 = [1 .5 -.2]';
[sS2,sA2] = atos(a2);
rS2 = roots(sS2);
rA2 = roots(sA2);
aS2 = atan(imag(rS2),real(rS2));
aA2 = atan(imag(rA2),real(rA2));

a2hat = lsptolpc(aS2,aA2);

a3 = [1 -1/3 -1/3 2/3]';
[sS3,sA3] = atos(a3);
rS3 = roots(sS3);
rA3 = roots(sA3);
aS3 = atan(imag(rS3),real(rS3));
aA3 = atan(imag(rA3),real(rA3));

a3hat = lsptolpc(aS3,aA3);

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/


//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (b): See what happens when a zero of
// sS is close to a zero of sA.
// It seems like the zeros of A(z) become
// closer and closer together for the
// second order case.  The result is a
// strong resonant peak in the frequency
// response of A(z).
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
a4_2 = [1 1.999 1]';
//a4_2 = [1 1 1]';
[sS4,sA4] = atos(a4_2);
rS4 = roots(sS4);
rA4 = roots(sA4);
aS4 = atan(imag(rS4),real(rS4));
aA4 = atan(imag(rA4),real(rA4));

a4_2hat = lsptolpc(aS4,aA4);

// Compute frequency response.
[h4_2,fr] = frmag(1,a4_2,1000);
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/






