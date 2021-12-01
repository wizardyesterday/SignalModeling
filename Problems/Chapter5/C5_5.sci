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
//  The vector sS represents Sp+1(z), and the vector sA represents the
//  vector S*p+1(z), where,
//
//  Sp+1(z) = Ap(z) + [1/z^(p+1)]*Ap(1/z)
//  S*p+1(z) = Ap(z) - [1/z^(p+1)]*Ap(1/z)
//
//  Here are the propertities of Sp+1(z) and S*p+1(z), noting that
//  Sp+1(z) is the symmetric polynomial, and S*p+1(z) is the
//  antisymmetric polynomial (from "Optimal Quantization of LSP Parameters"
//  by Frank K. Soonig and Biing-Hwang Jaung.
//
//  1. All zeros of LSP polynomials are on the unit circle.
//  2. Zeros of Sp+1(z) and S*P+1(z) are interlaced.
//  3. The minimum phase property of Ap(z) are preserved if the first
//  two properties are intact after quantization.
//
//  When p is even, S*p+1(z) has a trivial zero at z = 1, and Sp+1(z)
//  has a trivial zero at z = -1.  When p is odd, there are two trivial
//  zeros, z = 1 and -1, associated with S*p+1(z).
//
//  Note that items 1 and 2 occur if Ap(z) is minimum phase.  I have
//  constructed polynomials that have forced the roots of Sp+1(z) to
//  have equal phases that have magnitudes greater than 1 and less than 1.
//  We are not interested in these cases for application of line spectral
//  pairs though.  The rule to follow here is to make sure you start with
//  a minimum phase polynomial, otherwise, you might be chasing down
//  problems that don't really exist.
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
//    symmetric part of A(z).
//
//    sA - The singular predictor polynomial that represents the
//    antisymmetric part of A(z).
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
//  Note that it is quite simple to compute the frequency differences
//  between roots.  First, the vector of frequencies is sorted, and next,
//  frequency differences are computed.  To restore the original
//  frequencies, it is assummed that the first entry in the output
//  vector has a value of 0.  The remaining entries are computed as an
//  accumulation of their previous entries.  Here is the problem: When
//  the frequencies are sorted, the roots of S(z) and S*(z) are not
//  interleaved as dictated by the ascending order of frequencies.
//  More investigation needs to be carried out before an algorithm that
//  meets the demands of the problem statement can be devised.  For now,
//  I will be happy with *all* the frequencies separated on a per
//  polynomial basis.  Allocation of frequencies to each of the two
//  polynomial roots will be a future exercise.
//
//  Calling Sequence: [freqS,freqA] = lpctolsp(a)
//
//  Inputs:
//
//    a - The linear prediction coefficients.
//
//  Outputs:
//
//  freqS - The frequencies associated with the symmetric singular
//  predictor polynomial.  See the comment block for the atos()
//  function for a detailed description of what freqS represents.
//
//  freqS - The frequencies associated with the antisymmetric singular
//  predictor polynomial.  See the comment block for the atos()
//  function for a detailed description of what freqA represents.
//
//**********************************************************************
function [freqS,freqA] = lpctolsp(a)

  // Force a column vector.
  a = a(:);

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Compute singular predictor polynomials.  Note
  // that sS represents Sp+1(z) and sA represents
  // S*p+1(z), where,
  // Sp+1(z) = Ap(z) + [1/z^(p+1)]*Ap(1/z)
  // S*p+1(z) = Ap(z) - [1/z^(p+1)]*Ap(1/z)
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  [sS,sA] = atos(a);
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

  // Compute roots of singular predictor polynomials.
  rS = roots(sS);
  rA = roots(sA);

  // Compute frequency vectors.
  freqS = atan(imag(rS),real(rS));
  freqA = atan(imag(rA),real(rA));

endfunction

//**********************************************************************
//
//  Name: lsptolpc
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
//  freqA - The frequencies associated with the antisymmetric singular
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
  sS = poly(rS,"z");
  sA = poly(rA,"z");

  // Retrieve the coefficients of the polynomials.
  sS = fliplr(coeff(sS));
  sA = fliplr(coeff(sA));

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
//
//  Name: makeRandomProcess
//
//  Purpose: The purpose of this function is to construct the power
//  spectral density of a random process driving a filter with white
//  noise.
//
//  Calling Sequence: psd = makeRandomProcess(v,a)
//
//  Inputs:
//
//    v - A zero mean white noise sequence.
//
//    a - The filter coefficients.
//
//  Outputs:
//
//    psd - The power spectral density of the random process.
//
//**********************************************************************
function psd = makeRandomProcess(v,a)

// Filter the noise.
y = filterBlock(v,1,a);

// Compute the power spectral densith.
Y = fft(y,-1);
psd = Y .* conj(Y);

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

// System of order 4.
a4 = [1 -1/3 -1/3 2/3 3/4]';
[freq4S,freq4A] = lpctolsp(a4)

// System of order 5.
a5 = [1 -1/3 -1/3 2/3 1/2 1/4]';
[freq5S,freq5A] = lpctolsp(a5)

// System of order 6.
rx6 = [1 0.8 0.5 0.3 0.2 0.1 0.2];
a6 = rtoa(rx6);
[freq6S,freq6A] = lpctolsp(a6)
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
a4hat = lsptolpc(freq4S,freq4A);
a5hat = lsptolpc(freq5S,freq5A);
a6hat = lsptolpc(freq6S,freq6A);
a4_2hat = lsptolpc(freq4_2S,freq4_2A);
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (d): Investigate quantization.  A
// 12-pole filter will be driven by white
// noise in order to generate a random
// process.  Note: I went with a 6th-order
// filter since it was difficult to determine
// a 12-th order filter that would actually
// work.  I may revisit this.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

rx12 = [1 0.8 0.5 0.3 0.2 0.1 0.2];
ac12 = rtoa(rx6);
[freq12S,freq12A] = lpctolsp(ac12);

// Quantize the coefficients to 16 bits.
ac12q = ac12 * 32768;
ac12q = round(ac12q);

// Generate line spectral frequencies.
[freq12S,freq12A] = lpctolsp(ac12);

// Quantize the line spectral frequencies.
freq12Sq = freq12S * 32768;
freq12Sq = round(freq12Sq);
freq12Aq = freq12A * 32768;
freq12Aq = round(freq12Aq);

// Create 12-pole filter from line spectral pairs.
ac12r = lsptolpc(freq12Sq,freq12Aq);

//---------------------------------------------------------------
// Generate 400 samples of Gaussian noise with unity variance.
// Pass segments of the sequence to the 12-pole filter with
// quantized coefficients.  Next, use the line spectral pairs to
// construct its idea of the 12-pole filter. Pass blocks of 100
// samples to each filter so that results can be compared in
// plots of the random process that is created.
//---------------------------------------------------------------
// Run the noise generation function.
noisegen(1,400,1);

// Construct noise source.
v1 = feval([1:100],Noise) / 32768;
v2 = feval([101:200],Noise) / 32768;
v3 = feval([201:300],Noise) / 32768;
v4 = feval([301:400],Noise) / 32768;

// Create ransom processes from original filter.
psd1 = makeRandomProcess(v1,ac12);
psd2 = makeRandomProcess(v2,ac12);
psd3 = makeRandomProcess(v3,ac12);
psd4 = makeRandomProcess(v4,ac12);

// Create random processes from quantized filter.
psd1q = makeRandomProcess(v1/1000,ac12q);
psd2q = makeRandomProcess(v2/1000,ac12q);
psd3q = makeRandomProcess(v3/1000,ac12q);
psd4q = makeRandomProcess(v4/1000,ac12q);

// Create random processes from reconstructed filter.
psd1r = makeRandomProcess(v1,ac12r);
psd2r = makeRandomProcess(v2,ac12r);
psd3r = makeRandomProcess(v3,ac12r);
psd4r = makeRandomProcess(v4,ac12r);

//++++++++++++++++++++++++++++++++++++++++
// Plot the spectral densities.  Don't
// bother plotting the spectral densities
// that resulted from the quantized
// filter. There were overflows that
// resulted in useless data.
//++++++++++++++++++++++++++++++++++++++++
subplot(211);
title('Random Process Spectrum, Original Filter');
plot(psd1);
plot(psd2);
plot(psd3);
plot(psd4);

subplot(212)
title('Random Process Spectrum, Quantized Filter From LSP Frequencies');
plot(psd1r);
plot(psd2r);
plot(psd3r);
plot(psd4r);
//++++++++++++++++++++++++++++++++++++++++

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/






