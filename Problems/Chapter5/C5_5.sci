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
//  Name: lpcToLsp
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
//  Calling Sequence: [freqS,freqA] = lpcToLsp(a)
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
function [freqS,freqA] = lpcToLsp(a)

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
//  Name: lspToLpc
//
//  Purpose: The purpose of this function is to compute the linear
//  prediction coefficients given the line spectral pair frequencies
//  associated with the symmetric and antisymmetric singular predictor
//  polynomials.
//
//  Calling Sequence: a = lspToLpc(freqS,freqA)
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
function a = lspToLpc(freqS,freqA)

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
[freq1S,freq1A] = lpcToLsp(a1)

// System of order 2.
a2 = [1 .5 -.2]';
[freq2S,freq2A] = lpcToLsp(a2)

// System of order 3.
a3 = [1 -1/3 -1/3 2/3]';
[freq3S,freq3A] = lpcToLsp(a3)

// System of order 4.
a4 = [1 -1/3 -1/3 2/3 3/4]';
[freq4S,freq4A] = lpcToLsp(a4)

// System of order 5.
a5 = [1 -1/3 -1/3 2/3 1/2 1/4]';
[freq5S,freq5A] = lpcToLsp(a5)

// System of order 6.
rx6 = [1 0.8 0.5 0.3 0.2 0.1 0.2];
a6 = rtoa(rx6);
[freq6S,freq6A] = lpcToLsp(a6)
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
[freq4_2S,freq4_2A] = lpcToLsp(a4_2)

// Compute frequency response.
[h4_2,fr] = frmag(1,a4_2,1000);
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (c): compute the linear prediction
// filter coefficients given the frequencies
// associated with the singular predictor
// coefficients.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
a1hat = lspToLpc(freq1S,freq1A);
a2hat = lspToLpc(freq2S,freq2A);
a3hat = lspToLpc(freq3S,freq3A);
a4hat = lspToLpc(freq4S,freq4A);
a5hat = lspToLpc(freq5S,freq5A);
a6hat = lspToLpc(freq6S,freq6A);
a4_2hat = lspToLpc(freq4_2S,freq4_2A);
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (d): Investigate quantization.  A
// 12-pole filter will be driven by white
// noise in order to generate a random
// process.  This filter is defined as,
// A12(z) = z^12 + 0.8^12.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

a12 = zeros(1,13);
a12(1) = 1;
a12($) = (0.8)^12;
[freq12S,freq12A] = lpcToLsp(a12);

// Quantize the coefficients to 16 bits.
a12q = a12 * 32768;
a12q = round(a12q);
a12q = a12q / 32768;

// Quantize the line spectral frequencies.
freq12Sq = freq12S * 32768;
freq12Sq = round(freq12Sq);
freq12Sq = freq12Sq / 32768;
freq12Aq = freq12A * 32768;
freq12Aq = round(freq12Aq);
freq12Aq = freq12Aq / 32768;

// Create 12-pole filter from line spectral pairs.
a12r = lspToLpc(freq12Sq,freq12Aq);

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
v1 = feval([1:100],Noise);
v2 = feval([101:200],Noise);
v3 = feval([201:300],Noise);
v4 = feval([301:400],Noise);

// Create ransom processes from original filter.
psd1 = makeRandomProcess(v1,a12);
psd2 = makeRandomProcess(v2,a12);
psd3 = makeRandomProcess(v3,a12);
psd4 = makeRandomProcess(v4,a12);

// Create random processes from quantized filter.
psd1q = makeRandomProcess(v1,a12q);
psd2q = makeRandomProcess(v2,a12q);
psd3q = makeRandomProcess(v3,a12q);
psd4q = makeRandomProcess(v4,a12q);

// Create random processes from reconstructed filter.
psd1r = makeRandomProcess(v1,a12r);
psd2r = makeRandomProcess(v2,a12r);
psd3r = makeRandomProcess(v3,a12r);
psd4r = makeRandomProcess(v4,a12r);

//++++++++++++++++++++++++++++++++++++++++
// Plot the spectral densities.
//++++++++++++++++++++++++++++++++++++++++
subplot(311);
title('Random Process Spectrum, Original Filter');
plot(psd1);
plot(psd2);
plot(psd3);
plot(psd4);

subplot(312);
title('Random Process Spectrum, Quantized Filter');
plot(psd1q);
plot(psd2q);
plot(psd3q);
plot(psd4q);

subplot(313)
title('Random Process Spectrum, Quantized Filter From LSP Frequencies');
plot(psd1r);
plot(psd2r);
plot(psd3r);
plot(psd4r);
//++++++++++++++++++++++++++++++++++++++++

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/






