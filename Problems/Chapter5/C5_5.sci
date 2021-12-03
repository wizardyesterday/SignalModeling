//**********************************************************************
// File Name: C5_5.sci
//**********************************************************************

// Bring in all of the utility functions.
exec('utils.sci',-1);

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

  // Compute the power spectral density.
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

a12 = zeros(1,13)';
a12(1) = 1;
a12($) = (0.8)^12;
deltaFreq12 = lpcToDeltaLsp(a12);

// Quantize the coefficients to 16 bits.
a12q = round(a12 * 32768) / 32768;

// Quantize the line spectral frequencies.
deltaFreq12q = round(deltaFreq12 * 32768) / 32768;

// Create 12-pole filter from line spectral pairs.
a12r = deltaLspToLpc(deltaFreq12q);

//---------------------------------------------------------------
// Generate 4000 samples of Gaussian noise with unity variance.
// Pass segments of the sequence to the 12-pole filter with
// quantized coefficients.  Next, use the line spectral pair
// differences to construct its idea of the 12-pole filter.
// Pass blocks of 1090 samples to each filter so that results
// can be compared in plots of the random process that is
// created.
//---------------------------------------------------------------
// Run the noise generation function.
noisegen(1,4000,1);

// Construct noise source.
v1 = feval([1:100],Noise);
v2 = feval([1001:2000],Noise);
v3 = feval([2001:3000],Noise);
v4 = feval([3001:4000],Noise);

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
title('Random Process Spectrum, Quantized Filter From LSP Differences');
plot(psd1r);
plot(psd2r);
plot(psd3r);
plot(psd4r);
//++++++++++++++++++++++++++++++++++++++++

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/






