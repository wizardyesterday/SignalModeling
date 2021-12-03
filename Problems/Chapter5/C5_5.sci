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
//  Purpose: The purpose of this function is to compute the line
//  spectral pairs given the linear prediction coefficients.
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

  // Force column vectors.
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
//  Name: lpcToDeltaLsp
//
//  Purpose: The purpose of this function is to compute the
//  line spectral pair frequency differences associated with the
//  singular predictor polynomials given the linear prediction
//  coefficients.
//
//  Calling Sequence: lpcToDeltaLsp = lpcToLsp(a)
//
//  Inputs:
//
//    a - The linear prediction coefficients.
//
//  Outputs:
//
//  deltaFreq- The frequencies differences associated with the roots
//  of the symmetric singular predictor polynomial and the
//  antisymmetric singular predictor polynomial.
//
//**********************************************************************
function deltaFreq = lpcToDeltaLsp(a)

  // Force a column vector.
  a = a(:);

  // Compute the line spectral pair.
  [freqS,freqA] = lpcToLsp(a);

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Find and remove all negative frequencies.
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  n = find(freqS < 0);
  freqS(n) = [];

  n = find(freqA < 0)
  freqA(n) = [];
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

  // Combine the spectral frequencies.
  freqSA = [freqS; freqA];

  // Sort the frequencies in ascending order.
  freqSA = flipud(sort(freqSA));

  // Compute the frequency differences.
  deltaFreq = freqSA(2:$) - freqSA(1:$-1);
  
endfunction

//**********************************************************************
//
//  Name: deltaLspToLpc
//
//  Purpose: The purpose of this function is to compute the linear
//  prediction coefficients given the line spectral pair frequency
//  differences associated with the symmetric and antisymmetric
//  singular predictor polynomials.
//
//  Calling Sequence: a = deltaLspToLpc(deltaFreq)
//
//  Inputs:
//
//  deltaFreq- The frequencies differences associated with the roots
//  of the symmetric singular predictor polynomial and the
//  antisymmetric singular predictor polynomial.
//
//  Outputs:
//
//    a - The linear prediction coefficients.
//
//**********************************************************************
function a = deltaLspToLpc(deltaFreq)

  // Force a column vector.
  deltaFreq = deltaFreq(:);

  n = length(deltaFreq) + 1;

  // Preallocate vectors.
  freqS = [0 0]';
  freqA = [0 0]';
  f = [0 0]';

  // Construct line spectral frequencies.
  for i = 2:n
    f(i) = f(i-1) + deltaFreq(i-1);
  end

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Initialize control and index variables for the
  // loop that follows. The index, j, is used for
  // access to the frequencies associated with S*(z),
  // and the index, k, is used for access to the
  // roots associated with S(z).
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Start out filling in the frequencies for S*(z).
  antisymmetric = 1;

  // Initialize indices.
  j = 1;
  k = 1;
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/


  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Fill in all of the frequencies associated with
  // S(z) and S*(z).  The antisymmetric flag, when
  // set to 1, indicates that the frequencies,
  // associated with S*(z) are to be filled in,
  // whereas, when set to -1, the frequencies
  // associated with S(z) are to be filled in.
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  for i = 1:n
    if antisymmetric == 1
      // Working with the roots of S*(z).
      if i == 1 | i == n
        freqA(j) = f(i);
        j = j + 1;
      else
        // Working with the roots of S(z).
        freqA(j) = f(i);
        freqA(j+1) = -f(i);

        // Reference the next entry.
        j = j + 2;
      end
    else
      if i == n
        freqS(k) = f(i);
      else
        freqS(k) = f(i);
        freqS(k+1) = -f(i);

        // Reference the next entry.
        k = k + 2;
      end
    end

    // Toggle to process the other polynomial.
    antisymmetric = -antisymmetric;
  end
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

  // Compute roots from the spectra.
  rS = exp(%i*freqS);
  rA = exp(%i*freqA);

  // Construct S(z) and S*(z) in ascending order of z.
  sS = poly(rS,"z");
  sA = poly(rA,"z");

  // Retrieve coefficients of S(z).
  cS = fliplr(coeff(sS));

  // Retrieve coefficients of S*(z).
  cA = fliplr(coeff(sA));

  // Compute linear predictor coefficients.
  a = (cS + cA)/2;

  // Remove the imaginary part (due to roundoff error).
  a = real(a);

  // Remove extended portion.
  a($) = [];

endfunction

//**********************************************************************
//
//  Name: makeRandomProcess
//
//  Purpose: The purpose of this function is to construct the power
//  spectral density of a random process driving a filter with white
//  noise.  We're really cheating here since the proper way to compute
//  the power spectral density is to 
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

a12 = zeros(1,13);
a12(1) = 1;
a12($) = (0.8)^12;
deltaFreq12 = lpcToDeltaLsp(a12);

// Quantize the coefficients to 16 bits.
a12q = a12 * 32768;
a12q = round(a12q);
a12q = a12q / 32768;

// Quantize the line spectral frequencies.
deltaFreq12q = deltaFreq12 * 32768;
deltaFreq12q = round(deltaFreq12q);
deltaFreq12q = deltaFreq12q / 32768;

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






