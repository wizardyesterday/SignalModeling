//**********************************************************************
// File Name: Lsp.sci
//**********************************************************************

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

  // Create extended vector.
  aHat = [a; 0];

  // Construct symmetric singular predictor values.
  sS = aHat + flipud(aHat);

  // Construct antisymmetric singular predictor values.
  sA = aHat - flipud(aHat);

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
      // Working with the frequencies of S*(z).
      if i == 1 | i == n
        // Working with the trivial frequencies of S*(z).
        freqA(j) = f(i);

        // Reference the next entry.
        j = j + 1;
      else
        // Working with the frequency pairs of S*(z).
        freqA(j) = f(i);
        freqA(j+1) = -f(i);

        // Reference the next entry.
        j = j + 2;
      end
    else
      // Working with the frequencies of S(z).
      if i == n
        // Working with the trivial frequency of S(z).
        freqS(k) = f(i);
      else
        // Working with the frequency pairs of S(z).
        freqS(k) = f(i);
        freqS(k+1) = -f(i);

        // Reference the next entry.
        k = k + 2;
      end
    end

    // Toggle to process the other polynomial roots.
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

  // Force column vector.
  a = a(:);

endfunction

