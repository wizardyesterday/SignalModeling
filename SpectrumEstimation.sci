//**********************************************************************
// File Name: SpectrumEstimation.sci
//**********************************************************************

//**********************************************************************
//
//  Name: zeroPadSequence
//
//  The purpose of this function is to zero-pad a sequence by appending
//  zeros to the sequence. Note that if the input sequence is longer
//  than the desired length of the padded sequence, no processing will
//  occur.  This function can be used for the case that it is desired
//  to zero-pad a sequence to a length that is usable for input to an
//  FFT.  For example, the input sequence might have a length of 124
//  entries, and it might be desired to compute a 1024-point FFT.  Now,
//  Scilab is nice enough to round up to the next power of two for
//  computation of the FFT, but the user really has no control of
//  specifying the size of the FFT.  This function circumvents that
//  problem.
//
//  Calling Sequence: xpadded = zeroPadSequence(x,N,n1,n2)
//
//  Inputs:
//
//    x - The input sequence.
//
//    N - The desired padded sdquence length.
//
//  Outputs:
//
//    xpadded - The zero-padded sequence.
//
//**********************************************************************
function xpadded = zeroPadSequence(x,N)

  // Force a column vector.
  x = x(:);

  // Set to a default value.
  xpadded = x;

  // Compute the desired vector length of x(n).  
  vectorLength = length(x);

  if vectorLength < N
    // Compute the number of zeros to append to x(n).
    paddingLength = N - vectorLength;

    // Create the zeross
    zeroPadVector = zeros(paddingLength,1);

    // Create zero-padded vector.
    xpadded = [x; zeroPadVector];
  end

endfunction

//**********************************************************************
//
//  Name: truncatedAutocorrelation
//
//  The purpose of this function is to compute a truncated
//  autocorrelation sequence that such that for rx(k), -M <= k <= M.
//
//  Calling Sequence: r = truncatedAutocorrelation(x,M,n1,n2)
//
//  Inputs:
//
//    x - The input sequence.
//
//    M - The bounds of the lag window, -M <= k <= M.
//
//    n1 - The starting index such that processing starts x(n1)
//
//    n2 - The ending index, such that processing ends at x(n2).  If
//    n1 and n2 are not specified, the whole input sequence is
//    processed.
//
//  Outputs:
//
//    r - The autocorrelation sequence of length 2M +1 samples.
//
//**********************************************************************
function r = truncatedAutocorrelation(x,M,n1,n2)

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // argn(2) returns is the number of arguments passed to the
  // function.
  // If 1 argument was passed to the function, it is implied
  // that the last two parameters were not passed.  In this
  // case, the whole input sequence is processed.
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  if argn(2) == 2
    n1 = 1;
    n2 = length(x);
  end
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

  // Compute  autocorrelation sequence for positive lags.
  r = convol(x(n1:n2),conj(x(n2:-1:1)));
  kmax = find(abs(r) == max(abs(r)));
  r = r(kmax:$);

  // Truncate the autocorrelation sequence to M lags.
  r = r(1:M);

  // Make the autocorrelation sequence double-sided.
  r = [conj(r(M:-1:2)) r];

  // Force column vector.
  r = r(:);

endfunction

//**********************************************************************
//
//  Name: autocorrelationMatrix
//
//  The purpose of this function is to create an autocorrelation
//  matrix, given an autocorrelation sequence.  Given an input
//  autocorrelation sequence, a Toeplitz matrix will first be
//  constructed.  Next, the complex conjugate of the appropriate
//  entries of the matrix is taken.  Here's an example.
//
//    rx = [1+i1 2 3-i7]
//
//    Rx = [[1+i1 2 3+i7]
//          [2 1+i1] 2]
//          [3-i7 2 1+i]]
//
//  The result, R, is a Hermitian Toeplitz matrix.
//
//  Calling Sequence: R = autocorrelationMatrix(r)
//
//  Inputs:
//
//    r - The autocorrelation sequence.
//
//  Outputs:
//
//    R - The autocorrelation matrix.
//
//**********************************************************************
function R = autocorrelationMatrix(r)

  // Compute order of the autocorrelation sequence.
  p = length(r);

  // Construct the initial pxp matrix.
  R = toeplitz(r);

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Take the complex conjugate of the entries to
  // account for the conjugate symmetry of the
  // autocorrelation function.
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  for row = 1:p
    for column = row+1:p
      R(row,column) = conj(R(row,column));
    end
  end
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

endfunction

//**********************************************************************
//
//  Name: Covar
//
//  The purpose of this function is to compute an autocorrelation
//  matrix of order, p, from an input sequence, x(n).  The matrix is
//  created by the following two steps:
//
//    1. Construct a single-sided autocorrelation sequence of order, p.
//
//    2. Use the autocorrelationMatrix() function to construct the
//    autocorrelation matrix.
//
//  Calling Sequence: R = Covar(x,p)
//
//  Inputs:
//
//    x - The input sequence.
//
//    p - The order of the autocorrelation matrix.
//
//  Outputs:
//
//    R - The p x p autocorrelation matrix.
//
//**********************************************************************
function R = Covar(x,p)

  if p <= length(x)
    // Compute  autocorrelation sequence for positive lags.
    r = convol(x,conj(x($:-1:1)));
    kmax = find(abs(r) == max(abs(r)));
    r = r(kmax:$);

    // Truncate the autocorrelation sequence to p lags.
    r = r(1:p);

    // Create the autocorrelation matrix.
    R = autocorrelationMatrix(r);
  else
    error('Order too large');
  end

endfunction

//**********************************************************************
//
//  Name: matrixFft
//
//  The purpose of this function is to compute the FFT of a matrix
//  that is composed of column vectors.  Each column of the output
//  matrix will be the FFT of each column of the input matrix.  The
//  goal of this function is to strongarm the Scilab FFT of a matrix
//  to perform like the default behavior of a MATLAB FFT.  Scilab has
//  some interesting indexing with respect to FFT processing.
//
//  Calling Sequence: V =  matrixFft(v,N)
//
//  Inputs:
//
//    v - The time-domain input matrix.  Each column of this matrix
//    represents a separate time-domain sequence.
//
//    N - The desired Fast Fourier transform length for each column.
//
//  Outputs:
//
//    V - The frequency-domain output matrix. Each column of this
//    matrix represents the frequency-domain representation of each
//    column of the input matrix.
//
//**********************************************************************
function V = matrixFft(v,N)

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // We need to compute a N-point FFT, so all
  // of this magic to zero-pad each column of the
  // matrix so that they are at least of length N.
  // A matrix of zero-padded columns of will be
  // first created.  Next, the funky indexing stuff
  // will be used to compute the FFT of each
  // column of the matrix.
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Compute the number of columns of v.
  rowsColumns = size(v);
  columns = rowsColumns(2);

  // Zero-pad the first column of v.
  v1 = zeroPadSequence(v(:,1),N);

  // Zero-pad the rest of the columns in v1.
  for i = 2:columns
    v1(:,i) = zeroPadSequence(v(:,i),N);
  end

  // Compute the N-pointspectrum for each column.
  V = fft(v1,-1,N,1);
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

endfunction

//**********************************************************************
//
//  Name:  computeDtft
//
//  The purpose of this function is to compute the discrete Fourier
//  transform of a short sequence.
//
//  Calling Sequence: V = computeDtft(v,p,w)
//
//  Inputs:
//
//    v - The input sequence.  This is a (p+1) x 1 column vector.
//
//    p - The order of the DFT.
//
//    w - The radian frequency for which the DTFT is to be evaluated.
//
//  Outputs:
//
//    V - The discrete time Fourier transform of v.
//
//**********************************************************************
function V = computeDtft(v,p,w)

  // Clear sum.
  V = 0;

  // Compute the DTFT of v.
  for k = 0:p
    V = V + v(k+1) * exp(-%i*k*w);
  end

endfunction

//**********************************************************************
//
//  Name:  computeEigenPowers
//
//  The purpose of this function is to estimate the powers of complex
//  exponentials in a white noise process given the eigenvectors of
//  the autocorrelation matrix and their associated eigenvalues, the
//  variance of the random process, and the frequency estimates of
//  the complex exponentials.
//
//  Calling Sequence: P = computeEigenPowers(v,lamda,sigma2,omega)
//
//  Inputs:
//
//    v - A p x (p+1) matrix of column eigenvectors.
//
//    lamda - A (p+1) x 1 column vector of eigenvalues that are 
//    associated with V.
//
//    sigma2 - The noise power.
//
//    omega - A p x 1 column vector of frequency estimates.
//
//  Outputs:
//
//    P - The power of each sinusoid.
//
//**********************************************************************
function P = computeEigenPowers(v,lamda,sigma2,omega)

  // Determine size of V.
  s = size(v);
  p = s(2);

  for i = 1:p
    for k = 1:p
      // Compute V_i(omega_k).
      V(i,k) = computeDtft(v(:,i),p,omega(k));

      // Compute |V_i(omega_k)|^2.
      V(i,k) = V(i,k) .* conj(V(i,k));
    end
  end

  // Subtract out the noise term.
  lamda = lamda - sigma2;

  // Compute power estimates.
  P = V \ lamda;

  // Remove infinitesimal imaginary component.
  P = real(P);

endfunction

//**********************************************************************
//
//  Name: constructPowerSpectrum
//
//  Purpose: The purpose of this function is to generate a 1024-point
//  power spectrum representation.
//
//  Calling Sequence: Px = constructPowerSpectrum(x)
//
//  Inputs:
//
//    x - The input vector to be processed.

//
//  Outputs:
//
//  Px - The power spectrum in linear units.
//
//**********************************************************************
function Px = constructPowerSpectrum(x)

  // Force column vector.
  x = x(:);

  // Construct 1024-point FFT.
  X = matrixFft(x,1024);

  // Construct power spectrum.
  Px = X .* conj(X);

endfunction

//**********************************************************************
//
//  Name: blackman
//
//  The purpose of this function is to compute the Blackman window. 
//
//  Calling Sequence: w = blackman(n)
//
//  Inputs:
//
//    n - The number of samples of the Blackman window.
//
//  Outputs:
//
//    w = The sequence that represents the Blackman window function.
//
//**********************************************************************
function w = blackman(n)

  // Setup time indices.
  k = 0:n-1;

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Set the coefficients for the Blackman window.

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  a0 = 0.42;
  a1 = -0.50;
  a2 = 0.08;
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

  // Compute the window function.
  w = a0 + a1*cos(2*%pi*k / (n - 1)) + a2*cos(4*%pi*k / (n - 1));

endfunction

//**********************************************************************
//
//  Name: per
//
//  The purpose of this function is to estimate the spectrum of a
//   process using the periodogram.
//
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
//    Px - The periodogram estimate of the power spectrum of x(n)
//    using a linear scale.
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
  // Compute the periodogram.  A 1024-point FFT will be used
  // for the computation.
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Zero-pad if necessary.
  x = zeroPadSequence(x(n1:n2),1024);

  // Compute the FFT of the zero-padded sequence.
  X = fft(x,-1);

  // Compute power spectral density.
  Px = (X .* conj(X)) / (n2 - n1 + 1) ^2;

  Px(1) = Px(2);
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

endfunction

//**********************************************************************
//
//  Name: mper
//
//  The purpose of this function is to estimate the spectrum of a
//   process using the modified periodogram.
//
//  Calling Sequence: Px = mper(x,win,n1,n2)
//
//  Inputs:
//
//    x - The input sequence.
//
//    win - The window type.  Valid values are 1 (Rectangular),
//    2 (Hamming), 3 (Hanning), 4 (Bartlett), and 5 (Blackman).
//
//    n1 - The starting index such that processing starts x(n1)
//
//    n2 - The ending index, such that processing ends at x(n2).  If
//    n1 and n2 are not specified, the whole input sequence is
//    processed.
//
//  Outputs:
//
//    Px - The modified periodogram estimate of the power spectrum
//    of x(n) using a linear scale.
//
//**********************************************************************
function Px = mper(x,win,n1,n2)

  // Force a column vector.
  x = x(:);

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // argn(2) returns is the number of arguments passed to the
  // function.
  // If 2 arguments were  passed to the function, it is implied
  // that the last two parameters were not passed.  In this
  // case, the whole input sequence is processed.
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  if argn(2) == 2
    n1 = 1;
    n2 = length(x);
  end
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

  // Compute sequence length of interest.
  N  = n2 - n1 +1;

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Determine the desired window type.
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  select win

    case 1
      // Select Rectangular window.
      w  = window('re',N)';

    case 2
      // Select Hamming window.
      w = window('hm',N)';

    case 3
      // Select Hanning window.
      w = window('hn',N)';

    case 4
      // Select Bartlett window.
      w = window('tr',N)';

    case 5
      // Select Blackman window.
      w = blackman(N)';

  end // select
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

  // Compute windowed sequence.
  xw = x(n1:n2) .* w' / norm(w);

  // Compute the periodogram.
  Px = N * per(xw);

endfunction

//**********************************************************************
//
//  Name:  bart
//
//  The purpose of this function is to estimate the spectrum of a
//   process using Bartlett's method of periodogram averaging.
//
//  Calling Sequence: Px = bart(x,nsect)
//
//  Inputs:
//
//    x - The input sequence.
//
//    nsect - The number of subsequences to be used in the average.
//
//  Outputs:
//
//    Px - The Bartlett estimate of the power spectrum of x(n) using
//    a linear scale.
//
//**********************************************************************
function Px = bart(x,nsect)

  // Compute length of each subsection.
  L = floor(length(x) / nsect);

  // Clear average.
  Px = 0;

  // Start with the first subsequence.
  n1 = 1;

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Compute the average of the periodograms of each
  // subsequence. Note that division by the number of
  // subsequences occurs within the loop.
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  for i = 1:nsect
    // Add the next periodogram to the average.
    Px = Px + per(x(n1:n1+ L - 1)) / nsect;

    // Increment to the next subsection.
    n1 = n1 + L;
  end
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

endfunction

//**********************************************************************
//
//  Name:  welch
//
//  The purpose of this function is to estimate the spectrum of a
//   process using Welch's method of modified periodogram averaging.
//
//  Calling Sequence: Px = welch(x,L,over,win)
//
//  Inputs:
//
//    x - The input sequence.
//
//    over - The amount of overlap, where 0 < over < 1.
//
//    L - The section length.
//
//    win - The window type.  Valid values are 1 (Rectangular),
//    2 (Hamming), 3 (Hanning), 4 (Bartlett), and 5 (Blackman).
//
//  Outputs:
//
//    Px - The Welch estimate of the power spectrum of x(n) using
//    a linear scale.
//
//**********************************************************************
function Px = welch(x,L,over,win)

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // argn(2) returns is the number of arguments passed to the
  // function.
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  if argn(2) <= 3
    win = 1;
  end

  if argn(2) <= 2
   over = 0;
  end

  if argn(2) == 1
   L = length(x);
  end
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

  if over > 0 & over < 1
    // Set starting index for the first subsequence.
    n1 = 1;

    // Set the increment for the next subsequence.
    n0 = (1 - over) * L;

    // Compute the number of subsequences.
    nsect = 1+ floor((length(x) -L) / n0);

    // Clear average.
    Px = 0;

    for i = 1:nsect
      // Add the next periodogram to the average.
      Px = Px + mper(x,win,n1,n1 + L - 1) / nsect;

      // Reference the next subsequence.
      n1 = n1 + n0;  
    end
  else
    error('Overlap is invalid');
  end

endfunction

//**********************************************************************
//
//  Name: per_smooth
//
//  The purpose of this function is to estimate the spectrum of a
//   process using Blackman-Tukey method.
//
//  Calling Sequence: Px = per_smooth(x,win,M,n1,n2)
//
//  Inputs:
//
//    x - The input sequence.
//
//    win - The window type.  Valid values are 1 (Rectangular),
//    2 (Hamming), 3 (Hanning), 4 (Bartlett), and 5 (Blackman).
//
//    M - The bound of the lag window, such that -M <= k <= M for
//    rx(k) and w(k) (the lag window).
//
//    n1 - The starting index such that processing starts x(n1)
//
//    n2 - The ending index, such that processing ends at x(n2).  If
//    n1 and n2 are not specified, the whole input sequence is
//    processed.
//
//  Outputs:
//
//    Px - The Blackman-Tukey estimate of the power spectrum
//    of x(n) using a linear scale.
//
//**********************************************************************
function Px = per_smooth(x,win,M,n1,n2)

  // Force a column vector.
  x = x(:);

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // argn(2) returns is the number of arguments passed to the
  // function.
  // If only three arguments were supplied to this function,
  // this implies that n1 and n2 were not specified.
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  if argn(2) == 3
    n1 = 1;
    n2 = length(x);
  end
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

  // Create truncated double-sided autocorrelation sequence.
  r = truncatedAutocorrelation(x,M,n1,n2)

  // Adjust for the proper length of the lag window.
  M = 2 * M - 1;

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Determine the desired window type.
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  select win

    case 1
      // Select Rectangular window.
      w  = window('re',M)';

    case 2
      // Select Hamming window.
      w = window('hm',M)';

    case 3
      // Select Hanning window.
      w = window('hn',M)';

    case 4
      // Select Bartlett window.
      w = window('tr',M)';

    case 5
      // Select Blackman window.
      w = blackman(M)';

  end // select
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

  // Multiply the autocorrelation by the window.
  r = r .* w;

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Compute estimated spectrum.
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Zero-pad to the desired FFT size.
  r = zeroPadSequence(r,1024);

  // Compute the spectrum.
  Px = abs(fft(r,-1));

  Px(1) = Px(2);
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

endfunction

//**********************************************************************
//
//  Name:  minvar
//
//  The purpose of this function is to estimate the spectrum of a
//  process using the minimum variance method.
//
//  Calling Sequence: Px = minvar(x,p)
//
//  Inputs:
//
//    x - The input sequence.
//
//    p - The order of the minimum variance estimate.  For short
//    sequences, p is typically about length(x) / 3.
//
//  Outputs:
//
//    Px - The minimum variance estimate of the power spectrum of
//    x(n) using a decibel scale.
//
//**********************************************************************
function Px = minvar(x,p)

  // Enforce a column vector.
  x = x(:);

  // Compute the autocorrelation matrix of order, p.
  R = Covar(x,p);

  // Compute the eigenvectors of R and the diagonal matrix of eigenvalues.
  [v,d] = spec(R);

  // Construct a column vector of eigenvalue reciprocals.
  U = diag(inv(abs(d)+ %eps));

  // Compute the 1024-pointspectrum for each column.
  V = matrixFft(v,1024);
  V = V .* conj(V);

  // Estimate the power spectrum in decibels.
  Px = (10 * log10(p)) - (10 * log10(V * U));

endfunction

//**********************************************************************
//
//  Name:  mem
//
//  The purpose of this function is to estimate the spectrum of a
//  process using the maximum entropy method.  The autocorrelation
//  method is used to find a pth-order all-pole model for x(n), and
//  the spectral estimate is formed from the following equation:
//
//    Px = b^2(0) / |A(omega)|^2 
//
//  Calling Sequence: Px = mem(x,p)
//
//  Inputs:
//
//    x - The input sequence.
//
//    p - The order of the all-pole model.
//
//  Outputs:
//
//    Px - The maximum entropy estimate of the power spectrum of
//    x(n) using a decibel scale.
//
//**********************************************************************
function Px = mem(x,p)

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Construct an all-pole model using the autocorrelation 
  // matrix method.
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  [a,e] = acm(x,p);
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

  // Compute the FFT of the all-pole model.
  A = matrixFft(a,1024);

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Compute the estimated power spectral density.  Note that
  // e = b^(0), and that e is being normalized by the length
  // of the input vector.  Also, notice that the second log()
  // term is being multiplied by 20 to account that A
  // represents a voltage rather than a power (|A|^2),
  // Had we been working with |A|^2, the multiplier would have
  // been 10. 
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  Px = 10 * log10(e / length(x)) - 20 * log10(abs(A));
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

endfunction

//**********************************************************************
//
//  Name:  phd
//
//  The purpose of this function is to estimate the frequencies
//  and powers of a sum of complex sinusids in white noise using
//  Pisarenko's method.
//
//  The input sequence x is assumed to consist of p complex
//  exponentials in white noise.  The frequencies of the
//  complex exponentials and the variance of the white noise
//  are estimated using the Pisarenko harmonic decomposition.  
//
//  The frequency estimates are found from the peaks of the
//  pseudospectrum
//                               1
//              ------------------------------------
//              1 + a(1)exp(jw) + ... + a(p)exp(jpw)	
//
//   or from the roots of the polynomial formed from the 
//   vector a.  The estimate of the white noise variance is 
//   returned in sigma.
//
//  Calling Sequence: [vmin,sigma] = phd(x,p)
//
//  Inputs:
//
//    x - The input sequence.
//
//    p - The order of the all-pole model.
//
//  Outputs:
//
//    Px - The estimate of the power spectrum of
//    x(n) using a decibel scale.
//
//**********************************************************************
function [vmin,sigma] = phd(x,p)

  // Enforce a column vector.
  x = x(:);

  // Compute the autocorrelation matrix of order, p+1.
  R = Covar(x,p+1);

  // Compute the eigenvectors of R and the diagonal matrix of eigenvalues.
  [v,d] = spec(R);

  // Determine the variance of the noise.
  sigma = min(abs(diag(d)));

  // Find the eigenvector that spans the noise subspace.
  index = find(abs(diag(d)) == sigma);
 
  vmin = v(:,index);

endfunction

//**********************************************************************
//
//  Name:  ev
//
//  The purpose of this function is to estimate the spectrum of a
//  a process, x that consists of a sum of p complex exponentials
//  in white noise.  The frequencies of the complex exponentials and
//  the variance of the white noise are estimated using the
//  eigenvector method.
//  autocorrelation matrix.  The model for the process is such that
//  x(n) consists of a sum of complex exponentials in white noise.
//  After a principle components analysis is carried out, the
//  principal eigenvectors are used in the Blackman-Tukey estimate.
//
//  Calling Sequence: Px = ev(x,p,M)
//
//  Inputs:
//
//    x - The input sequence.
//
//    p - number of complex exponentials in x.
//
//    M - The size of the autocorrelation matrix to use in
//    estimating the complex exponential frequencies
//
//  Outputs:
//
//    Px - The pseudospectrum estimate using a dB scale.
//
//**********************************************************************
function Px = ev(x,p,M)

  // Enforce a column vector.
  x = x(:);

  if M >= p+1
    // Compute the autocorrelation matrix of order, M.
    R = Covar(x,M);

    // Compute the eigenvectors of R and the diagonal matrix of eigenvalues.
    [v,d] = spec(R);

    //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
    // Note that the sort() command sorts items in
    // descending order, therefore we have to provide
    // code to produce a sort in ascending order.
    //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
    // Sort the eigenvalues in decending order.
    [y,i] = sort(diag(d));

    // Reverse the sort to ascending order.
    y = y($:-1:1);
    i = i($:-1:1);
    //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

    // Clear initial sum.
    Px = 0;

    for j=1:M-p
      Px = Px + abs(matrixFft(v(:,i(j)),1024)).^2 / abs(y(j));
    end

    // Convert the pseudospectrum to a decibel scale.
    Px = -10 * log10(Px);
  else
    error('Specified size of R is too small');
  end

endfunction

//**********************************************************************
//
//  Name:  min_norm
//
//  The purpose of this function is to estimate the spectrum of a
//  a process, x that consists of a sum of p complex exponentials
//  in white noise.  The frequencies of the complex exponentials and
//  the variance of the white noise are estimated using the
//  minimum norm method.
//
//  Calling Sequence: Px = min_norm(x,p,M)
//
//  Inputs:
//
//    x - The input sequence.
//
//    p - number of complex exponentials in x.
//
//    M - The size of the autocorrelation matrix to use in
//    estimating the complex exponential frequencies
//
//  Outputs:
//
//    Px - The pseudospectrum estimate using a dB scale.
//
//**********************************************************************
function Px = min_norm(x,p,M)

  // Enforce a column vector.
  x = x(:);

  if M >= p+1
    // Compute the autocorrelation matrix of order, M.
    R = Covar(x,M);

    // Compute the eigenvectors of R and the diagonal matrix of eigenvalues.
    [v,d] = spec(R);

    //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
    // Note that the sort() command sorts items in
    // descending order, therefore we have to provide
    // code to produce a sort in ascending order.
    //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
    // Sort the eigenvalues in decending order.
    [y,i] = sort(diag(d));

    // Reverse the sort to ascending order.
    y = y($:-1:1);
    i = i($:-1:1);
    //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

    // We have to kludge this for Scilab.
    V = v(:,i(1));

    for j = 2:M-p
      V = [V,v(:,i(j))];
    end

    // Compute minimum norm solution.
    a = V * V (1,:)';

    // Estimate the pseudospectrum on a decibel scale.
    Px = -20 * log10(abs(matrixFft(a,1024)));
  else
    error('Specified size of R is too small');
  end

endfunction

//**********************************************************************
//
//  Name:  music
//
//  The purpose of this function is to estimate the frequencies
//  and powers of a sum of complex sinusids in white noise, using
//  the MUSIC algorithm.  Additionally, the variance of the white
//  noise is estimated.
//
//  Calling Sequence: Px = music(x,p,M)
//
//  Inputs:
//
//    x - The input sequence.
//
//    p - The number of complex exponentials to find.
//
//    M - The number of noise eigenvectors to use.
//
//  Outputs:
//
//    Px - The pseudospectrum.
//
//**********************************************************************
function Px = music(x,p,M)

  // Enforce a column vector.
  x = x(:);

  if M >= p+1 & length(x) >= M
    // Compute the autocorrelation matrix of order, M.
    R = Covar(x,M);

    // Compute the eigenvectors of R and the diagonal matrix of eigenvalues.
    [v,d] = spec(R);

    //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
    // Note that the sort() command sorts items in
    // descending order, therefore we have to provide
    // code to produce a sort in ascending order.
    //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
    // Sort the eigenvalues in decending order.
    [y,i] = sort(diag(d));

    // Reverse the sort to ascending order.
    y = y($:-1:1);
    i = i($:-1:1);
    //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

    // Clear initial sum.
    Px = 0;

    for j = 1:M-p
      Px = Px + abs(matrixFft(v(:,i(j)),1024));
    end

    // Convert to a decibel scale.
    Px = -20 * log10(Px);
  else
    error('Size of R is inappropriate');
  end

endfunction

//**********************************************************************
//
//  Name:  bt_pc
//
//  The purpose of this function is to estimate the spectrum of a
//  a process, x using a principal components analysis of the
//  autocorrelation matrix.  The model for the process is such that
//  x(n) consists of a sum of complex exponentials in white noise.
//  After a principle components analysis is carried out, the
//  principal eigenvectors are used in the Blackman-Tukey estimate.
//
//  Calling Sequence: Px = bt_pc(x,p,M)
//
//  Inputs:
//
//    x - The input sequence.
//
//    p - number of complex exponentials in x.
//
//    M - The size of autocorrelation matrix.
//
//  Outputs:
//
//    Px - The spectrum estimate using a dB scale.
//
//**********************************************************************
function Px = bt_pc(x,p,M)

  // Enforce a column vector.
  x = x(:);

  if M >= p+1
    // Compute the autocorrelation matrix of order, M.
    R = Covar(x,M);

    // Compute the eigenvectors of R and the diagonal matrix of eigenvalues.
    [v,d] = spec(R);

    //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
    // Note that the sort() command sorts items in
    // descending order, therefore we have to provide
    // code to produce a sort in ascending order.
    //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
    // Sort the eigenvalues in decending order.
    [y,i] = sort(diag(d));

    // Reverse the sort to ascending order.
    y = y($:-1:1);
    i = i($:-1:1);
    //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

    // Clear initial sum.
    Px = 0;

    for j = M-p+1:M
      Px = Px + abs(matrixFft(v(:,i(j)),1024)) * sqrt(real(y(j)));
    end

    // Convert to a decibel scale.
    Px = 20 * log10(Px) - 10 * log10(M);
  else
    error('Specified size of R is too small');
  end

endfunction

//**********************************************************************
//
//  Name: Aic
//
//  The purpose of this function is to determine the error criteria of
//  a filter order using the Akaike Information Criterion.
//
//  Calling Sequence: c = Aic(ep,p,N)
//
//  Inputs:
//
//    ep - The modeling error for a given filter order.
//
//    p - The filter order that is getting tested.
//
//    N - The data record length of the input sequence being
//    evaluated.
//
//  Outputs:
//
//    c - C(p), the error criteria which includes the penalty term.
//
//**********************************************************************
function c = Aic(ep,p,N)

  // Compute the error criteria.
  c = N * log(ep) + 2 * p;

endfunction

//**********************************************************************
//
//  Name: Mdl
//
//  The purpose of this function is to determine the error criteria of
//  a filter order using the minimum description length.
//
//  Calling Sequence: c = Mdl(ep,p,N)
//
//  Inputs:
//
//    ep - The modeling error for a given filter order.
//
//    p - The filter order that is getting tested.
//
//    N - The data record length of the input sequence being
//    evaluated.
//
//  Outputs:
//
//    c - C(p), the error criteria which includes the penalty term.
//
//**********************************************************************
function c = Mdl(ep,p,N)

  // Compute the error criteria.
  c = N * log(ep) + log(N) * p;

endfunction

//**********************************************************************
//
//  Name: Fpe
//
//  The purpose of this function is to determine the error criteria of
//  a filter order using Akaike's Final Prediction Error.
//
//  Calling Sequence: c = Fpe(ep,p,N)
//
//  Inputs:
//
//    ep - The modeling error for a given filter order.
//
//    p - The filter order that is getting tested.
//
//    N - The data record length of the input sequence being
//    evaluated.
//
//  Outputs:
//
//    c - C(p), the error criteria which includes the penalty term.
//
//**********************************************************************
function c = Fpe(ep,p,N)

  // Compute the error criteria.
  c = ep * ((N + p + 1) / (N - p - 1));

endfunction

//**********************************************************************
//
//  Name: Cat
//
//  The purpose of this function is to determine the error criteria of
//  a filter order using Parzen's Criterion Autoregressive Transfer
//  function.
//
//  Calling Sequence: c = Cat(ep,p,N)
//
//  Inputs:
//
//    ep - The modeling error vector for a given filter order.
//
//    p - The filter order that is getting tested.
//
//    N - The data record length of the input sequence being
//    evaluated.
//
//  Outputs:
//
//    c - C(p), the error criteria which includes the penalty term.
//
//**********************************************************************
function c = Cat(ep,p,N)

  // Initialize sum.
  c = 0;

  // Compute the error criteria.
  for j = 1:p
    c = c + ((N - j) / (N * ep(j))); 
  end

  // Compute the average.
  c = c / N;

  // Include final term.
  c = c - (N - p) / (N * ep(p));

endfunction

