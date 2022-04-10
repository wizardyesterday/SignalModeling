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
  r = convol(x(n1:n2),x(n2:-1:1));
  kmax = find(r == max(r));
  r = r(kmax:$);

  // Truncate the autocorrelation sequence to M lags.
  r = r(1:M);

  // Make the autocorrelation sequence double-sided.
  r = [r(M:-1:2) r];

  // Force column vector.
  r = r(:);

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
//    2. Use the toeplitz function to construct the autocorrelation
//    matrix.
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
    r = convol(x,x($:-1:1));
    kmax = find(r == max(r));
    r = r(kmax:$);

    // Truncate the autocorrelation sequence to p lags.
    r = r(1:p);

    // Create the autocorrelation matrix.
    R = toeplitz(r);
  else
    error('Order too large');
  end

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
//   process using the minimum variance method.
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
//    Px - The Welch estimate of the power spectrum of x(n) using
//    a linear scale.
//
//**********************************************************************
function [v,V1] = minvar(x,p)

  // Enforce a column vector.
  x = x(:);

  // Compute the autocorrelation matrix of order, p.
  R = Covar(x,p);

  // Compute the eigenvectors of R and the diagonal matrix of eigenvalues.
  [v,d] = spec(R);

  U = diag(inv(abs(d)+ %eps));

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Compute the number of colums in the matrix
  // of eigenvectors.
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  rowsColumns = size(v);
  columns = rowsColumns(2);
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/


  V1 = zeroPadSequence(v(:,1),1024);

  for i = 2:columns
    V1(:,i) = zeroPadSequence(v(:,i),1024);
  end

//  V  = abs(fft(v,1024)).^2;
//  V  = abs(fft(v,-1)).^2;

//  Px = 10*log10(p) - 10*log10(V * U);

endfunction

