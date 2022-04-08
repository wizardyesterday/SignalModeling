//**********************************************************************
// File Name: SpectrumEstimation.sci
//**********************************************************************

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
//  The purpose of this function is to compute the spectrum of a
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
  // Compute the periodogram.
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Compute Fourier transform.
  X = fft(x(n1:n2),-1);

  // Compute power spectral density.
  Px = (X .* conj(X)) / (n2 - n1 + 1) ^2;

  Px(1) = Px(2);
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

endfunction

//**********************************************************************
//
//  Name: mper
//
//  The purpose of this function is to compute the spectrum of a
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
      w  = window('re',N);

    case 2
      // Select Hamming window.
      w = window('hm',N);

    case 3
      // Select Hanning window.
      w = window('hn',N);

    case 4
      // Select Bartlett window.
      w = window('tr',N);

    case 5
      // Select Blackman window.
      w = blackman(N);

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
//  The purpose of this function is to compute the spectrum of a
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
//  Name:  bart
//
//  The purpose of this function is to compute the spectrum of a
//   process using Welch's method of modified periodogram averaging.
//
//  Calling Sequence: x = welch(x,L,over,win)
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
    // Set initial limits for the first subsequence.
    n1 = 1;
    n2 = L;

    // Set the increment for the next subsequence.
    n0 = (1 - over) * L;

    // Compute the number of subsequences.
    nsect = 1+ floor((length(x) -L ) / n0);

    // Clear average.
    Px=0;

    for i = 1:nsect
      // Add the next periodogram to the average.
      Px = Px + mper(x,win,n1,n2) / nsect;

      // Reference the next subsequence.
      n1 = n1 + n0;  
      n2 = n2 + n0;
    end
  else
    error('Overlap is invalid');
  end

endfunction

