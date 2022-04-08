//**********************************************************************
// File Name: SpectrumEstimation.sci
//**********************************************************************

//**********************************************************************
//
//  Name: per
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
//    Px - The periodogram estimate of the power spectrum of x(n).
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
//    of x(n).
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

  // Compute L2 norm of w(n).
  U  = norm(w)^2 / N;

  // Compute windowed sequence.
  xw = x(n1:n2) .* w' / norm(w);

  // Compute the periodogram.
  Px = N * per(xw);

endfunction






