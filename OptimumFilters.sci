//**********************************************************************
// File Name: OptimumFilters.sci
//**********************************************************************

//**********************************************************************
//
//  Name: FirWienerFilter
//
//  Purpose: The purpose of this function is to design a pth-order FIR
//  Wiener filter.
//
//  Calling Sequence: [w,e] = FirWienerFilter(rd,sigmav2,p)
//
//  Inputs:
//
//    rd - A sequence of autocorrelation values for the desired
//    signal, d(n).
//
//    sigmav2 - The variance of the noise signal that perturbs d(n).
//
//    p - The order of the desired Wiener filter.
//
//  Outputs:
//
//    w - The Wiener filter coefficients.
//
//    e - The mean-square error associated with the filter.
//
//**********************************************************************
function [w,e] = FirWienerFilter(rd,sigmav2,p)

  // Force column a vector.
  rd = rd(:);

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Construct the autocorrelation sequence of
  // x(n) = d(n) + v(n).  The result will be
  // rx(k) = rd(k) + sigmav2 * delta(k),
  // where delta(k) is the delta function.
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  rx = rd(1:p);
  rx(1) = rx(1) + sigmav2;
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

  // Construct autocorrelation matrix.
  Rx = toeplitz(rx);

  // Compute Wiener filter coefficients noting that rdx(k) = rd(k)
  w = Rx \ rd(1:p);

  // Compute mean-square error.
  e = rd(1) - w' * rd(1:p);

endfunction

//**********************************************************************
//
//  Name: WienerNoiseCanceller
//
//  Purpose: The purpose of this function is to design a pth-order FIR
//  Wiener filter.  Here's how things work.  There exists a desired
//  signal, d(n), that is measured in the presence of noise, this
//  signal being x(n) = d(n) + v(n), such that v1(n) is white noise.
//  This tainted signal might represent a voice signal, d(n) in a
//  cockpit of a fighter jet, and v1(n) might represent wind noise and
//  engine noise.  A secondary sensor is placed away from the desired
//  signal with the purpose of measuring the background noise. Let it
//  be supposed that this signal is called v2(n), and although v2(n)
//  is correlated with v1(n), the signals are not the same.  Now,
//  v2(n) is passed through a Wiener filter, and the output of the filter
//  is an estimate of v1(n).  Let's calle it v1hat(n).  Finally,
//  v2hat(n) is subtracted from x(n), and this provides an estimate of
//  d(n).  Let's call this estimate dhat(n).
//  noise signal from the secondary sensor is pas
//
//  Calling Sequence: [w,e] = WienerNoiseCanceller(x,v1,v2,p)
//
//  Inputs:
//
//    x - An observed signal represented by, x(n) = d(n) + v1(n).
//
//    v1 - The noise, v1(n), that is combined with the observed signal.
//
//    v2 - A signal, v2(n), that is measured by a secondary sensor.
//
//    p - The order of the desired Wiener filter.
//
//  Outputs:
//
//    dhat - The estimate, dhat(n), of the desired signal.
//
//    e - The mean-square error associated with the filter.
//
//**********************************************************************
function [dhat,v1hat,e] = WienerNoiseCanceller(x,v1,v2,p)

  //Force column vectors.
  x = x(:);
  v1 = v1(:);
  v2 = v2(:);

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Compute autocorrelation sequences, and cross-
  // correlation sequences.  Note that v1(n) is
  // not really available in a real system.  We're
  // keeping it available only so that the mean-
  // square error can be computed.  Note also, that
  // v2(n) is measured from a secondary sensor
  // whose sole purpose is for the measurement of
  // noise.
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Compute autocorrelation of v1(n).
  rv1 = convol(v1,v1($:-1:1));
  kmax = find(rv1 == max(rv1));
  rv1 = rv1(kmax:kmax+p-1);

  // Compute autocorrelation of v2(n).
  rv2 = convol(v2,v2($:-1:1));
  kmax = find(rv2 == max(rv2));
  rv2 = rv2(kmax:kmax+p-1);

  // Compute autocorrelation of d(n);
  rd = convol(d,d($:-1:1));
  kmax = find(rd == max(rd));
  rd = rd(kmax:$);

  // Construct cross-correlation of x(n) and v2(n).
  rv1v2 = convol(x,v2($:-1:1));
  kmax = find(rv1v2 == max(rv1v2));
  rv1v2 = rv1v2(kmax:kmax+p-1);
  rv1v2 = rv1v2(:);
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

  // Construct autocorrelation matrix.
  Rv2 = toeplitz(rv2);

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // In this block of code, the Wiener filter
  // coefficients are first computed.  Once
  // computed, the noise signal, v2(n), from the
  // secondary sensor is filtered by the Weiner
  // filter.  The result is an estimate of v1(n),
  // called v1hat(n).
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Compute the Wiener filter coefficients.
  w = Rv2 \ rv1v2(1:p);

  // Compute estimate of v1(n).
  v1hat = filterBlock(v2,w,0);
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // The estimate of d(n) is next computed by
  // subtracting the estimate of d1(n) from x(n).
  // The result is the desired signal with noise
  // removed.  This signal is called dhat(n).
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  dhat = x - v1hat;
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

  // Compute mean-square error.
  e = rv1(1) - w' * rv1v2(1:p);

endfunction

//**********************************************************************
//
//  Name: kalmanFilter
//
//  Purpose: The purpose of this function is to perform Kalman
//  filtering.
//
//  Calling Sequence: xhat = kalmanFilter(N,A,C,Qv,Qw,P_0_0)
//
//  Inputs:
//
//    Y - A signal that contains x(n) + noise, v(n).
//
//    N - The number of iterations for the Kalman gain computation.
//
//    A - The state transition matrix.
//
//    C - The observation matrix.
//
//    Qv - The measurement noise.
//
//    Qw - The process noise.
//
//    P_0_0 - The initial covariance estimate, P(0|0).
//
//  Outputs:
//
//    xhat - The estimated values.
//
//**********************************************************************
function xhat = kalmanFilter(y,N,A,C,Qv,Qw,x_0_0,P_0_0)

  // Retrieve the initial values.
  x = x_0_0;
  P = P_0_0;

  for n = 2:N
    // Update x(n|n-1).
    x_n_nm1 = A * x;

    // Update P(n|n-1).
    P_n_nm1 = (A * P * A') + Qw;

    // Update the Kalman gain.
    K = (P_n_nm1 * C') / ((C * P_n_nm1 *C') + Qv);

    // Update x(n|n).
    x = x_n_nm1 + K * (y(n) - (C * x_n_nm1));

    // Update P(n|n).
    P = P_n_nm1 - (K * C * P_n_nm1);

    // Construct returned estimate.
    xhat(n-1) = x(1);
  end

endfunction

//**********************************************************************
//
//  Name: stateTransitionMatrix
//
//  Purpose: The purpose of this function is to create the state
//  transition matrix that corresponds to a set of recursive filter
//  coefficients.
//
//  Calling Sequence: A = stateTransitionMatrix(a)
//
//  Inputs:
//
//    a - An input vector of recursive coefficients of a difference
//    equation.
//
//  Outputs:
//
//    A - The state transision matrix that is constructed from the
//    input vector, a.
//
//**********************************************************************
function A = stateTransitionMatrix(a)

  // Force a colunn vector.
  a = a(:);

  // Convert to a row vector.
  a = a';

  n = length(a) - 1;

  // Create an identify submatrix.
  Ia = eye(n,n);

  // Add a column of zeros to the identity matrix.
  Ia = [Ia zeros(n,1)];

  // Construct the final matrix.
  A = [a; Ia];

endfunction

