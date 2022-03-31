//**********************************************************************
// File Name: C7_3.sci
//**********************************************************************

exec('utils.sci',-1);

//**********************************************************************
//
//  Name: computeKalmanGain
//
//  Purpose: The purpose of this function is to create lowpass filters
//  using the Pade approximation.  Plots are created for various values
//  of p and q.  The heavy lifting is performed by the pade() function.
//
//  Calling Sequence: K = computeKalmanGain(N,A,C,Qv,Qw,P_0_0)
//
//  Inputs:
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
//    K - The Kalman gain matrix.
//
//**********************************************************************
function K = computeKalmanGain(N,A,C,Qv,Qw,P_0_0)

  // Retrieve the initial value.
  P = P_0_0;

  for n = 1:N
    // Update P(n|n-1).
    P_n_nm1 = A * P * A' + Qw;

    // Update the Kalman gain.
    K(:,:,n) = (P_n_nm1 * C') / ((C * P_n_nm1 * C') + Qv);

    // Update P(n|n).
    P = P_n_nm1 - (K(:,:,n) * C * P_n_nm1);
  end

endfunction

//**********************************************************************
// Mainline code.
//**********************************************************************
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (b): Compute the Kalman gain for for the process,
// x(n) = -0.1x(n-1) - 0.09x(n-2) + 0.648x(n-3) + w(n)
// Both w(n) is a white noise process withsigmaW^2 = 1.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Filter coefficient for AR(3) process.
a3 = [0.1 0.09 -0.648];

// Set parameters.
A3 = stateTransitionMatrix(a3);
C3 = [1 0 0];
Qv = 0.64;
Qw = 1;
x3_0_0 = [1 0 0]';
P3_0_0 = [1 0 0; 0 0 0; 0 0 0];

K = computeKalmanGain(11,A3,C3,Qv,Qw,P3_0_0);

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (c): Determine steady-state value for the Kalman gain.
// Try different values of P(0|0).
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
K2 = computeKalmanGain(11,A3,C3,Qv,Qw,2*P3_0_0);
K3 = computeKalmanGain(11,A3,C3,Qv,Qw,4*P3_0_0);

disp('K K2 K3');
disp([K K2 K3]);

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (d): Generate:
// x(n) = -0.1x(n-1) - 0.09x(n-2) + 0.648x(n-3) + w(n)
// y(n) = x(n) + v(n).
// Both w(n) and v(n) are white noise sequences with
// sigmaW^2 = 1 and sigmaV^2 = 0.64, respectively.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Generate time vector.
// Generate time vector.
n = 0:499;

// Generate white noise sequence with unit variance.
noisegen(1,500,1);
w = feval([1:500],Noise);
w = w(:);

// Generate white noise sequence with a variance of 0.64;
noisegen(1,500,sqrt(0.64));
v = feval([1:500],Noise);
v = v(:);

// Generate third-order autoregressive process.
a1 = [0.1 0.09 -0.648];
x = filterBlock(w,1,a1) + w;

// Compute observations.
y = x + v;

