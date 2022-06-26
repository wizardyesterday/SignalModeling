//**********************************************************************
// File Name: C8_13.sci
//**********************************************************************

// Bring in all of the utility functions.
exec('utils.sci',-1);

//**********************************************************************
//
//  Name:  aicEigen
//
//  The purpose of this function is to compute the AIC metric for
//  the eigenvalues of a random process that contains a set of
//  complex exponentials.
//
//  Calling Sequence: metric = aicEigen(k,M,N,lamda)
//
//  Inputs:
//
//    k - An estimate of the number of complex exponentials in a
//   random process.
//
//    M - The number of eigenvalues.
//
//    N - The length of the data record for which the autocorrelation
//    matrix was generated.   
//
//    lamda - A vector of eigenvalues associated with the MxM
//    autocorrelation matrix for a random process.
//
//  Outputs:
//
//    V - The discrete time Fourier transform of v.
//
//**********************************************************************
function metric = aicEigen(k,M,N,lamda)

  // Compute the harmonic mean, alpha(k).
  Alpha = prod(lamda(k+1:M))^(1/(M - k));

  // Compute arithmetric mean, beta(k).
  Beta = sum(lamda(k+1:M)) / (M - k);

  // We can now compute the metric.
  metric = -N*(M - k)*log(Alpha / Beta) + k*(2*M - k);

endfunction

//**********************************************************************
// Mainline code.
//**********************************************************************
// Initialize unit variance noise source.
noisegen(1,100,1);

// Set up parameters
A1 = 2;
A2 = 2;
A3 = 2;
w1 = 0.4 * %pi;
w2 = 0.5 * %pi;
w3 = 1.2 * %pi;

// Generate time vector.
n = 0:99;

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (a): Generate x(n), construct a 6x6 autocorrelation
// matrix, and compute the eigenvalues.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Generate white noise process.
w = feval([1:100],Noise);

// Generate random phases.
phi1 = generateBipolarNoise(%pi,1);
phi2 = generateBipolarNoise(%pi,1);
phi3 = generateBipolarNoise(%pi,1);

phi1 = 0;
phi2 = 0;
phi3 = 0;

// Generate realization.
x = A1*exp(%i*w1*n + phi1) + + A2*exp(%i*w2*n + phi2) ...
    + A3*exp(%i*w3*n + phi3) + w;

// Construct 6x6 autocorrelation matrix.
R6 = Covar(x,6);

// Compute eigenvalues.
lamda6 = spec(R6);
lamda6 = clean(lamda6);

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (b): Repeat part (a) for autocorrelation matrices of
// size 15x15 and 30x30.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Construct autocorrelation matrices.
R15 = Covar(x,15);
R30 = Covar(x,30);

// Compute eigenvalues.
lamda15 = spec(R15);
lamda15 = clean(lamda15);
lamda30 = spec(R30);
lamda30 = clean(lamda30);

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (c): Experiment with different variances.  The goal is
// to see if the variance of white noise affects the estimate
// for the number of complex exponentials in the process.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/


//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (d): Perform AIC analysis.  The goal is to see how well
// the AIC metric is in determining how many complex exponentials
// are contained in the process.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
for k = 1:4
  // Compute the vector of metrics.
  metric6(k) = aicEigen(k,6,100,lamda6);
end

// Eliminate any residual imaginary part.
metric6 = real(metric6);

// Estimate the number of complex exponentials.
p6 = find(metric6 == min(metric6));

for k = 1:10
  // Compute the vector of metrics.
  metric15(k) = aicEigen(k,15,100,lamda15);
end

// Eliminate any residual imaginary part.
metric15 = real(metric15);

// Estimate the number of complex exponentials.
p15 = find(metric15 == min(metric15));

for k = 1:20
  // Compute the vector of metrics.
  metric30(k) = aicEigen(k,30,100,lamda30);
end

// Eliminate any residual imaginary part.
metric30 = real(metric30);

// Estimate the number of complex exponentials.
p30 = find(metric30 == min(metric30));

// Display results.
disp([p6; p15; p30]);



