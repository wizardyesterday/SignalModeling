//**********************************************************************
// File Name: C8_9.sci
//**********************************************************************

// Bring in all of the utility functions.
exec('utils.sci',-1);

//**********************************************************************
// Mainline code.
//**********************************************************************
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (b): Generate the 128x128 autocorrelation matrix for the
// process with spectral density Px(w) = 1; -PI/2 < |w| < PI/2.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Generate time vector.
k = 0:127;
k = k(:);
k(1) = 1e-6;

num = sin(%pi*k/2);
den = %pi*k/2;

// Compute autocorrelation sequence.
rx_a = 2*num ./ den;

// Construct autocorrelation matrix.
Rx_a = toeplitz(rx_a);

// Compute eigenvalues
lamda_a = spec(Rx_a);

// Compute sum of log of eigenvalues.
s_a = sum(log(lamda_a));
s_a = real(s_a);
s_a = s_a / 128;

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (b): Repeat parts (a) and (b) for the harmonic process,
// x(n = A*cos(n*w0 + phi) + w(n), where w(n) is unit variance
// white noise.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Initialize white noise generator with unit variance.
noisegen(1,1024,1);

// Generate white noise sequence.
w = feval([1:128],Noise);
w = w(:);

// Choose parameters.
A = 4;
w0 = %pi/16;

// Generate signal.
x = A*cos(k*w0) + w;

// Generate autocorrelation sequence.
rx_b = convol(x,x($:-1:1));
nMax = find(rx_b == max(rx_b));
rx_b = rx_b(nMax:$);

// Construct autocorrelation matrix.
Rx_b = toeplitz(rx_b);

// Compute eigenvalues
lamda_b = spec(Rx_b);

// Compute sum of log of eigenvalues.
s_b = sum(log(lamda_b));
s_b = real(s_b);
s_b = s_b / 128;

sideal_b = log(%pi*A*A/2)/%pi;


//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Plot results.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

scf(1);

// Part (a).
subplot(211);
title('Eigenvalues For Constant Spectrum');
plot(lamda_a);

// Part (b).
subplot(212);
title('Eigenvalues For Harmonic Process');
plot(lamda_b);
