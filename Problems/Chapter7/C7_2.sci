//**********************************************************************
// File Name: C7_2.sci
//**********************************************************************

// Bring in all of the utility functions.
exec('utils.sci',-1);

//**********************************************************************
// Mainline code.
//**********************************************************************

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (a): Generate 500 samples of the processes, x(n) and
// v2(n).
// x(n = d(n) + g(n),
// d(n) = sin(n*w0 + phi),
// where w0 = 0.05*pi, and phi is a random variable that is
// uniformly distributed between -pi and pi.  Assume that g(n) is
// unit variance white noise.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Generate time vector.
n = 0:499;
n = n(:);

// Generate random phase angla sequence.
phi = generateBipolarNoise(%pi,500);
phi = phi(:);

// Generate white noise sequence with unit variance.
noisegen(1,500,1);
g = feval([1:500],Noise);
g = g(:);

// Construct the desired signal.
w0 = 0.05 * %pi;
d = sin(n*w0 + phi);
d = d(:);

// Construct measured signal.
x = d + g;

