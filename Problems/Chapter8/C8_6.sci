//**********************************************************************
// File Name: C8_6.sci
//**********************************************************************

// Bring in all of the utility functions.
exec('utils.sci',-1);

//**********************************************************************
// Mainline code.
//**********************************************************************
// Construct filter for white noise process.
b = [1 -1];
b = convol(b,[1 1]);
b = convol(b,[1 1]);
b = convol(b,[1 1]);

// Generate white noise sequence with unit variance.
noisegen(1,1024,1);
v = feval([1:1024],Noise);
v = v(:);

// Generate MA(4) process.
w = filterBlock(v,b,0);

// Generate random phases uniformly distributed between 0 and 2PI.
phi1 = rand() * 2*%pi;
phi2 = rand() * 2*%pi;

// Generate time vector.
n = 0:length(v) - 1;
n = n(:);

// Set frequencies of sinusoids.
w1 = %pi / 2;
w2 = 1.1 * %pi/2;

// Construct random process.
x = 2*cos(w1*n + phi1) + 2*cos(w2*n + phi2) + w;

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (a): The power spectrum is first plotted.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Construct exact power spectrum.
Px = constructPowerSpectrum(x);
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (b): We find the variance of w(n) and compare to the
// power in each of the sinusoids.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
sigmaW2 = variance(w);

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (c): No code needed.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (d): 
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (e): Generate a sample realization of x(n) of length,
// N = 256, and plot the periodogram.  Repeat for 20 different
// realizations of the process.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Generate time vector.
n_e = 0:255;
n_e = n_e(:);

// Generate random phases.
phi1_e = rand() * 2*%pi;
phi2_e = rand() * 2*%pi;

v_e = feval([1:256],Noise);
v_e = v_e(:);

// Generate MA(4) process.
w_e = filterBlock(v_e,b,0);

// Construct realization of random process.
x_e = 2*cos(w1*n_e + phi1_e) + 2*cos(w2*n_e + phi2_e) + w_e;

// Generate initial periodogram.
Pxper = per(x_e);

// Allocate the matrix.
Pxper_mat = zeros(1024,1);

// Generate 20 different realiations.
for j = 1:20
  // Generate random phases.
  phi1_e = rand() * 2*%pi;
  phi2_e = rand() * 2*%pi;

  // Generate white noise sequence with unit variance.
  v_e = feval([1:256],Noise);
  v_e = v_e(:);

  // Generate MA(4) process.
  w_e = filterBlock(v_e,b,0);

  // Construct realization of random process.
  x_e = 2*cos(w1*n_e + phi1_e) + 2*cos(w2*n_e + phi2_e) + w_e;

  Pxper_mat(:,j) = per(x_e);
end

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (f):
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Plot results.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

scf(1);

// Part (a).
subplot(311);
title('Exact spectrum');
plot(20*log10(Px));

// Part (e).
subplot(312);
title('Periodogram of Length 256');
plot(20*log10(Pxper));

subplot(313);
title('20 Realizations of Periodogram of Length 256');
plot(20*log10(Pxper_mat));


