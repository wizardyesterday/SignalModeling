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
n_256 = 0:255;
n_256 = n_256(:);

// Generate random phases.
phi1_e = rand() * 2*%pi;
phi2_e = rand() * 2*%pi;

v_e = feval([1:256],Noise);
v_e = v_e(:);

// Generate MA(4) process.
w_e = filterBlock(v_e,b,0);

// Construct realization of random process.
x_e = 2*cos(w1*n_256 + phi1_e) + 2*cos(w2*n_256 + phi2_e) + w_e;

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
  x_e = 2*cos(w1*n_256 + phi1_e) + 2*cos(w2*n_256 + phi2_e) + w_e;

  Pxper_mat(:,j) = per(x_e);
end

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (f): Estimate the expected value of the periodogram by
// averaging the periodograms that are formed from an ensemble
// of 50 realizations of x(n).
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Allocate the matrix.
Pxper_mat_f = zeros(1024,1);

// Generate 50 different realiations.
for j = 1:50
  // Generate random phases.
  phi1_f = rand() * 2*%pi;
  phi2_f = rand() * 2*%pi;

  // Generate white noise sequence with unit variance.
  v_f = feval([1:256],Noise);
  v_f = v_f(:);

  // Generate MA(4) process.
  w_f = filterBlock(v_e,b,0);

  // Construct realization of random process.
  x_e = 2*cos(w1*n_256 + phi1_f) + 2*cos(w2*n_256 + phi2_f) + w_f;

  Pxper_mat_f(:,j) = per(x_e);
end

// Compute the average power spectrum.
Pxperavg = sum(Pxper_mat_f,'c') / 50;

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (g): Repeat parts (e) and (f) using Bartlett's method
// with K = 2, 4, and 8 sections..
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Allocate the matrices
PxB_20_2 = zeros(1024,1);
PxB_20_4 = zeros(1024,1);
PxB_20_8 = zeros(1024,1);
PxB_50_2 = zeros(1024,1);
PxB_50_4 = zeros(1024,1);
PxB_50_8 = zeros(1024,1);

// Generate 20 different realiations.
for j = 1:20
  // Generate random phases.
  phi1_ge = rand() * 2*%pi;
  phi2_ge = rand() * 2*%pi;

  // Generate white noise sequence with unit variance.
  v_ge = feval([1:256],Noise);
  v_ge = v_ge(:);

  // Generate MA(4) process.
  w_ge = filterBlock(v_ge,b,0);

  // Construct realization of random process.
  x_ge = 2*cos(w1*n_256 + phi1_ge) + 2*cos(w2*n_256 + phi2_ge) + w_ge;

  PxB_20_2(:,j) = bart(x_ge,2);
  PxB_20_4(:,j) = bart(x_ge,4);
  PxB_20_8(:,j) = bart(x_ge,8);
end

// Generate 50 different realiations.
for j = 1:50
  // Generate random phases.
  phi1_gf = rand() * 2*%pi;
  phi2_gf = rand() * 2*%pi;

  // Generate white noise sequence with unit variance.
  v_gf = feval([1:256],Noise);
  v_gf = v_gf(:);

  // Generate MA(4) process.
  w_gf = filterBlock(v_gf,b,0);

  // Construct realization of random process.
  x_gf = 2*cos(w1*n_256 + phi1_gf) + 2*cos(w2*n_256 + phi2_gf) + w_gf;

  PxB_50_2(:,j) = bart(x_gf,2);
  PxB_50_4(:,j) = bart(x_gf,4);
  PxB_50_8(:,j) = bart(x_gf,8);
end

// Compute the average power spectra.
PxBavg_2 = sum(PxB_50_2,'c') / 50;
PxBavg_4 = sum(PxB_50_4,'c') / 50;
PxBavg_8 = sum(PxB_50_8,'c') / 50;

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Plot results.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

scf(1);

// Part (a).
subplot(411);
title('Exact spectrum');
plot(20*log10(Px));

// Part (e).
subplot(412);
title('Periodogram of Length 256');
plot(20*log10(Pxper));

subplot(413);
title('20 Realizations of Periodogram of Length 256');
plot(20*log10(Pxper_mat));

// Part (f).
subplot(414);
title('50 Realizations of Periodogram of Length 256 (Average)');
plot(20*log10(Pxperavg));

// Part (g).
scf(2);

subplot(321);
title('Overlay Bartlett, Ensemble: 20, 2 Sections');
plot(20*log10(PxB_20_2));

subplot(323);
title('Overlay Bartlett, Ensemble: 20, 4 Sections');
plot(20*log10(PxB_20_4));

subplot(325);
title('Overlay Bartlett, Ensemble: 20, 8 Sections');
plot(20*log10(PxB_20_8));

subplot(322);
title('Average Bartlett, Ensemble: 50, 2 Sections');
plot(20*log10(PxBavg_2));

subplot(324);
title('Average Bartlett, Ensemble: 50, 4 Sections');
plot(20*log10(PxBavg_4));

subplot(326);
title('Average Bartlett, Ensemble: 50, 8 Sections');
plot(20*log10(PxBavg_8));





