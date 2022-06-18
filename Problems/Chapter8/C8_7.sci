//**********************************************************************
// File Name: C8_7.sci
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

// Initialize white noise generator with unit variance.
noisegen(1,1024,1);

// Generate time vector.
n_256 = 0:255;
n_256 = n_256(:);

// Set frequencies of sinusoids.
w1 = %pi / 2;
w2 = 1.1 * %pi/2;

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (a): Generate a sample realization of x(n) of length
// N = 256, and estimate the spectrum using the maximum entropy
// method with p = 4, 8, 16, and 32.  Repeat for 20 different
// realizations of the process, and generate an overlay plot of
// the estimates along with the ensemble average.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Preallocate matrices.
Pmem_20_4 = zeros(1024,1);
Pmem_20_8 = zeros(1024,1);
Pmem_20_16 = zeros(1024,1);
Pmem_20_32 = zeros(1024,1);

for j = 1:20
  // Generate random phases.
  phi1 = rand() * 2*%pi;
  phi2 = rand() * 2*%pi;

  // Generate white noise sequence.
  v = feval([1:256],Noise);
  v = v(:);

  // Generate MA(4) process.
  w = filterBlock(v,b,0);

  // Construct random process.
  xmem = 2*cos(w1*n_256 + phi1) + 2*cos(w2*n_256 + phi2) + w;

  // Construct spectral estimates.
  Pmem_20_4(:,j) = mem(xmem,4);
  Pmem_20_8(:,j) = mem(xmem,8);
  Pmem_20_16(:,j) = mem(xmem,16);
  Pmem_20_32(:,j) = mem(xmem,32);
end

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (b): Repeat part (a) using the Burg algorithm and the
// modified covariance method.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

//----------------------------------------------------------------
// Use Burg method.
//----------------------------------------------------------------
// Preallocate matrices.
Pburg_20_4 = zeros(1024,1);
Pburg_20_8 = zeros(1024,1);
Pburg_20_16 = zeros(1024,1);
Pburg_20_32 = zeros(1024,1);

for j = 1:20
  // Generate random phases.
  phi1 = rand() * 2*%pi;
  phi2 = rand() * 2*%pi;

  // Generate white noise sequence.
  v = feval([1:256],Noise);
  v = v(:);

  // Generate MA(4) process.
  w = filterBlock(v,b,0);

  // Construct random process.
  xburg = 2*cos(w1*n_256 + phi1) + 2*cos(w2*n_256 + phi2) + w;

  // Generate spectral estimate using the Burg method.
  [gamm,e] = burg(xburg,4);
  aburg_4 = gtoa(gamm);
  [gamm,e] = burg(xburg,8);
  aburg_8 = gtoa(gamm);
  [gamm,e] = burg(xburg,16);
  aburg_16 = gtoa(gamm);
  [gamm,e] = burg(xburg,32);
  aburg_32 = gtoa(gamm);

  // Construct spectral estimates.
  Pburg_20_4(:,j) = constructPowerSpectrum(aburg_4);
  Pburg_20_8(:,j) = constructPowerSpectrum(aburg_8);
  Pburg_20_16(:,j) = constructPowerSpectrum(aburg_16);
  Pburg_20_32(:,j) = constructPowerSpectrum(aburg_32);
end

// Convert to decibels.
Pburg_20_4 = 20*log10(Pburg_20_4);
Pburg_20_8 = 20*log10(Pburg_20_8);
Pburg_20_16 = 20*log10(Pburg_20_16);
Pburg_20_32 = 20*log10(Pburg_20_32);

//----------------------------------------------------------------
// Use modified covariance method.
//----------------------------------------------------------------
// Preallocate matrices.
Pmcov_20_4 = zeros(1024,1);
Pmcov_20_8 = zeros(1024,1);
Pmcov_20_16 = zeros(1024,1);
Pmcov_20_32 = zeros(1024,1);

for j = 1:20
  // Generate random phases.
  phi1 = rand() * 2*%pi;
  phi2 = rand() * 2*%pi;

  // Generate white noise sequence.
  v = feval([1:256],Noise);
  v = v(:);

  // Generate MA(4) process.
  w = filterBlock(v,b,0);

  // Construct random process.
  xmcov = 2*cos(w1*n_256 + phi1) + 2*cos(w2*n_256 + phi2) + w;

  // Generate spectral estimate using the Burg method.
  [amcov_4,e] = mcov(xmcov,4);
  [amcov_8,e] = mcov(xmcov,8);
  [amcov_16,e] = mcov(xmcov,16);
  [amcov_32,e] = mcov(xmcov,32);

  // Construct spectral estimates.
  Pmcov_20_4(:,j) = constructPowerSpectrum(amcov_4);
  Pmcov_20_8(:,j) = constructPowerSpectrum(amcov_8);
  Pmcov_20_16(:,j) = constructPowerSpectrum(amcov_16);
  Pmcov_20_32(:,j) = constructPowerSpectrum(amcov_32);
end

// Convert to decibels.
Pmcov_20_4 = 20*log10(Pmcov_20_4);
Pmcov_20_8 = 20*log10(Pmcov_20_8);
Pmcov_20_16 = 20*log10(Pmcov_20_16);
Pmcov_20_32 = 20*log10(Pmcov_20_32);

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Plot results.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

scf(1);

// Part (a), maximum entropy method.
subplot(4,3,1);
title('Px Using Mem, p:4');
plot(Pmem_20_4(:,1));

subplot(4,3,4);
title('Px Using Mem, p:8');
plot(Pmem_20_8(:,1));

subplot(4,3,7);
title('Px Using Mem, p:16');
plot(Pmem_20_16(:,1));

subplot(4,3,10);
title('Px Using Mem, p:32');
plot(Pmem_20_32(:,1));

subplot(4,3,2);
title('Px Using Mem - Overlay, p:4');
plot(Pmem_20_4);

subplot(4,3,5);
title('Px Using Mem - Overlay, p:8');
plot(Pmem_20_8);

subplot(4,3,8);
title('Px Using Mem - Overlay, p:16');
plot(Pmem_20_16);

subplot(4,3,11);
title('Px Using Mem - Overlay, p:32');
plot(Pmem_20_32);

subplot(4,3,3);
title('Px Avg Using Mem p:4');
plot(sum(Pmem_20_4,'c') / 20);

subplot(4,3,6);
title('Px Avg Using Mem p:8');
plot(sum(Pmem_20_8,'c') / 20);

subplot(4,3,9);
title('Px Avg Using Mem p:16');
plot(sum(Pmem_20_16,'c') / 20);

subplot(4,3,12);
title('Px Avg Using Mem p:32');
plot(sum(Pmem_20_32,'c') / 20);


scf(2);

// Part (b), Burg method.
subplot(4,3,1);
title('Px Using Burg, p:4');
plot(Pburg_20_4(:,1));

subplot(4,3,4);
title('Px Using Burg, p:8');
plot(Pburg_20_8(:,1));

subplot(4,3,7);
title('Px Using Burg, p:16');
plot(Pburg_20_16(:,1));

subplot(4,3,10);
title('Px Using Burg, p:32');
plot(Pburg_20_32(:,1));

subplot(4,3,2);
title('Px Using Burg - Overlay, p:4');
plot(Pburg_20_4);

subplot(4,3,5);
title('Px Using Burg - Overlay, p:8');
plot(Pburg_20_8);

subplot(4,3,8);
title('Px Using Burg - Overlay, p:16');
plot(Pburg_20_16);

subplot(4,3,11);
title('Px Using Burg - Overlay, p:32');
plot(Pburg_20_32);

subplot(4,3,3);
title('Px Avg Using Burg p:4');
plot(sum(Pburg_20_4,'c') / 20);

subplot(4,3,6);
title('Px Avg Using Burg p:8');
plot(sum(Pburg_20_8,'c') / 20);

subplot(4,3,9);
title('Px Avg Using Burg p:16');
plot(sum(Pburg_20_16,'c') / 20);

subplot(4,3,12);
title('Px Avg Using Burg p:32');
plot(sum(Pburg_20_32,'c') / 20);

scf(3);

// Part (b), modified covariance method.
subplot(4,3,1);
title('Px Using Mcov, p:4');
plot(Pmcov_20_4(:,1));

subplot(4,3,4);
title('Px Using Mcov, p:8');
plot(Pmcov_20_8(:,1));

subplot(4,3,7);
title('Px Using Mcov, p:16');
plot(Pmcov_20_16(:,1));

subplot(4,3,10);
title('Px Using Mcov, p:32');
plot(Pmcov_20_32(:,1));

subplot(4,3,2);
title('Px Using Mcov - Overlay, p:4');
plot(Pmcov_20_4);

subplot(4,3,5);
title('Px Using Mcov - Overlay, p:8');
plot(Pmcov_20_8);

subplot(4,3,8);
title('Px Using Mcov - Overlay, p:16');
plot(Pmcov_20_16);

subplot(4,3,11);
title('Px Using Mcov - Overlay, p:32');
plot(Pmcov_20_32);

subplot(4,3,3);
title('Px Avg Using Mcov p:4');
plot(sum(Pmcov_20_4,'c') / 20);

subplot(4,3,6);
title('Px Avg Using Mcov p:8');
plot(sum(Pmcov_20_8,'c') / 20);

subplot(4,3,9);
title('Px Avg Using Mcov p:16');
plot(sum(Pmcov_20_16,'c') / 20);

subplot(4,3,12);
title('Px Avg Using Mcov p:32');
plot(sum(Pmcov_20_32,'c') / 20);

