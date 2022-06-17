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


//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (c): No code needed.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (d): 
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/


//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Plot results.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

scf(1);

// Part (a).
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

