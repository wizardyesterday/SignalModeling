//**********************************************************************
// File Name: C9_9.sci
//**********************************************************************

// Bring in all of the utility functions.
exec('utils.sci',-1);

//**********************************************************************
// Mainline code.
//**********************************************************************
// Open output file.
fd = mopen('C9_9_output.txt','w');

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (a): Paper exercise.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (b): Create NLMS noise canceller.  Given an input signal,
// x(n), the reference signal is x(n-k0), where k0 is the
// specified delay.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (c): Generate 1000 samples of v(n) and x(n), and use the
// NLMS noise canceller to estimate d(n) = sin(0.01PI*n), for
// filter orders p = 5, 10, 15, 20.  Use values for k0 that
// range from the minimum value determined in part (a) to
// k0 = 25.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Generate unit variance white Gaussian noise. random process.
g = generateGaussianProcess(1,1000,1);

// Coefficients for difference equation v(n) = g(n) + 0.5g(n-2).
b = [1 0 0.5];

// Generate MA(2) process.
v = filterBlock(g,b,0);

// Construct desired signal.
d = sin(0.01*%pi*[0:999]);

// Construct input signal.
x = d + v';

Beta = 0.1;

// Set filter orders.
p5 = 5;
p10 = 10;
p15 = 15;
p20 = 20;

// Set delays.
n5 = 5;
n10 = 10;
n25 = 25;

[W_n5_p5,dhat_n5_p5] = nlms_noiseCanceller(x,n5,Beta,p5);
[W_n5_p10,dhat_n5_p10] = nlms_noiseCanceller(x,n5,Beta,p10);
[W_n5_p15,dhat_n5_p15] = nlms_noiseCanceller(x,n5,Beta,p15);
[W_n5_p20,dhat_n5_p20] = nlms_noiseCanceller(x,n5,Beta,p20);

[W_n10_p5,dhat_n10_p5] = nlms_noiseCanceller(x,n10,Beta,p5);
[W_n10_p10,dhat_n10_p10] = nlms_noiseCanceller(x,n10,Beta,p10);
[W_n10_p15,dhat_n10_p15] = nlms_noiseCanceller(x,n10,Beta,p15);
[W_n10_p20,dhat_n10_p20] = nlms_noiseCanceller(x,n10,Beta,p20);

[W_n25_p5,dhat_n25_p5] = nlms_noiseCanceller(x,n25,Beta,p5);
[W_n25_p10,dhat_n25_p10] = nlms_noiseCanceller(x,n25,Beta,p10);
[W_n25_p15,dhat_n25_p15] = nlms_noiseCanceller(x,n25,Beta,p15);
[W_n25_p20,dhat_n25_p20] = nlms_noiseCanceller(x,n25,Beta,p20);

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (d):
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Compute cross-correlation sequence.
for n = 1:900
  rdx(n) = crosscorrelate(d,x,n-1);
end

// Compute the variance of v(n).
sigmavSq = variance(v);

// Run the signal through the Wiener filter.
[w,e] = FirWienerFilter(rdx,sigmavSq,10);

for i = 1:10
  mfprintf(fd,"Tap: %d ",i);
  mfprintf(fd,"Adaptive Filter: %f  Wiener Filter: %f\n",W_n25_p10($,i),w(i));
end

// Let's test the Wiener filter.
dhatWiener = filterBlock(x,w,0);

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// We're done with the output file.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
mclose(fd);

////_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Plot results.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
scf(1);

subplot(421);
title('Desired Signal, d(n)');
plot(d);

subplot(423);
title('Noise Signal, v(n)');
plot(v);

subplot(425);
title('Noise Corrupted Signal, x(n)');
plot(x);

subplot(427);
title('Wiener Filter Output, p = 10');
plot(dhatWiener);

// n0 = 5.

subplot(422);
title('Estimated Signal, dhat(n), n0 = 5, p = 5');
plot(dhat_n5_p5);

subplot(424);
title('Estimated Signal, dhat(n), n0 = 5, p = 10');
plot(dhat_n5_p10);

subplot(426);
title('Estimated Signal, dhat(n), n0 = 5, p = 15');
plot(dhat_n5_p15);

subplot(428);
title('Estimated Signal, dhat(n), n0 = 5, p = 20');
plot(dhat_n5_p20);

scf(2);

// n0 = 10.
subplot(421);
title('Estimated Signal, dhat(n), n0 = 10, p = 5');
plot(dhat_n10_p5);

subplot(423);
title('Estimated Signal, dhat(n), n0 = 10, p = 10');
plot(dhat_n10_p10);

subplot(425);
title('Estimated Signal, dhat(n), n0 = 10, p = 15');
plot(dhat_n10_p15);

subplot(427);
title('Estimated Signal, dhat(n), n0 = 10, p = 20');
plot(dhat_n10_p20);

// n0 = 25.
subplot(422);
title('Estimated Signal, dhat(n), n0 = 25, p = 5');
plot(dhat_n25_p5);

subplot(424);
title('Estimated Signal, dhat(n), n0 = 25, p = 10');
plot(dhat_n25_p10);

subplot(426);
title('Estimated Signal, dhat(n), n0 = 25, p = 15');
plot(dhat_n25_p15);

subplot(428);
title('Estimated Signal, dhat(n), n0 = 25, p = 20');
plot(dhat_n25_p20);



