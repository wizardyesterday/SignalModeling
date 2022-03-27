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
// The noise signal, measured from the secondary sensor is,
// v2(n) = 0.8*v2(n) + g(n).
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

// Generate the noise signal that is measured by the secondary sensor.
v2 = filterBlock(g,1,-0.8);

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (c): Using filters of order p = 2,4, and 6, design and
// implement the Wiener noise cancellation filters.  Make plots
// of the estimated process ghat(n), and compare the average
// squared errors for each filter.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Generate corrected signal.
[dhat2,ghat2,e2] = WienerNoiseCanceller(x,g,v2,2);
[dhat4,ghat4,e4] = WienerNoiseCanceller(x,g,v2,4);
[dhat6,ghat6,e6] = WienerNoiseCanceller(x,g,v2,6);

//*********************************************
// Plot results.
//*********************************************
// Select window 1.
scf(1);

subplot(421);
title('Desired Signal, d(n)');
plot(d);

subplot(422);
title('Desired Signal Corrupted by Noise, x(n) = d(n) + g(n)');
plot(x);

subplot(423);
title('Noise Estimate with Wiener Filter of Order 2');
plot(ghat2);

subplot(425);
title('Noise Estimate with Wiener Filter of Order 4');
plot(ghat4);

subplot(427);
title('Noise Estimate with Wiener Filter of Order 6');
plot(ghat6);


subplot(424);
title('Signal Estimate with Wiener Filter of Order 2');
plot(dhat2);

subplot(426);
title('Signal Estimate with Wiener Filter of Order 4');
plot(dhat4);

subplot(428);
title('Signal Estimate with Wiener Filter of Order 6');
plot(dhat6);
//*********************************************

//*********************************************
// Let's see how things look with a desired signal
// of the form: x(n) = sin(n*w0).
//*********************************************
// Construct the desired signal.
w0 = 0.05 * %pi;
d1 = sin(n*w0);
d1 = d1(:);

// Construct measured signal.
x1 = d1 + g;

// Generate corrected signal.
[dhat2_1,ghat2_1,e2_] = WienerNoiseCanceller(x1,g,v2,2);
[dhat4_1,ghat4_1,e4_1] = WienerNoiseCanceller(x1,g,v2,4);
[dhat6_1,ghat6_1,e6_1] = WienerNoiseCanceller(x1,g,v2,6);

//*********************************************
// Plot results.
//*********************************************
// Select window 2.
scf(2);

subplot(511);
title('Desired Signal, d(n)');
plot(d1);

subplot(512);
title('Desired Signal Corrupted by Noise, x(n) = d(n) + g(n)');
plot(x1);

subplot(513);
title('Signal Estimate with Wiener Filter of Order 2');
plot(dhat2_1);

subplot(514);
title('Signal Estimate with Wiener Filter of Order 4');
plot(dhat4_1);

subplot(515);
title('Signal Estimate with Wiener Filter of Order 6');
plot(dhat6_1);

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (d): For this part, the desired signal has leaked into
// the secondary sensor.  The input to the Wiener filter is now,
// v0(n) = v2(n) + alpha*d(n).  The performance is to be evaluated
// for alpha = {0.1,0.5} for filter orders of p = 2, 4, and
// 6.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Compute the three perturbed noise measurements.
v01 = v2 + 0.1*d;
v02 = v2 + 0.7*d;

// Alpha = 0.1.
[dhat2_d1,ghat2_d1,e2_d1] = WienerNoiseCanceller(x1,g,v01,2);
[dhat4_d1,ghat4_d1,e4_d1] = WienerNoiseCanceller(x1,g,v01,4);
[dhat6_d1,ghat6_d1,e6_d1] = WienerNoiseCanceller(x1,g,v01,6);

// Alpha = 0.7.
[dhat2_d2,ghat2_d2,e2_d2] = WienerNoiseCanceller(x1,g,v02,2);
[dhat4_d2,ghat4_d2,e4_d2] = WienerNoiseCanceller(x1,g,v02,4);
[dhat6_d2,ghat6_d2,e6_d2] = WienerNoiseCanceller(x1,g,v02,6);
//*********************************************

// Plot results.
//*********************************************
// Select window 3.
scf(3);

subplot(321);
title('Tainted Noise, Alpha = 0.1, with Wiener Filter of Order 2');
plot(dhat2_d1);

subplot(322);
title('Tainted Noise, Alpha = 0.1, with Wiener Filter of Order 4');
plot(dhat4_d1);

subplot(323);
title('Tainted Noise, Alpha = 0.1, with Wiener Filter of Order 6');
plot(dhat6_d1);

subplot(324);
title('Tainted Noise, Alpha = 0.7, with Wiener Filter of Order 2');
plot(dhat2_d2);

subplot(325);
title('Tainted Noise, Alpha = 0.7, with Wiener Filter of Order 4');
plot(dhat4_d2);

subplot(326);
title('Tainted Noise, Alpha = 0.7, with Wiener Filter of Order 6');
plot(dhat6_d2);


