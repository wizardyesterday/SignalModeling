//**********************************************************************
// File Name: C8_4.sci
//**********************************************************************

// Bring in all of the utility functions.
exec('utils.sci',-1);

//**********************************************************************
// Mainline code.
//**********************************************************************

// In this exercise, we look at what happens if an autoregressive
// spectrum estimation technique is used on a moving average process.

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (a): The process, x(n), is a second-order moving average
// process that is formed by filtering unit variance white
// Gaussian noise w(n) as follows: x(n) = w(n) - w(n-2).
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Generate white noise sequence with unit variance.
noisegen(1,256,1);
w = feval([1:256],Noise);
w = w(:);


// Generate process.
xa = filterBlock(w,[1 0 -1],0);

// Generate exact spectra.
Pxa = constructPowerSpectrum(xa);

// Perform spectral estimates.
Pxa_2 = mem(xa,2);
Pxa_8 = mem(xa,8);
Pxa_20 = mem(xa,20);

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (b): Repeat part (a) for the MA(3) process that is formed
// by filtering white Gaussian noise with the filter,
// H(z) = {1 - 0.98z^(-1)}{1 - 0.96z^(-2)}.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Form filter coefficients.
b = convol([1 -0.98],[1 0 0.96]);

// Generate process.
xb = filterBlock(w,b,0);

// Generate exact spectra.
Pxb = constructPowerSpectrum(xb);

// Perform spectral estimates.
Pxb_3 = mem(xb,3);
Pxb_8 = mem(xb,8);
Pxb_20 = mem(xb,20);

//*********************************************
// Plot results.
//*********************************************
// Select window 1.
//scf(1);

// Part (a).
subplot(421);
title('MA(2) Spectrum');
plot(Pxa);

subplot(423);
title('Mem Spectrum for MA(2) Process, p: 2');
plot(Pxa_2);

subplot(425);
title('Mem Spectrum for MA(2) Process, p: 8');
plot(Pxa_8);

subplot(427);
title('Mem Spectrum for MA(2) Process, p: 20');
plot(Pxa_20);

// Part (b).
subplot(422);
title('MA(3) Spectrum');
plot(Pxb);

subplot(424);
title('Mem Spectrum for MA(3) Process, p: 3');
plot(Pxb_3);

subplot(426);
title('Mem Spectrum for MA(3) Process, p: 8');
plot(Pxb_8);

subplot(428);
title('Mem Spectrum for MA(3) Process, p: 20');
plot(Pxb_20);
