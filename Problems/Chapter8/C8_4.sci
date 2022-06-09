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


// Generate processe.
x = filterBlock(w,[1 0 -1],0);

// Generate exact spectra.
Px = constructPowerSpectrum(x);

Pxa_2 = mem(x,2);
Pxa_4 = mem(x,4);
Pxa_40 = mem(x,40);

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (b): For sigmaW^2 = 0.5, 1, 2, 5, generate N = 100
// samples of the process y(n), and estimate the power spectrum
// of x(n) from y(n) using the maximum entropy method with p = 2.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

//*********************************************
// Plot results.
//*********************************************
// Select window 1.
//scf(1);

