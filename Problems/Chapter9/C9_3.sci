//**********************************************************************
// File Name: C9_3.sci
//**********************************************************************

// Bring in all of the utility functions.
exec('utils.sci',-1);

//**********************************************************************
// Mainline code.
//**********************************************************************
// Set number of realizations.
N = 100;

// Set variance.
sigmavSq = 1;

// Set scale factors.
Gamma1 = 0.1;
Gamma2 = 1;

//------------------------------------------
// Generate Gaussian random process.
x = generateGaussianProcess(1,1000,sigmavSq);

// Set MA(2) filter coefficients.
g = [1 1.8 0.8]';

// Set step size.
Beta = 0.1;

// Generate measurement noise.
vmeas = generateGaussianProcess(1,1000,sigmavSq);

// Generate reference signal, d(n).
d = filterBlock(x,g,0);

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (a): Implement adaptive filters of order, p = 4 and
// p = 5, using the normalized LMS algorithm.  Use a step size of
// Beta = 0.1.  The reference is d(n) + 0.1vmeas(n).
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Add noise to reference.
d1 = d + (Gamma1 * vmeas);

// Run normalized LMS adaptive filters.
[W_nlms_gamma1_p4,E_nlms_gamma1_p4] = nlms(x,d1,0.1,4);
[W_nlms_gamma1_p5,E_nlms_gamma1_p5] = nlms(x,d1,0.1,5);

// Output result.
printf("\nW_nlms_gamma1_p4, 0.1vmeasnoise, p = 4");
disp(W_nlms_gamma1_p4($,1:4));
printf("\nW_nlms_gamma1_p4, 0.1vmeasnoise, p = 5");
disp(W_nlms_gamma1_p5($,1:5));

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (b): Implement adaptive filters of order, p = 4 and
// p = 5, using the normalized LMS algorithm.  Use a step size of
// Beta = 0.1.  The reference is d(n) + vmeas(n).
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Add noise to reference.
d2 = d + (Gamma2 * vmeas);

// Run normalized LMS adaptive filters.
[W_nlms_gamma2_p4,E_nlms_gamma1_p4] = nlms(x,d2,0.1,4);
[W_nlms_gamma2_p5,E_nlms_gamma1_p5] = nlms(x,d2,0.1,5);

// Output result.
printf("\nW_nlms_gamma2_p4, vmeasnoise, p = 4");
disp(W_nlms_gamma2_p4($,1:4));
printf("\nW_nlms_gamma2_p4, vmeasnoise, p = 5");
disp(W_nlms_gamma2_p5($,1:5));

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Plot results.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

subplot(411);
title('Steady-state Coefficients (Noisy Reference), Gamma = 0.1, P = 4');
plot(W_nlms_gamma1_p4);

subplot(412);
title('Steady-state Coefficients (Noisy Reference), Gamma = 0.1, P = 5');
plot(W_nlms_gamma1_p5);

subplot(413);
title('Steady-state Coefficients (Noisy Reference), Gamma = 1, P = 4');
plot(W_nlms_gamma2_p4);

subplot(414);
title('Steady-state Coefficients (Noisy Reference), Gamma = 1, P = 5');
plot(W_nlms_gamma2_p5);
