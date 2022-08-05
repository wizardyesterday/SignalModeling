//**********************************************************************
// File Name: C9_4.sci
//**********************************************************************

// Bring in all of the utility functions.
exec('utils.sci',-1);

//**********************************************************************
// Mainline code.
//**********************************************************************
// Set number of realizations.
N = 30;

// Set number of samples.
numberOfSamples = 500;

// Set variance.
sigmavSq = 1;

//------------------------------------------
// Generate Gaussian random process.  We
// want 30 realizations of length 5=0.
// The number of realizations may be scaled
// back due to memory limitations.
//------------------------------------------
X = generateGaussianProcess(N,numberOfSamples,sigmavSq);

//-------------------------------------------------------------
// Generate reference noise.  The noise
// sources have variances of 0.01, 0.1, and
// 1, respectively
//-------------------------------------------------------------
v1 = generateGaussianProcess(1,numberOfSamples,sqrt(0.01));
v2 = generateGaussianProcess(1,numberOfSamples,sqrt(0.1));
v3 = generateGaussianProcess(1,numberOfSamples,1);

// Coefficients for g1.
b1 = [1 0.5];
a1 = [1 0.9];

// Coefficients for g2.
b2 = [1 0.5];
a2 = [1 0.2];

// Open output file.
fd = mopen('C9_4_output.txt','w');

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (a): Use the normalized LMS algorithm, with p = 4, and
// Beta = 0.1 to model G1(z).  Record the values of the filter
// coefficients, and make a plot of the learning curve as
// described in Exercise C9.2.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (b): Repeat part (a) for different values of p.  Let's
// choose p = 6, p = 8, p = 10.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

// Generate error sequences.
for j = 1:N
  // Generate reference signal, d(n).
  d = filterBlock(X(:,j),b1,a1(2:$));

  // Run normalized LMS adaptive filter.
  [W_g1_p4,E1_g1_p4(:,j)] = nlms(X(:,j),d,0.1,4);

  // Run normalized LMS adaptive filter.
  [W_g1_p6,E1_g1_p6(:,j)] = nlms(X(:,j),d,0.1,6);

  // Run normalized LMS adaptive filter.
  [W_g1_p8,E1_g1_p8(:,j)] = nlms(X(:,j),d,0.1,8);

  // Run normalized LMS adaptive filter.
  [W_g1_p10,E1_g1_p10(:,j)] = nlms(X(:,j),d,0.1,10);
end

// Compute squared errors.
ESq_g1_p4 = E1_g1_p4 .* E1_g1_p4;
ESq_g1_p6 = E1_g1_p6 .* E1_g1_p6;
ESq_g1_p8 = E1_g1_p8 .* E1_g1_p8;
ESq_g1_p10 = E1_g1_p10 .* E1_g1_p10;

// Compute ensemble averages.  This is the learning curve.
for j = 1:numberOfSamples
  ESqAvg_g1_p4(j) = mean(ESq_g1_p4(j,:)); 
  ESqAvg_g1_p6(j) = mean(ESq_g1_p6(j,:)); 
  ESqAvg_g1_p8(j) = mean(ESq_g1_p8(j,:)); 
  ESqAvg_g1_p10(j) = mean(ESq_g1_p10(j,:)); 
end

// Output results.
mfprintf(fd,"\nW, g = g1, p = 4\n");
mfprintf(fd,"%f\n",W_g1_p4($,1:$)');

mfprintf(fd,"\nW, g = g1, p = 6\n");
mfprintf(fd,"%f\n",W_g1_p6($,1:$)');

mfprintf(fd,"\nW, g = g1, p = 8\n");
mfprintf(fd,"%f\n",W_g1_p8($,1:$)');

mfprintf(fd,"\nW, g = g1, p = 10");
mfprintf(fd,"%f\n",W_g1_p10($,1:$)');

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (c): Repeat parts (a) and (b) when the denominator of
// G(z) is replaced by A(z) = 1 - 0.2z^(-1).
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

// Generate error sequences.
for j = 1:N
  // Generate reference signal, d(n).
  d = filterBlock(X(:,j),b2,a2(2:$));

  // Run normalized LMS adaptive filter.
  [W_g2_p4,E1_g2_p4(:,j)] = nlms(X(:,j),d,0.1,4);

  // Run normalized LMS adaptive filter.
  [W_g2_p6,E1_g2_p6(:,j)] = nlms(X(:,j),d,0.1,6);

  // Run normalized LMS adaptive filter.
  [W_g2_p8,E1_g2_p8(:,j)] = nlms(X(:,j),d,0.1,8);

  // Run normalized LMS adaptive filter.
  [W_g2_p10,E1_g2_p10(:,j)] = nlms(X(:,j),d,0.1,10);
end

// Compute squared errors.
ESq_g2_p4 = E1_g2_p4 .* E1_g2_p4;
ESq_g2_p6 = E1_g2_p6 .* E1_g2_p6;
ESq_g2_p8 = E1_g2_p8 .* E1_g2_p8;
ESq_g2_p10 = E1_g2_p10 .* E1_g2_p10;

// Compute ensemble averages.  These are the learning curve.
for j = 1:numberOfSamples
  ESqAvg_g2_p4(j) = mean(ESq_g2_p4(j,:)); 
  ESqAvg_g2_p6(j) = mean(ESq_g2_p6(j,:)); 
  ESqAvg_g2_p8(j) = mean(ESq_g2_p8(j,:)); 
  ESqAvg_g2_p10(j) = mean(ESq_g2_p10(j,:)); 
end

// Output results.
mfprintf(fd,"\nW, g = g2, p = 4\n");
mfprintf(fd,"%f\n",W_g2_p4($,1:$)');

mfprintf(fd,"\nW, g = g2, p = 6\n");
mfprintf(fd,"%f\n",W_g2_p6($,1:$)');

mfprintf(fd,"\nW, g = g2, p = 8\n");
mfprintf(fd,"%f\n",W_g2_p8($,1:$)');

mfprintf(fd,"\nW, g = g2, p = 10");
mfprintf(fd,"%f\n",W_g2_p10($,1:$)');

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (d): Add white noise, v(n), to the output of the
// unknown system.  Repeat part (a) with noise variances of
// sigmavSq = 0.01, 0.1, and 1.0.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

// Generate error sequences.
for j = 1:N
  // Generate reference signal, d(n).
  d = filterBlock(X(:,j),b1,a1(2:$));

  // Run normalized LMS adaptive filters wiht noisy reference.
  [W_g1_v1_p4,E1_g1_v1_p4(:,j)] = nlms(X(:,j),d+v1,0.1,4);
  [W_g1_v2_p4,E1_g1_v2_p4(:,j)] = nlms(X(:,j),d+v2,0.1,4);
  [W_g1_v3_p4,E1_g1_v3_p4(:,j)] = nlms(X(:,j),d+v3,0.1,4);

end

// Compute squared errors.
ESq_g1_v1_p4 = E1_g1_v1_p4 .* E1_g1_v1_p4;
ESq_g1_v2_p4 = E1_g1_v2_p4 .* E1_g1_v2_p4;
ESq_g1_v3_p4 = E1_g1_v3_p4 .* E1_g1_v3_p4;

// Compute ensemble averages.  This is the learning curve.
for j = 1:numberOfSamples
  ESqAvg_g1_v1_p4(j) = mean(ESq_g1_v1_p4(j,:)); 
  ESqAvg_g1_v2_p4(j) = mean(ESq_g1_v2_p4(j,:)); 
  ESqAvg_g1_v3_p4(j) = mean(ESq_g1_v3_p4(j,:)); 
end

// Output results.
mfprintf(fd,"\nW, g = g1, Reference noise = v1, p = 4\n");
mfprintf(fd,"%f\n",W_g1_v1_p4($,1:$)');

mfprintf(fd,"\nW, g = g1, Reference noise = v2, p = 4\n");
mfprintf(fd,"%f\n",W_g1_v2_p4($,1:$)');

mfprintf(fd,"\nW, g = g1, Reference noise = v3, p = 4\n");
mfprintf(fd,"%f\n",W_g1_v3_p4($,1:$)');

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (e): Write a MATLAB file to implement the filtered
// signal approach for IIR adaptive filtering.  With an adaptive
// recursive filter of the form,
//
// y(n) = a_n(1)y(n-1) + b_n(0)x(n) + b_n(1)x(n-1),
//
// determine an appropriate value for the step size mu, and use
// your m-file to model the system G(z).
// Note that mu = 0.0088 was found by experiment.  Additionally,
// 1500 samples are needed to illustrate convergence of the
// filter.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

// Generate a 1500-sample unit variance white noise process.
Xiir = generateGaussianProcess(1,1500,1);

// Generate reference signal, d(n).
d = filterBlock(Xiir,b1,a1(2:$));

// Run IIR adaptive filter.
[A_noNoise,B_noNoise,E_noNoise] = lms_iirFilteredSignal(Xiir,d,1,1,0.0088);

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Add white noise, v(n), to the output of the unknown system.
// Repeat part (e) with noise variances of sigmaV^2 = 0.01, 0.1,
// and 1.0.  For stability with noise, we set mu = 0.003.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

v1_iir = generateGaussianProcess(1,1500,sqrt(0.01));
v2_iir = generateGaussianProcess(1,1500,sqrt(0.1));
v3_iir = generateGaussianProcess(1,1500,1);

// Run IIR adaptive filter with different noise values.
[A_v1,B_v1,E_v1] = lms_iirFilteredSignal(Xiir,d+v1_iir,1,1,0.003);
[A_v2,B_v2,E_v2] = lms_iirFilteredSignal(Xiir,d+v2_iir,1,1,0.003);
[A_v3,B_v3,E_v3] = lms_iirFilteredSignal(Xiir,d+v3_iir,1,1,0.003);

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// We're done with the output file.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
mclose(fd);

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Plot results.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
scf(1);

// Part (a), part (b).
subplot(621);
title('FIR Learning curve, g2, p = 4');
plot(ESqAvg_g2_p4);

subplot(622);
title('FIR Learning curve, g2, p = 6');
plot(ESqAvg_g2_p6);

subplot(623);
title('Learning curve, g2, p = 8');
plot(ESqAvg_g2_p8);

subplot(624);
title('FIR Learning curve, g2, p = 10');
plot(ESqAvg_g2_p10);

// Part (c).
subplot(625);
title('FIR Learning curve, g2, p = 4');
plot(ESqAvg_g2_p4);

subplot(626);
title('FIR Learning curve, g2, p = 6');
plot(ESqAvg_g2_p6);

subplot(627);
title('FIR Learning curve, g2, p = 8');
plot(ESqAvg_g2_p8);

subplot(628);
title('FIR Learning curve, g2, p = 10');
plot(ESqAvg_g2_p10);

// Part (d).
subplot(629);
title('FIR Learning curve, g1, sigmaV^2 = 0.01, p = 4');
plot(ESqAvg_g1_v1_p4);

subplot(6,2,10);
title('FIR Learning curve, g1, sigmaV^2 = 0.1, p = 4');
plot(ESqAvg_g1_v2_p4);

subplot(6,2,11);
title('FIR Learning curve, g1, sigmaV^2 = 1, p = 4');
plot(ESqAvg_g1_v3_p4);

scf(2);

// Part (d), part (e).
subplot(621);
title('IIR Coefficients, a, mu = 0.0088, Noiseless');
plot(A_noNoise);

subplot(622);
title('IIR Coefficients, b, mu = 0.0088, Noiseless');
plot(B_noNoise);

subplot(623);
title('IIR Mean-Square Error, b, mu = 0.0088, Noiseless');
plot(E_noNoise);

subplot(624);
title('IIR Coefficients, a, mu = 0.003, sigmaV^2 = 0.01');
plot(A_v1);

subplot(625);
title('IIR Coefficients, b, mu = 0.003, signaV^2 = 0.01');
plot(B_v1);

subplot(626);
title('IIR Mean-Square Error, b, mu = 0.003, sigmaV^2 = 0.01');
plot(E_v1);

subplot(627);
title('IIR Coefficients, a, mu = 0.003, sigmaV^2 = 0.1');
plot(A_v2);

subplot(628);
title('IIR Coefficients, b, mu = 0.003, signaV^2 = 0.1');
plot(B_v2);

subplot(629);
title('IIR Mean-Square Error, b, mu = 0.003, sigmaV^2 = 0.1');
plot(E_v2);

subplot(6,2,10);
title('IIR Coefficients, a, mu = 0.003, sigmaV^2 = 1');
plot(A_v3);

subplot(6,2,11);
title('IIR Coefficients, b, mu = 0.003, signaV^2 = 1');
plot(B_v3);

subplot(6,2,12);
title('IIR Mean-Square Error, b, mu = 0.003, sigmaV^2 = 1');
plot(E_v3);

