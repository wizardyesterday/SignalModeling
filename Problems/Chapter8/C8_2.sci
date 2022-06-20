//**********************************************************************
// File Name: C8_2.sci
//**********************************************************************

// Bring in all of the utility functions.
exec('utils.sci',-1);

//**********************************************************************
// Mainline code.
//**********************************************************************

// The process of interest is created by driving a filter,
// H(z) = 1/{1 -  1.585*z^(-1) + 0.96z^(-2)}{1 -  1.152*z^(-1) + 0.96z^(-2)}
// with unit variance white Gaussian noise.

// Create the filter.
a1 = [1 -1.585 0.96];
a2 = [1 -1.152 0.96];
a = convol(a1,a2);
a = a(:);

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (a): Generate 256 samples of the process random process.
// Make a spectrum estimate using the autocorrelation method, and
// compare to the exact spectrum.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Generate white noise sequence with unit variance.
noisegen(1,256,1);
w = feval([1:256],Noise);

// Generate process.
x = filterBlock(w,1,a(2:$));

// Generate exact spectrum.
Pa = constructPowerSpectrum(1 ./ a);


// Generate spectral estimate using the autocorrelation method.
[ahatAcm,e] = acm(x,4);
PaAcm = constructPowerSpectrum(1 ./ ahatAcm);

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (b): Repeat part (a) for 20 different realizations of the
// process x(n).  Generate an overlay plot of the 20 estimates,
// and plot the ensemble average.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
for i = 1:20
  noisegen(1,256,1);
  w = feval([1:256],Noise);

  // Create a realization.
  x20 = filterBlock(w,1,a(2:$));

  // Generate spectral estimate using the autocorrelation method.
  [abhatAcm,e] = acm(x20,4);
  PbAcm(:,i) = constructPowerSpectrum(1 ./ abhatAcm);  
end

// Compute the average power spectrum.
PbAcmAvg = sum(PbAcm,'c') / 20;

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (c): Repeat part (b) for model orders of p = 6, 8, and
// 12. Describe what happens when the model order becomes too
// large.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
for i = 1:20
  noisegen(1,256,1);
  w = feval([1:256],Noise);

  // Create a realization.
  x20 = filterBlock(w,1,a(2:$));

  // Generate spectral estimate using the autocorrelation method.
  [achatAcm_6,e] = acm(x20,6);
  [achatAcm_8,e] = acm(x20,8);
  [achatAcm_12,e] = acm(x20,12);

  // Autocorrelation method spectrum.
  PcAcm_6(:,i) = constructPowerSpectrum(1 ./ achatAcm_6);  
  PcAcm_8(:,i) = constructPowerSpectrum(1 ./ achatAcm_8);  
  PcAcm_12(:,i) = constructPowerSpectrum(1 ./ achatAcm_12);  
end

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (d): Repeat parts (b) and (c) for the covariance,
// modified covariance, and Burg algorithms.  Which approach
// seems to be the best for a broadband autoregressive process?
// Note that 20 different realizations will be generated and
// estimates for p = 4, 6, 8, and 12 will be used.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
for i = 1:20
  noisegen(1,256,1);
  w = feval([1:256],Noise);

  // Create a realization.
  x20 = filterBlock(w,1,a(2:$));

  // Generate spectral estimate using the covariance method.
  [adhatCovm_4,e] = covm(x20,4);
  [adhatCovm_6,e] = covm(x20,6);
  [adhatCovm_8,e] = covm(x20,8);
  [adhatCovm_12,e] = covm(x20,12);

  // Generate spectral estimate using the modified covariance method.
  [adhatMcov_4,e] = mcov(x20,4);
  [adhatMcov_6,e] = mcov(x20,6);
  [adhatMcov_8,e] = mcov(x20,8);
  [adhatMcov_12,e] = mcov(x20,12);

  // Generate spectral estimate using the Burg method.
  [gamm,e] = burg(x,4);
  adhatBurg_4 = gtoa(gamm);
  [gamm,e] = burg(x,6);
  adhatBurg_6 = gtoa(gamm);
  [gamm,e] = burg(x,8);
  adhatBurg_8 = gtoa(gamm);
  [gamm,e] = burg(x,12);
  adhatBurg_12 = gtoa(gamm);

  // Covariance method spectrum.
  PdCovm_4(:,i) = constructPowerSpectrum(1 ./ adhatCovm_4);  
  PdCovm_6(:,i) = constructPowerSpectrum(1 ./ adhatCovm_6);  
  PdCovm_8(:,i) = constructPowerSpectrum(1 ./ adhatCovm_8);  
  PdCovm_12(:,i) = constructPowerSpectrum(1 ./ adhatCovm_12);

  // Modified covariance method spectrum.
  PdMcov_4(:,i) = constructPowerSpectrum(1 ./ adhatMcov_4);  
  PdMcov_6(:,i) = constructPowerSpectrum(1 ./ adhatMcov_6);  
  PdMcov_8(:,i) = constructPowerSpectrum(1 ./ adhatMcov_8);  
  PdMcov_12(:,i) = constructPowerSpectrum(1 ./ adhatMcov_12);  

  // Burg method spectrum. 
  PdBurg_4(:,i) = constructPowerSpectrum(1 ./ adhatBurg_4);  
  PdBurg_6(:,i) = constructPowerSpectrum(1 ./ adhatBurg_6);  
  PdBurg_8(:,i) = constructPowerSpectrum(1 ./ adhatBurg_8);  
  PdBurg_12(:,i) = constructPowerSpectrum(1 ./ adhatBurg_12);   
end

//*********************************************
// Plot results.
//*********************************************
// Select window 1.
scf(1);

// Part (a).
subplot(421);
title('Desired Signal');
plot(Pa);

subplot(422);
title('ACM Method, p: 4');
plot(PaAcm);

// Part (b).
subplot(423);
title('ACM Method, p: 4');
plot(PbAcm);

subplot(424);
title('ACM Method Average, p: 4');
plot(PbAcmAvg);

// Part (c).
subplot(425);
title('ACM Method, p: 6');
plot(PcAcm_6);

subplot(426);
title('ACM Method, p: 8');
plot(PcAcm_8);

subplot(427);
title('ACM Method, p: 12');
plot(PcAcm_12);

// Select window 2.
scf();

// Part (d).
subplot(431);
title('Covm Method, p:4');
plot(PdCovm_4);

subplot(434);
title('Covm Method, p:6');
plot(PdCovm_6);

subplot(437);
title('Covm Method, p:8');
plot(PdCovm_8);

subplot(4,3,10);
title('Covm Method, p:12');
plot(PdCovm_12);

subplot(432);
title('Mcov Method, P:4');
plot(PdMcov_4);

subplot(435);
title('Mcov Method, P:6');
plot(PdMcov_6);

subplot(438);
title('Mcov Method, P:8');
plot(PdMcov_8);

subplot(4,3,11);
title('Mcov Method, P:12');
plot(PdMcov_12);

subplot(433);
title('Burg Method, p:4');
plot(PdBurg_4);

subplot(436);
title('Burg Method, p:6');
plot(PdBurg_6);

subplot(439);
title('Burg Method, p:8');
plot(PdBurg_8);

subplot(4,3,12);
title('Burg Method, p:12');
plot(PdBurg_12);
