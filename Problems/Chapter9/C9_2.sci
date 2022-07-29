//**********************************************************************
// File Name: C9_2.sci
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

//------------------------------------------
// Generate Gaussian random process.  We
// want 100 realizations of length 1000.
// The number of realizations may be scaled
// back due to memory limitations.
//------------------------------------------
X = generateGaussianProcess(N,1000,sigmavSq);

// Set MA(2) filter coefficients.
g = [1 1.8 0.8]';

//+++++++++++++++++++++++++++++++++++++++++
// Generate autocorrelation sequence.
//+++++++++++++++++++++++++++++++++++++++++
rg = convol(g,g($:-1:1));
k = find(g == max(g));
rg = rg(k:$);

// Construct autocorrelation matrix.
Rg = toeplitz(rg);

// Compute eigenvalues.
lamda = spec(Rg);
lamdaMax = max(lamda);

// Set step size
muMax = 2/lamdaMax;
mu1 = 0.1 * muMax;
mu2 = 0.01 * muMax;
mu3 = 0.2 * muMax;

// Compute LMS misadjustments.
M1 = computeLmsMisadjustment(mu1,lamda);
M2 = computeLmsMisadjustment(mu2,lamda);
M3 = computeLmsMisadjustment(mu3,lamda);

// Generate reference signal, d(n).
d = filterBlock(X(:,1),g,0);

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (a): Paper problem.  The range of values for mu have
// already been computed.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (b): Implement an adaptive filter of order, p = 4, using
// the LMS algorithm.  Use a step size of mu = 0.1muMax.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Run LMS adaptive filter.
[W_lms_mu1,E_lms_mu1] = lms(X(:,1),d,mu1,4);

// Output result.
printf("\nW_lms_mu1, p = 4");
disp(W_lms_mu1($,1:4));

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (c): Repeat part (b) using the normalized LMS algorithm
// with beta = 0.1
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Run normalized LMS adaptive filter.
[W_nlms_mu1,E_nlms_mu1] = nlms(X(:,1),d,0.1,4);

// Output result.
printf("\nW_nlms_mu1, p = 4");
disp(W_nlms_mu1($,1:4));

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (d): Make a plot of the learning curve by repeating the
// experiment described in part (b) for 100 different realizations
// of d(n), and plotting the average of the plots of e^2(n) versus
// n.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Generate error sequences.
for j = 1:N
  // Generate reference signal, d(n).
  d = filterBlock(X(:,j),g,0);

  [W_dummy,E_lms_mu1(:,j)] = lms(X(:,j),d,mu1,4);
end

// Compute squared error.
ESq_lms_mu1 = E_lms_mu1 .* E_lms_mu1;

// Compute ensemble averages.  This is the learning curve.
for j = 1:1000
  ESqAvg_lms_mu1(j) = mean(ESq_lms_mu1(j,:)); 
end

// Compute steady-state error.
Einf_lms_mu1 = mean(ESqAvg_lms_mu1(901:1000));

// Output result.
printf("\nEinf_lms_mu1, p = 4");
disp(Einf_lms_mu1);

printf("\nTheoretical Excess Mean Square Value, mu1, p = 4");
disp(trace(Rg) * mu1 * sigmavSq / 2);

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (e): Repeat parts (b) and (d) for mu = 0.01muMax and
// mu = 0.2muMax.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Generate error sequences.
for j = 1:N
  // Generate reference signal, d(n).
  d = filterBlock(X(:,j),g,0);

  [W_dummy,E_lms_mu2(:,j)] = lms(X(:,j),d,mu2,4);
  [W_dummy,E_lms_mu3(:,j)] = lms(X(:,j),d,mu3,4);
end

// Compute squared errors.
ESq_lms_mu2 = E_lms_mu2 .* E_lms_mu2;
ESq_lms_mu3 = E_lms_mu3 .* E_lms_mu3;

// Compute ensemble averages.  This is the learning curve.
for j = 1:1000
  ESqAvg_lms_mu2(j) = mean(ESq_lms_mu2(j,:)); 
  ESqAvg_lms_mu3(j) = mean(ESq_lms_mu3(j,:)); 
end

// Compute steady-state errors.
Einf_lms_mu2 = mean(ESqAvg_lms_mu2(901:1000));
Einf_lms_mu3 = mean(ESqAvg_lms_mu3(901:1000));

// Output results.
printf("\nE1inf_lms_mu2, p = 4");
disp(Einf_lms_mu2);

printf("\nE1inf_lms_mu3, p = 4");
disp(Einf_lms_mu3);

printf("\nTheoretical Excess Mean Square Value, mu2, p = 4");
disp(trace(Rg) * mu2 * sigmavSq / 2);

printf("\nTheoretical Excess Mean Square Value, mu3, p = 4");
disp(trace(Rg) * mu3 * sigmavSq / 2);

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (f): Repeat parts (b) and (d) for p = 3 and p = 2.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Generate error sequences.
for j = 1:N
  // Generate reference signal, d(n).
  d = filterBlock(X(:,j),g,0);

  [W_dummy,E_lms_mu1_p3(:,j)] = lms(X(:,j),d,mu2,3);
  [W_dummy,E_lms_mu1_p2(:,j)] = lms(X(:,j),d,mu3,2);
end

// Compute squared errors.
ESq_lms_mu1_p3 = E_lms_mu1_p3 .* E_lms_mu1_p3;
ESq_lms_mu1_p2 = E_lms_mu1_p2 .* E_lms_mu1_p2;

// Compute ensemble averages.  This is the learning curve.
for j = 1:1000
  ESqAvg_lms_mu1_p3(j) = mean(ESq_lms_mu1_p3(j,:)); 
  ESqAvg_lms_mu1_p2(j) = mean(ESq_lms_mu1_p2(j,:)); 
end

// Compute steady-state errors.
Einf_lms_mu1_p3 = mean(ESqAvg_lms_mu1_p3(901:1000));
Einf_lms_mu1_p2 = mean(ESqAvg_lms_mu1_p2(901:1000));

// Output results.
printf("\nEinf_lms_mu1_p3, p = 3");
disp(Einf_lms_mu1_p3);

printf("\nEinf_lms_mu1_p2, p = 4");
disp(Einf_lms_mu1_p2);

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Plot results.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

// Part (d).
subplot(321);
title('Learning Curve, mu = 0.1muMax, p = 4');
plot(ESqAvg_lms_mu1);

// Part (e).
subplot(323);
title('Learning Curve mu = 0.01muMax, p = 4');
plot(ESqAvg_lms_mu2);

subplot(324);
title('Learning Curve mu = 0.2muMax, p = 4');
plot(ESqAvg_lms_mu3);

// Part (f).
subplot(325);
title('Learning Curve, mu = 0.1muMax, p = 3');
plot(ESqAvg_lms_mu1_p3);

subplot(326);
title('Learning Curve, mu = 0.1muMax, p = 2');
plot(ESq_lms_mu1_p2);


