//**********************************************************************
// File Name: C9_1.sci
//**********************************************************************

// Bring in all of the utility functions.
exec('utils.sci',-1);

//**********************************************************************
// Mainline code.
//**********************************************************************
// Set number of realizations.
N = 100;

// Set variance.
sigmavSq = 0.25;

// Initialize white noise generator.
noisegen(1,500,sqrt(sigmavSq));

// Set AR(2) filter coefficients.
a1 = [1 .1 .8]';
a2 = [1 .1 -.8]';

// Set step sizes.
mu1 = 0.05
mu2 = 0.01;

//----------------------------------------
// Generate AR(2) processes.
//----------------------------------------
for j = 1:N
  // Generate the noise.
  v = feval([1:500],Noise);

  X1(:,j) = filterBlock(v,1,a1(2:$));
  X2(:,j) = filterBlock(v,1,a2(2:$));
end

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (a): Paper problem.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (b): Plot the learning curve of the adaptive filter for
// 100 different realizations of x(n).  Estimate the misadjustment
// by time-averaging over the final iterations of the ensemble-
// average learning curves.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Initialize running sums.
W1avg_1 = 0;
W1avg_2 = 0;

// Generate error sequences.
for j = 1:N
  // Generate x(n-1).
  xnm1 = filterBlock(X1(:,j),[0 1],0);

  [W1_1,E1_1(:,j)] = lms(xnm1,X1(:,j),mu1,2);
  [W1_2,E1_2(:,j)] = lms(xnm1,X1(:,j),mu2,2);

  // Update running sums.
  W1avg_1 = W1avg_1 + W1_1($,1:2);
  W1avg_2 = W1avg_2 + W1_2($,1:2);
end

// Compute squared error.
E1Sq_1 = E1_1 .* E1_1;
E1Sq_2 = E1_2 .* E1_2;

// Compute ensemble averages.  This is the learning curve.
for j = 1:500
  E1SqAvg_1(j) = mean(E1Sq_1(j,:)); 
  E1SqAvg_2(j) = mean(E1Sq_2(j,:)); 
end

// Compute steady-state error.
E1inf_1 = mean(E1SqAvg_1(401:500));
E1inf_2 = mean(E1SqAvg_2(401:500));

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (c): Estimate the steady-state values of the adaptive
// filter coefficients for step sizes of mu = 0.05 and 0.01.
// Note that running sums were computed in Part (b).
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Compute average frequency estimates.
W1avg_1 = W1avg_1 / N;
W1avg_2 = W1avg_2 / N;

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (d): Repeat the experiments in parts (b) and (c) with
// a(1) = 0.1, a(2) = -0.8, and sigmaV^2 = 0.25.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Initialize running sums.
W2avg_1 = 0;
W2avg_2 = 0;

// Generate error sequences.
for j = 1:N
  // Generate x(n-1).
  xnm1 = filterBlock(X2(:,j),[0 1],0);

  [W2_1,E2_1(:,j)] = lms(xnm1,X2(:,j),mu1,2);
  [W2_2,E2_2(:,j)] = lms(xnm1,X2(:,j),mu2,2);

  // Update running sums.
  W2avg_1 = W2avg_1 + W2_1($,1:2);
  W2avg_2 = W2avg_2 + W2_2($,1:2);
end

// Compute squared error.
E2Sq_1 = E2_1 .* E2_1;
E2Sq_2 = E2_2 .* E2_2;

// Compute ensemble averages.  This is the learning curve.
for j = 1:500
  E2SqAvg_1(j) = mean(E2Sq_1(j,:)); 
  E2SqAvg_2(j) = mean(E2Sq_2(j,:)); 
end

// Compute steady-state error.
E2inf_1 = mean(E2SqAvg_1(401:500));
E2inf_2 = mean(E2SqAvg_2(401:500));

// Compute average frequency estimates.
W2avg_1 = W2avg_1 / N;
W2avg_2 = W2avg_2 / N;

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (e): Repeat the above exercises using the sign-error,
// sign-data, and sign-sign algorithms.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

//-------------------------------------------------------------
// Sign-error.
//-------------------------------------------------------------
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// a(1) = 0.1, a(2) = -0.8, and sigmaV^2 = 0.25.
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Initialize running sums.
Wse1avg_1 = 0;
Wse1avg_2 = 0;

// Generate error sequences.
for j = 1:N
  // Generate x(n-1).
  xnm1 = filterBlock(X1(:,j),[0 1],0);

  [Wse1_1,Ese1_1(:,j)] = lms_signError(xnm1,X1(:,j),mu1,2);
  [Wse1_2,Ese1_2(:,j)] = lms_signError(xnm1,X1(:,j),mu2,2);

  // Update running sums.
  Wse1avg_1 = Wse1avg_1 + Wse1_1($,1:2);
  Wse1avg_2 = Wse1avg_2 + Wse1_2($,1:2);
end

// Compute squared error.
Ese1Sq_1 = Ese1_1 .* Ese1_1;
Ese1Sq_2 = Ese1_2 .* Ese1_2;

// Compute ensemble averages.  This is the learning curve.
for j = 1:500
  Ese1SqAvg_1(j) = mean(Ese1Sq_1(j,:)); 
  Ese1SqAvg_2(j) = mean(Ese1Sq_2(j,:)); 
end

// Compute steady-state error.
Ese1inf_1 = mean(Ese1SqAvg_1(401:500));
Ese1inf_2 = mean(Ese1SqAvg_2(401:500));

// Compute average frequency estimates.
Wse1avg_1 = Wse1avg_1 / N;
Wse1avg_2 = Wse1avg_2 / N;

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// a(1) = 0.1, a(2) = -0.8, and sigmaV^2 = 0.25.
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Initialize running sums.
Wse2avg_1 = 0;
Wse2avg_2 = 0;

// Generate error sequences.
for j = 1:N
  // Generate x(n-1).
  xnm1 = filterBlock(X2(:,j),[0 1],0);

  [Wse2_1,Ese2_1(:,j)] = lms_signError(xnm1,X2(:,j),mu1,2);
  [Wse2_2,Ese2_2(:,j)] = lms_signError(xnm1,X2(:,j),mu2,2);

  // Update running sums.
  Wse2avg_1 = Wse2avg_1 + Wse2_1($,1:2);
  Wse2avg_2 = Wse2avg_2 + Wse2_2($,1:2);
end

// Compute squared error.
Ese2Sq_1 = Ese2_1 .* Ese2_1;
Ese2Sq_2 = Ese2_2 .* Ese2_2;

// Compute ensemble averages.  This is the learning curve.
for j = 1:500
  Ese2SqAvg_1(j) = mean(Ese2Sq_1(j,:)); 
  Ese2SqAvg_2(j) = mean(Ese2Sq_2(j,:)); 
end

// Compute steady-state error.
Ese2inf_1 = mean(Ese2SqAvg_1(401:500));
Ese2inf_2 = mean(Ese2SqAvg_2(401:500));

// Compute average frequency estimates.
Wse2avg_1 = Wse2avg_1 / N;
Wse2avg_2 = Wse2avg_2 / N;

//-------------------------------------------------------------
// Sign-data.
//-------------------------------------------------------------
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// a(1) = 0.1, a(2) = -0.8, and sigmaV^2 = 0.25.
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Initialize running sums.
Wsd1avg_1 = 0;
Wsd1avg_2 = 0;

// Generate error sequences.
for j = 1:N
  // Generate x(n-1).
  xnm1 = filterBlock(X1(:,j),[0 1],0);

  [Wsd1_1,Esd1_1(:,j)] = lms_signData(xnm1,X1(:,j),mu1,2);
  [Wsd1_2,Esd1_2(:,j)] = lms_signData(xnm1,X1(:,j),mu2,2);

  // Update running sums.
  Wsd1avg_1 = Wsd1avg_1 + Wsd1_1($,1:2);
  Wsd1avg_2 = Wsd1avg_2 + Wsd1_2($,1:2);
end

// Compute squared error.
Esd1Sq_1 = Esd1_1 .* Esd1_1;
Esd1Sq_2 = Esd1_2 .* Esd1_2;

// Compute ensemble averages.  This is the learning curve.
for j = 1:500
  Esd1SqAvg_1(j) = mean(Esd1Sq_1(j,:)); 
  Esd1SqAvg_2(j) = mean(Esd1Sq_2(j,:)); 
end

// Compute steady-state error.
Esd1inf_1 = mean(Esd1SqAvg_1(401:500));
Esd1inf_2 = mean(Esd1SqAvg_2(401:500));

// Compute average frequency estimates.
Wsd1avg_1 = Wsd1avg_1 / N;
Wsd1avg_2 = Wsd1avg_2 / N;

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// a(1) = 0.1, a(2) = -0.8, and sigmaV^2 = 0.25.
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Initialize running sums.
Wsd2avg_1 = 0;
Wsd2avg_2 = 0;

// Generate error sequences.
for j = 1:N
  // Generate x(n-1).
  xnm1 = filterBlock(X2(:,j),[0 1],0);

  [Wsd2_1,Esd2_1(:,j)] = lms_signData(xnm1,X2(:,j),mu1,2);
  [Wsd2_2,Esd2_2(:,j)] = lms_signData(xnm1,X2(:,j),mu2,2);

  // Update running sums.
  Wsd2avg_1 = Wsd2avg_1 + Wsd2_1($,1:2);
  Wsd2avg_2 = Wsd2avg_2 + Wsd2_2($,1:2);
end

// Compute squared error.
Esd2Sq_1 = Esd2_1 .* Esd2_1;
Esd2Sq_2 = Esd2_2 .* Esd2_2;

// Compute ensemble averages.  This is the learning curve.
for j = 1:500
  Esd2SqAvg_1(j) = mean(Esd2Sq_1(j,:)); 
  Esd2SqAvg_2(j) = mean(Esd2Sq_2(j,:)); 
end

// Compute steady-state error.
Esd2inf_1 = mean(Esd2SqAvg_1(401:500));
Esd2inf_2 = mean(Esd2SqAvg_2(401:500));

// Compute average frequency estimates.
Wsd2avg_1 = Wsd2avg_1 / N;
Wsd2avg_2 = Wsd2avg_2 / N;

//-------------------------------------------------------------
// Sign-sign.
//-------------------------------------------------------------
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// a(1) = 0.1, a(2) = -0.8, and sigmaV^2 = 0.25.
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Initialize running sums.
Wss1avg_1 = 0;
Wss1avg_2 = 0;

// Generate error sequences.
for j = 1:N
  // Generate x(n-1).
  xnm1 = filterBlock(X1(:,j),[0 1],0);

  [Wss1_1,Ess1_1(:,j)] = lms_signSign(xnm1,X1(:,j),mu1,2);
  [Wss1_2,Ess1_2(:,j)] = lms_signSign(xnm1,X1(:,j),mu2,2);

  // Update running sums.
  Wss1avg_1 = Wss1avg_1 + Wss1_1($,1:2);
  Wss1avg_2 = Wss1avg_2 + Wss1_2($,1:2);
end

// Compute squared error.
Ess1Sq_1 = Ess1_1 .* Ess1_1;
Ess1Sq_2 = Ess1_2 .* Ess1_2;

// Compute ensemble averages.  This is the learning curve.
for j = 1:500
  Ess1SqAvg_1(j) = mean(Ess1Sq_1(j,:)); 
  Ess1SqAvg_2(j) = mean(Ess1Sq_2(j,:)); 
end

// Compute steady-state error.
Ess1inf_1 = mean(Ess1SqAvg_1(401:500));
Ess1inf_2 = mean(Ess1SqAvg_2(401:500));

// Compute average frequency estimates.
Wss1avg_1 = Wss1avg_1 / N;
Wss1avg_2 = Wss1avg_2 / N;

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// a(1) = 0.1, a(2) = -0.8, and sigmaV^2 = 0.25.
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Initialize running sums.
Wss2avg_1 = 0;
Wss2avg_2 = 0;

// Generate error sequences.
for j = 1:N
  // Generate x(n-1).
  xnm1 = filterBlock(X2(:,j),[0 1],0);

  [Wss2_1,Ess2_1(:,j)] = lms_signSign(xnm1,X2(:,j),mu1,2);
  [Wss2_2,Ess2_2(:,j)] = lms_signSign(xnm1,X2(:,j),mu2,2);

  // Update running sums.
  Wss2avg_1 = Wss2avg_1 + Wss2_1($,1:2);
  Wss2avg_2 = Wss2avg_2 + Wsd2_2($,1:2);
end

// Compute squared error.
Ess2Sq_1 = Ess2_1 .* Ess2_1;
Ess2Sq_2 = Ess2_2 .* Ess2_2;

// Compute ensemble averages.  This is the learning curve.
for j = 1:500
  Ess2SqAvg_1(j) = mean(Ess2Sq_1(j,:)); 
  Ess2SqAvg_2(j) = mean(Ess2Sq_2(j,:)); 
end

// Compute steady-state error.
Ess2inf_1 = mean(Ess2SqAvg_1(401:500));
Ess2inf_2 = mean(Ess2SqAvg_2(401:500));

// Compute average frequency estimates.
Wss2avg_1 = Wss2avg_1 / N;
Wss2avg_2 = Wss2avg_2 / N;

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Plot results.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

//------------------------------------------------------
// Parts (b) and (d).
//------------------------------------------------------
scf(1);

subplot(421);
title('LMS Coefficients, a = [1 0.1 0.8], mu = 0.05');
plot(W1_1);

subplot(422);
title('LMS Coefficients, a = [1 0.1 0.8], mu = 0.01');
plot(W1_2);

subplot(423);
title('Mean-square error, a = [1 0.1 0.8], mu = 0.05');
plot(E1Sq_1(:,N));

subplot(424);
title('Mean-square error, a = [1 0.1 0.8], mu = 0.01');
plot(E1Sq_2(:,N));

subplot(425);
title('LMS Coefficients, a = [1 0.1 -0.8], mu = 0.05');
plot(W2_1);

subplot(426);
title('LMS Coefficients, a = [1 0.1 -0.8], mu = 0.01');
plot(W2_2);

subplot(427);
title('Mean-square error, a = [1 0.1 -0.8], mu = 0.05');
plot(E2Sq_1(:,N));

subplot(428);
title('Mean-square error, a = [1 0.1 -0.8], mu = 0.01');
plot(E1Sq_2(:,N));

//------------------------------------------------------
// Part (e), Sign-Error.
//------------------------------------------------------
scf(2);

subplot(421);
title('LMS Coefficients, Sign-Error, a = [1 0.1 0.8], mu = 0.05');
plot(Wse1_1);

subplot(422);
title('LMS Coefficients, Sign-Error, a = [1 0.1 0.8], mu = 0.01');
plot(Wse1_2);

subplot(423);
title('Mean-square error, Sign-Error, a = [1 0.1 0.8], mu = 0.05');
plot(Ese1Sq_1(:,N));

subplot(424);
title('Mean-square error, Sign-Error, a = [1 0.1 0.8], mu = 0.01');
plot(Ese1Sq_2(:,N));

subplot(425);
title('LMS Coefficients, Sign-Error, a = [1 0.1 -0.8], mu = 0.05');
plot(Wse2_1);

subplot(426);
title('LMS Coefficients, Sign-Error, a = [1 0.1 -0.8], mu = 0.01');
plot(Wse2_2);

subplot(427);
title('Mean-square error, Sign-Error, a = [1 0.1 -0.8], mu = 0.05');
plot(Ese2Sq_1(:,N));

subplot(428);
title('Mean-square error, Sign-Error a = [1 0.1 -0.8], mu = 0.01');
plot(Ese2Sq_2(:,N));

//------------------------------------------------------
// Part (e), Sign-Data.
//------------------------------------------------------
scf(3);

subplot(421);
title('LMS Coefficients, Sign-Data, a = [1 0.1 0.8], mu = 0.05');
plot(Wsd1_1);

subplot(422);
title('LMS Coefficients, Sign-Data, a = [1 0.1 0.8], mu = 0.01');
plot(Wsd1_2);

subplot(423);
title('Mean-square error, Sign-Data, a = [1 0.1 0.8], mu = 0.05');
plot(Esd1Sq_1(:,N));

subplot(424);
title('Mean-square error, Sign-Data, a = [1 0.1 0.8], mu = 0.01');
plot(Esd1Sq_2(:,N));

subplot(425);
title('LMS Coefficients, Sign-Data, a = [1 0.1 -0.8], mu = 0.05');
plot(Wsd2_1);

subplot(426);
title('LMS Coefficients, Sign-Data, a = [1 0.1 -0.8], mu = 0.01');
plot(Wsd2_2);

subplot(427);
title('Mean-square error, Sign-Data, a = [1 0.1 -0.8], mu = 0.05');
plot(Esd2Sq_1(:,N));

subplot(428);
title('Mean-square error, Sign-Data, a = [1 0.1 -0.8], mu = 0.01');
plot(Esd2Sq_2(:,N));

//------------------------------------------------------
// Part (e), Sign-Sign.
//------------------------------------------------------
scf(4);

subplot(421);
title('LMS Coefficients, Sign-Sign, a = [1 0.1 0.8], mu = 0.05');
plot(Wss1_1);

subplot(422);
title('LMS Coefficients, Sign-Sign, a = [1 0.1 0.8], mu = 0.01');
plot(Wss1_2);

subplot(423);
title('Mean-square error, Sign-Sign, a = [1 0.1 0.8], mu = 0.05');
plot(Ess1Sq_1(:,N));

subplot(424);
title('Mean-square error, Sign-Sign, a = [1 0.1 0.8], mu = 0.01');
plot(Ess1Sq_2(:,N));

subplot(425);
title('LMS Coefficients, Sign-Sign, a = [1 0.1 -0.8], mu = 0.05');
plot(Wss2_1);

subplot(426);
title('LMS Coefficients, Sign-Sign, a = [1 0.1 -0.8], mu = 0.01');
plot(Wss2_2);

subplot(427);
title('Mean-square error, Sign-Sign, a = [1 0.1 -0.8], mu = 0.05');
plot(Ess2Sq_1(:,N));

subplot(428);
title('Mean-square error, Sign-Sign, a = [1 0.1 -0.8], mu = 0.01');
plot(Ess2Sq_2(:,N));

//------------------------------------------------------
// Plot learning curves.
//------------------------------------------------------
scf(5)

subplot(421)
title('Learning curve, LMS, a = [1 0.1 0.8], mu = 0.05');
plot(E1SqAvg_1);

subplot(422)
title('Learning curve, LMS, a = [1 0.1 0.8], mu = 0.01'); 
plot(E1SqAvg_2);

subplot(423)
title('Learning curve, LMS, a = [1 0.1 -0.8], mu = 0.05'); 
plot(E2SqAvg_1);

subplot(424)
title('Learning curve, LMS, a = [1 0.1 -0.8], mu = 0.01');
plot(E2SqAvg_2);

subplot(425);
title('Learning curve, Sign-Error, a = [1 0.1 0.8], mu = 0.05'); 
plot(Ese1SqAvg_1);

subplot(426);
title('Learning curve, Sign-Error, a = [1 0.1 0.8], mu = 0.01'); 
plot(Ese1SqAvg_2);

subplot(427);
title('Learning curve, Sign-Error, a = [1 0.1 -0.8], mu = 0.05'); 
plot(Ese2SqAvg_1);

subplot(428);
title('Learning curve, Sign-Error, a = [1 0.1 -0.8], mu = 0.01');
plot(Ese2SqAvg_2);

scf(6);

subplot(421)
title('Learning curve, Sign-Data, a = [1 0.1 0.8], mu = 0.05'); 
plot(Esd1SqAvg_1);

subplot(422)
title('Learning curve, Sign-Data, a = [1 0.1 0.8], mu = 0.01'); 
plot(Esd1SqAvg_2);

subplot(423)
title('Learning curve, Sign-Data, a = [1 0.1 -0.8], mu = 0.05'); 
plot(Esd2SqAvg_1);

subplot(424)
title('Learning curve, Sign-Data, a = [1 0.1 -0.8], mu = 0.01');
plot(Esd2SqAvg_2);

subplot(425);
title('Learning curve, Sign-Sign, a = [1 0.1 0.8], mu = 0.05'); 
plot(Ess1SqAvg_1);

subplot(426);
title('Learning curve, Sign-Sign, a = [1 0.1 0.8], mu = 0.01'); 
plot(Ess1SqAvg_2);

subplot(427);
title('Learning curve, Sign-Sign, a = [1 0.1 -0.8], mu = 0.05'); 
plot(Ess2SqAvg_1);

subplot(428);
title('Learning curve, Sign-Sign, a = [1 0.1 -0.8], mu = 0.01');
plot(Ess2SqAvg_2);

