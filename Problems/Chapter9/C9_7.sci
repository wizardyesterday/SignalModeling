//**********************************************************************
// File Name: C9_7.sci
//**********************************************************************

// Bring in all of the utility functions.
exec('utils.sci',-1);

//**********************************************************************

//  Name: crosscorrelate
//
//  Purpose: The purpose of this function is to compute the
//  estimated cross-correlation function of two input sequences.
//
//  Calling Sequence: r = crosscorrelate(x,lag)
//
//  Inputs:
//
//    x - The input vector to be processed.
//
//    y - The input vector to be cross-correlated with x.
//
//    numberOfSamples - The number of samples for which to carry out
//    the cross-correlation.
//
//    lag - The autocorrelation lag.
//
//  Outputs:
//
//    r - The estimated autocorrelation function.
//
//**********************************************************************
function r = crosscorrelate(x,y,numberOfSamples,lag);

  // Start with a clean slate.
  accumulator = 0;

  // Avoid constant reevaluation during loop execution.
  upperLimit = numberOfSamples - lag;

  for i = 1:upperLimit
    // Perform correlation processing.
    accumulator = accumulator + x(i) * y(i + lag);
  end

  // Scale as appropriate.
  r = accumulator / upperLimit;

endfunction

//**********************************************************************
//
//  Name: lms_pvector
//
//  The purpose of this function is to perform adaptive filtering using
//  the Widrow-Hoff LMS algorithm with a p-vector variant.
//
//  Calling Sequence: [W,E] = lms_pvector(x,Rdx,mu,nord,w0)
//
//  Inputs:
//
//    x - The input data to the adaptive filter.
//
//    Rdx - The cross-correlation between d(n) and x(n).  This is
//    a matrix of nord columns per row such that a cross-correlation
//    was constructed nord elements at a time.
//    Note: this parameter may have to be reworked.
//
//    mu - The adaptive filtering update (step-size) parameter.
//
//    nord - The number of filter coefficients.
//
//    w0 - An optional row vector that serves as the initial guess
//    for FIR filter coefficients.  If w0 is omitted, then w0 = 0 is
//    assumed.
//    
//  Outputs:
//
//    W - A matrix of filter coefficients as they evolve over time.
//    Each row of this matrix contains the coefficients at the
//    iteration that is associated with the row.  For example, row 1
//    contains the coefficients for iteration 1, row 2 contains the
//    coefficients for iteration 2, and so on.
//
//    E - A vector of errors as they evolve over time.  Each entry
//    contains the error that is associated with each iteration.  For
//    example, the first entry contains the error for iteration 1,
//    the second entry contains the error for iteration 2, and so on.
//    Note that E(n) = d(n) - dhat(n).
//
//**********************************************************************
function W = lms_pvector(x,Rdx,mu,nord,w0)

  // Construct the data matrix.
  X = convm(x,nord);

  // Retrieve the size of the data matrix.
  [M,N] = size(X);

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // argn(2) returns is the number of arguments passed to the 
  // function.
  // If 4 arguments were passed to the function, it is implied
  // that the last two parameter was not passed.  In this
  // case, the initial condition for the filter coefficients
  // is set to a default value of all zeros. 
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  if argn(2) < 5
    w0 = zeros(1,N);
  end

  // Force a row vector without altering the values.
  w0   = w0(:).';

  // Perform the first iteration for the filter vector.
  W(1,:) = w0 + mu * Rdx(1,:) - mu * (w0 * X(1,:).') * X(1,:);

  if M > 1
    for k = 2:M - nord + 1

      // Do this to simplify notation.      
      correction = mu * Rdx(k,:) - mu * (W(k-1,:) * X(k,:).') * X(k,:);

      // Update the filter vector.
      W(k,:) = W(k-1,:) + correction;
    end
  end

endfunction

//**********************************************************************
// Mainline code.
//**********************************************************************
// Set number of realizations.
N = 1;

// Set number of samples.
numberOfSamples = 500;

// Set variance.
sigmavSq = 1;

//------------------------------------------
// Generate Gaussian random process.  We
// want 30 realizations of length 500.
// The number of realizations may be scaled
// back due to memory limitations.
//------------------------------------------
X = generateGaussianProcess(N,numberOfSamples,sigmavSq);

// Coefficients for a1.
a1 = [1 0.8 -0.9]';

// Compute autocorrelation.
r = ator(a1,1);

// Construct autocorrelation matrix.
R = toeplitz(r);

// Compute eigenvalues.
lamda = spec(R);

// Compute maximum LMS step size.
muMax = 2 / max(lamda);

// Choise required step size.
mu1 = muMax / 800;

//----------------------------------------
// Generate AR(2) processes.
//----------------------------------------
for j = 1:N
  X1(:,j) = filterBlock(X(:,j),1,-a1(2:$));
end

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (a): Generate a sequence, x(n), of length N = 500, using
// the difference equation, x(n) = 0.8x(n-1) - 0.9x(n-2) + v(n),
// where v(n) is unit variance white Gaussian noise.  Determine
// a value for the LMS step size so that w_n converges in the
// mean to w = [0.8 -0.9]' within about 200 iterations.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (b): Implement the LMS adaptive predictor, and plot
// w_n(k) versus n for k = 1,2.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (c): Estimate the mean-square error, E(infinity), and
// compare it to the theoretical value.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

// Generate multiple realizations.
for j = 1:N
  // Generate x(n-1).
  xnm1 = filterBlock(X1(:,j),[0 1],0);

  // run LMS adaptive filter.
  [W_a1_mu1_p2,E_a1_mu1_p2(:,j)] = lms(xnm1,X1(:,j),mu1,2);
end

// Compute squared error.
ESq_a1_mu1_p2 = E_a1_mu1_p2 .* E_a1_mu1_p2;

// Compute ensemble averages.  This is the learning curve.
for j = 1:numberOfSamples
 ESqAvg_a1_mu1_p2(j) = mean(ESq_a1_mu1_p2(j,:)); 
end

// Estimate the steady-state mean-square error.
EInf_a1_mu1_p2 = mean(ESqAvg_a1_mu1_p2(numberOfSamples-99:numberOfSamples));

// Compute the theoretical steady-state mean-square error
EInfTheor_a1_mu1_p2 = sigmavSq / (1 - mu1 * trace(R) / 2);

mprintf("\nLMS E(infinity): %f\n",EInf_a1_mu1_p2);
mprintf("\nLMS Theoretical E(infinity): %f\n",EInfTheor_a1_mu1_p2);

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (d): Paper exercise.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (e): Find the vector, rdx = E{d(n)x(n)}, for a one-step
// linear predictor, and implement the p-vector algorithm for
// x(n).  Plot w_n(k) versus n, and compare your results to the
// LMS adaptive linear predictor.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Set filter order.
p = 2;

// Generate unit variance white Gaussian noise.
v = generateGaussianProcess(1,2 * numberOfSamples,1);

// Set desired signal, the AR(2) process.
d = filterBlock(v,1,-a1(2:$));

// Generate x(n-1).
xnm1 = filterBlock(d,[0 1],0);

//----------------------------------------------------
// Perform cross-correlation using p values at a time.
// This creates p columns in each row of the cross-
// correlation matrix.  The number of elements, in
// each row, is the same as the number of elements
// in the coefficient vector.
//----------------------------------------------------
for j = 1:numberOfSamples

  // Perform cross-correlation.
  for k = 1:p
    rdx(k) = crosscorrelate(d(j:$),xnm1(j:$),numberOfSamples,k-1);
  end

  // Set the current row to the cross-correlation.
  Rdx(j,:) = rdx';
end
//----------------------------------------------------

// Pass a truncated sequence to the p-vector algorithm to maintain consistancy,
xnm1_t = xnm1(1:numberOfSamples);

// run p-vector adaptive filter.
W = lms_pvector(xnm1_t,Rdx,mu1,p);

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Plot results.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
//scf(1);

