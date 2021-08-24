//**********************************************************************
// File Name: C4_5.sci
//**********************************************************************

// Bring in all of the utility functions.
exec('utils.sci',-1);

//**********************************************************************
//
//  Name: modifiedYuleWalker
//
//  Purpose: The purpose of this function is to solve for the
//  coefficients of a filter that models a signal given the number of
//  poles and zeros of the model.
//  The signal is modeled as the white noise response of a system
//  represented by, H(z) = B(z) / A(z), such that the coefficients of
//  B(z) and A(z) are contained in the vectors, [b(0), b(1),.. b(q)],
//  and [1 a(1), a(2),... a(p)] respectively.
//
//  Calling Sequence: [a,b] = pade(x,p,q)
//
//  Inputs:
//
//    x - The input vector to be processed.
//
//    p - The number of poles in the model.
//
//    q - The number of zeros in the model.
//
//  Outputs:
//
//    a - The denominator coefficients in the model.
//
//    b - The numerator coefficients in the model.
//
//**********************************************************************
function [a,b] = modifiedYuleWalker(x,p,q)

  // Compute autocorrelation sequence.
  rx = convol(x,x($:-1:1));

  // Find rx(0).
  n = find(rx == max(rx));

  // Use only rx(0) and the positive lags.
  rx = rx(n:$);

  // Let Pade' do the work with b being a dummy value.
  [a,b] = pade(rx,p,q);

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Determine the numerator coefficients.
  // Note that if the system function is
  // H(z) = B(z) / A(z), we can pass the
  // x(n) through a filter with the system
  // function described by H1(z) = A(z).
  // That is, define y(n) = h1(n) * x(n).  The
  // result will be an MA(q) process.  We can
  // then present this process to the Durbin
  // algorithm to estimate b(n).
  // Note that when invoking Durbin's method,
  // the all-pole model needs to be at least
  // four times the order of the estimated
  // all-zero model.  We choose 20q.
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  y = filterBlock(x,a,0);

  // Durbin's method does the rest of the work.  
  b = durbin(y,20*q,q);
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

endfunction

//**********************************************************************
// Mainline code.
//**********************************************************************

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (a): The writing of our function,
// modifiedYuleWalker().
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (b)
// Create x(n).  The spectrum will be plotted
// later.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Generate 100 samples of white Gaussian noise with unit variance.
noisegen(1,100,1);
w = feval([1:100],Noise);

// Set the numerator coefficients of H(z).
b = [1 -0.9 0.81];

// Set the denominator coefficients of H(z).
a = [-1.978 2.853 -1.877 0.904];

// Filter the noise.
x = filterBlock(w,b,a);

// Compute the autocorrelation function.
rx = convol(x,x($:-1:1));

// Compute the power spectrum.
Px = fft(rx,-1);
Px = Px(1:length(Px)/2 + 1);
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (c)
// Use the modified Yule-Walker equation to
// find an ARMA(4.2) model for x(n).
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
[ahat,bhat] = modifiedYuleWalker(x,4,2);

xhat = filterBlock(w,bhat,ahat(2:$));

//-----------------------------------------
// Now create 10 realizations of the model
// so that some statistical analysis can be
// carried out.
//-----------------------------------------
as = zeros(5,10);
bs = zeros(2,10);

for i = 1:10
 // Generate 100 samples of white Gaussian noise with unit variance.
  noisegen(1,100,1);
  ws = feval([1:100],Noise);

  // Filter the noise.
  xs = filterBlock(ws,b,a);

  // Generate the model coefficients.
  [as bs] = modifiedYuleWalker(xs,4,2);

  // Save in matrix.
  am(:,i) = as;
  bm(:,i) = bs;
end

// Perform statistical analysis.
for i = 1:5
  meanA(i) = mean(am(i,:));
  stA(i) = stdev(am(i,:));
end

for i = 1:3
  meanB(i) = mean(bm(i,:));
  stB(i) = stdev(bm(i,:));
end
//-----------------------------------------

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/


//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Plot the output so that comparasions can
// be made.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
subplot(311);
title('Original Spectrum, Px(w)');
plot(abs(Px));

subplot(312);
title('Original, x(n) (normalized)');
plot(x / max(x));

subplot(313);
title('Modified Yule-Walker Equations, xhat(n) (normalized)');
plot(xhat / max(xhat));

