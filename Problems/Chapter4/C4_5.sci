//**********************************************************************
// File Name: C4_5.sci
//**********************************************************************

// Bring in all of the utility functions.
exec('utils.sci',-1);

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

