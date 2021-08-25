//**********************************************************************
// File Name: C4_3.sci
//**********************************************************************

// Bring in all of the utility functions.
exec('utils.sci',-1);

//**********************************************************************
//
//  Name: generateWhiteNoise
//
//  Purpose: The purpose of this function is to generate a white noise
//  sequence for which the random variable is uniformly distributed
//  between -0.005 and 0.005.
//
//  Calling Sequence: w = generateWhiteNoise(numberOfSamples)
//
//  Inputs:
//
//    numberOfSampels - The number of random numbers to generate.
//
//  Outputs:
//
//    w - A uniformly distributed white noise sequence.
//
//**********************************************************************
function w = generateWhiteNoise(numberOfSamples)

  // Ensure that we have a uniform distribution.
  rand('uniform');

  for i = 1:numberOfSamples
    // Generate the next random number.
    a = rand();

    // Translate for symmetry about 0.
    a = a - .5;

    // Scale so that -0.005 <= a <= .005.
    a = a * .01;

    // Save in returned vector.
    w(i) = a;
  end

endfunction

//**********************************************************************
//
//  Name: createImpulseTrain
//
//  Purpose: The purpose of this function is to a stream of weighted
//  impulses for testing purposes.
//
//  Calling Sequence: x = createImpulseTrain()
//
//  Inputs:
//
//    None.
//
//  Outputs:
//
//    x - The impulse train..
//
//**********************************************************************
function x = createImpulseTrain()

  // Generate locations of impulses.
  nk = [25 40 55 65 85 95 110 130 140 155];

  // Generate impulse weights.
  xk = [1 .8 .7 .5 .7 .2 .9 .5 .6 .3];

  // This index is used for lookup purposes.
  k = 1:10;

  // Create impulse train.
  x(nk(k)) = xk(k);

  // Force column vector.
  x = x(:);

endfunction

//**********************************************************************
//
//  Name: createdBlurredSignal
//
//  Purpose: The purpose of this function is to pass a signal through
//  a channel so as to introduce distortion of the signal.
//
//  Calling Sequence: y = createdBlurredSignal(x,g,sigma2)
//
//  Inputs:
//
//    x - The input signal.
//
//    g - The unit sample response of the channel to be modeled.
//
//    sigma2 - The variance of the noise to be added to the channel.
//
//
//  Outputs:
//
//    y - The distorted signal.
//
//**********************************************************************
function y = createdBlurredSignal(x,g,sigma2)

  // Force column vectors.
  x = x(:);
  g = g(:);

  // Create distorted signal.
  y = convol(x,g);


  N = length(y);

  if sigma2 == 0
    v = zeros(1,N);
  else
    noisegen(1,N,sqrt(sigma2));

    // Generate Gaussian noise sequence.
    v = feval([1:N],Noise);
  end

  // Force column vector.
  y = y(:);
  v = v(:);

  // Add in noise.
  y = y + v;

endfunction

//**********************************************************************
//
//  Name: recoverSignal
//
//  Purpose: The purpose of this function is to recover a signal by
//  passing a distorted signal through an inverse filter of that of the
//  channel.  Depending whether or not alpha is passed to this function,
//  one of two spike() routines is invoked.
//
//  Calling Sequence: xhar = recoverSignal(y,g,n0,alpha)
//
//  Inputs:
//
//    y - The input signal.
//
//    g - The unit sample response of the channel to be modeled.
//
//    n0 - The delay that is used to optimize the inverse filter.
//
//    alpha - The prewhitening parameter.
//
//  Outputs:
//
//    xhat - An estimate of the original signal.
//
//**********************************************************************
function xhat = recoverSignal(y,g,n0,sigma2,alpha)

  N = length(g);

  if sigma2 == 0
    v = zeros(1,N);
  else
    noisegen(1,N,sqrt(sigma2));

    // Generate Gaussian noise sequence.
    v = feval([1:N],Noise);
  end

  // Force column vectors.
  g = g(:);
  v = v(:);

  // Introduce noise into the channel.
  g0 = g + v;

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // argn(2) returns is the number of arguments passed to the
  // function.
  // If 5 arguments were not passed to the function, it is
  // implied that the last parameter was not passed.
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  if argn(2) == 5
   // Create spiking filter with prewhitening.
    [hn,err] = spikeWithPrewhitening(g0,n0,50,alpha); 
  else
    // Create spiking filter.
    [hn,err] = spike(g0,n0,50);
  end
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

  // Recover x(n).
  xhat = convol(hn,y);

endfunction

//**********************************************************************
//
//  Name: spikeWithPrewhitening
//
//  Purpose: The purpose of this function is to compute a least squares
//  inverse FIR filter that approximates a delayed unit sample.  This
//  function differs in the implementation of the original spike()
//  function since it solves the system of equations using the
//  autocorrelation matrix for the vector g.  Additionally, this
//  function accepts a whitening parameter.  The equations are of the
//  form,
//
//  (Rg + alpha * I)h = g0,
//
//  where, Rg is the autocorrelation matrix of g, alpha is a
//  prewhitening parameter, I is an identify matrix of the same size
//  as Rg, h is the vector of coefficients, and g0 is the vector of
//  delayed values of g.
//
//  Calling Sequence: [h,err] = spikeWithPrewhitening(g,n0,n,alpha)
//
//  Inputs:
//
//    g - The filter for which an inverse it to be computed.
//
//    n0 - The delay of the unit sample to be approximated.
//
//    n - The order of the inverse filter.
//
//    alpha - The whitening parameter.
//
//  Outputs:
//
//    h - The inverse filter coefficients.
//
//    err - The mean square error of the approximation to the delayed
//    unit sample.
//
//**********************************************************************
function [h,err] = spikeWithPrewhitening(g,n0,n,alpha)

  // Force column vector.
  g = g(:);

  if alpha < 0
    // Nuke the prewhitening parameter.
    alpha = 0;
  end

  if n > n0
    // Construct data matrix.
    G = convm(g,n);

    // Construct autocorrelation matrix.
    Rg = G' * G;

    //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
    // Extract the number of columns in Rg.  This is
    // necessary since other vectors, used to compute
    // the solution to the system of equations need
    // to be of sized consistant for the computations.
    //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
    // Determine the number of rows and columns in Rg.
    k = size(Rg);

    // Extract the number of rows of Rg.
    numberOfColumnsInRg = k(2);
    //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

    // Add in the whitening parameter.
    Rg = Rg + alpha * eye(numberOfColumnsInRg);

    //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
    // Construct g0(n).
    //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
    // Allocate space for desired vector, d(n).
    d = zeros(numberOfColumnsInRg,1);

    // Allocate space for expanded g(n).
    gx = zeros(numberOfColumnsInRg,1);

    // Copy g(n) to expanded version of g(n).
    gx(1:length(g)) = g(1:$);

    // Construct the g0(n) vector.
    g0(1:n0+1) = gx(n0+1:-1:1);
    //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

    // Compute coefficients of inverse filter.
    h = Rg \ g0;

    // Compute error.
    err = 1 - G(n0+1,:) * h;
  else
    error('Delay too large');
  end

endfunction

//**********************************************************************
// Mainline code.
//**********************************************************************

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (a):
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Generate time vector.
n = 0:50;

// Generate blurring filter.
g = cos(0.2 * (n - 25)) .* exp(-0.01 * (n - 25)^2);

// Create original signal.
x = createImpulseTrain();

// Create noise-free distorted signal.
y = createdBlurredSignal(x,g,0);
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (b):
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Create noise-free distorted signal.
y = createdBlurredSignal(x,g,0);
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (c):
// Recover the original signal with an
// optimum delay of 49 samples.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Recover the signal assumming a noise-free g(n).
xhat = recoverSignal(y,g,49,0);
xhat = xhat(50:50+154);
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (d):
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Create noisy distorted signals.
y_0001 = createdBlurredSignal(x,g,0.0001);
y_001 = createdBlurredSignal(x,g,0.001);

// Recover original signal in both cases.
xhat_0001 = recoverSignal(y_0001,g,49,0);
xhat_001 = recoverSignal(y_001,g,49,0);

// Keep only the portions of interest.
xhat_0001 = xhat_0001(50:50+154);
xhat_001 = xhat_001(50:50+154);
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (e):
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
white_0001_50 = recoverSignal(y_0001,g,49,0,50);
white_0001_500 = recoverSignal(y_0001,g,49,0,500);
white_001_50 = recoverSignal(y_001,g,49,0,50);
white_001_500 = recoverSignal(y_001,g,49,0,500);

white_0001_50 = white_0001_50(50:50+154);
white_0001_500 = white_0001_500(50:50+154);

white_001_50 = white_001_50(50:50+154);
white_001_500 = white_001_500(50:50+154);

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (f):
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Generate uniformly distributed white noise.
w = generateWhiteNoise(length(g));

// Create noisy blurring filter.
gNoisy = g + w';

xhatNoisyG = recoverSignal(y,gNoisy,49,0);
xhatNoisyG = xhatNoisyG(50:50+154);
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Make some plots.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
subplot(331);
title('Original x(n)');
plot2d3([1:length(x)],x,style=[color("red")]);

subplot(332);
title('Channel-impaired x(n)');
plot(y);

subplot(333);
title('Recovered, sigma^2: 0.0001');
plot(xhat_0001 / max(xhat_0001));

subplot(334);
title('Recovered, sigma^2: 0.001');
plot(xhat_001 / max(xhat_001));

subplot(335);
title('Recovered, sigma^2: 0.0001, alpha: 50');
plot(white_0001_50 / max(white_0001_50));

subplot(336);
title(' Recovered, sigma^2: 0.0001, alpha: 500');
plot(white_0001_500 / max(white_0001_500));

subplot(337);
title('Recovered, sigma^2: 0.001, alpha: 50');
plot(white_001_50 / max(white_001_50));

subplot(338);
title('Recovered, sigma^2: 0.001, alpha: 500');
plot(white_001_500 / max(white_001_500));

subplot(339);
title('Recovered, Noisy g(n)');
plot(xhatNoisyG / max(xhatNoisyG));



