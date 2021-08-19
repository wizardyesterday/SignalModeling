//**********************************************************************
// File Name: C4_3.sci
//**********************************************************************

// Bring in all of the utility functions.
exec('utils.sci',-1);

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

  // Generate location of impulses.
  nk = [25 40 55 65 85 95 110 130 140 155];

  // Generate template vector
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
    noisegen(1,N,sigma2);

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
//  channel.
//
//  Calling Sequence: xhar = recoverSignal(y,g,n0,sigma2)
//
//  Inputs:
//
//    y - The input signal.
//
//    g - The unit sample response of the channel to be modeled.
//
//    n0 - The delay that is used to optimize the inverse filter.
//
//    sigma2 - The variance of the noise to be added to the channel.
//
//
//  Outputs:
//
//    xhat - An estimate of the original signal.
//
//**********************************************************************
function xhat = recoverSignal(y,g,n0,sigma2)

  N = length(g);

  if sigma2 == 0
    v = zeros(1,N);
  else
    noisegen(1,N,sigma2);

    // Generate Gaussian noise sequence.
    v = feval([1:N],Noise);
  end

  // Force column vectors.
  g = g(:);
  v = v(:);

  // Introduce noise into the channel.
  g0 = g + v;

  // Create spiking filter.
  [hn,err] = spike(g0,n0,50);

  // Recover x(n).
  xhat = convol(hn,y);

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

// Create noise-free signal.
x = createImpulseTrain();

// Create noise-free distorted signal.
y = createdBlurredSignal(x,g,0)

// Recover the signal assumming a noise-free g(n).
xhat = recoverSignal(y,g,49,0);

xhat49 = recoverSignal(y,g,49,0);
xhat49 = xhat49 / max(xhat49);

xhat50 = recoverSignal(y,g,50,0);
xhat50 = xhat50 / max(xhat50);

subplot(211);
plot(x);
subplot(212);
plot(xhat / max(xhat));

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (b):
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (c):
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (d):
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (e):
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (f):
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

