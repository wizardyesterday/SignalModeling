//**********************************************************************
// File Name: C4_4.sci
//**********************************************************************

// Bring in all of the utility functions.
exec('utils.sci',-1);

//**********************************************************************
//
//  Name: constructMaModel
//
//  Purpose: The purpose of this function is to construct a moving
//  average model of order 2 using the method of spectral factorization.
//
//  Calling Sequence: b = constructMaModel(x)
//
//  Inputs:
//
//   x - The input sequence that is to be modeled.
//
//  Outputs:
//
//    b - The model parameters.
//
//**********************************************************************
function b = constructMaModel(x)

  // Compute autocorrelation sequence.
  rx = convol(x,x($:-1:1));

  // Find rx(0).
  n = find(rx == max(rx));

  // Only a second-order model is of interest.
  r = [rx(n-2) rx(n-1) rx(n) rx(n+1) rx(n+2)];

  // Perform our "z-transform".
  p = poly(r,'z','coeff');

  // Factor our polynomial.
  c = factors(p);

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // This is a kludge, but it's the easiest way to
  // deal with the fact that sometimes the
  // polynomial is factored into first-order
  // equations, and sometimes it is factored into
  // second-order equations.
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  if length(c) == 4
    // Expand to a second-order polynomial.
    m1 = c(1) * c(2);
  else
    // This is already a second-order polynomial.
    m1 = c(1);
  end
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

  // Retrieve the coefficients.
  b = coeff(m1);

  // Place in increasing powers of z^(-1).
  b = b($:-1:1);

  // Force column vectors.
  b = b(:);

endfunction

//**********************************************************************
//
//  Name: flexibleDurbin
//
//  Purpose: The purpose of this function is to solve for the
//  moving average model parameters for an input sequence using
//  Durbin's method.
//  The signal is modeled as the unit sample response of a system
//  represented by, H(z) = B(z) , such that the coefficients of
//  B(z)  are contained in the vector, [b(0), b(1),.. b(q)].
//  This function differs from the original dubin function due to the
//  fact that, rather than use the acm (autocorrelation) method in
//  the two steps, the caller can supply a function for each step,
//  individually.  This provides some flexiblity so that different
//  correlation functions can be used for experimentation.
//
//  Calling Sequence: b = flexibleDurbin(f1,f2,x,p,q)
//
//  Inputs:
//
//    f1 - A function that performs the first autocorrelation.
//
//    f2 - A function that performs the second correlation.
//
//    x - The input vector to be processed.
//
//    p - The order of the all-pole model used to approximate 1/B(z).
//    This parameter should be at least 4*q.
//
//    q - The order of the moving average model.
//
//  Outputs:
//
//    b - The parameters of the moving average model.
//
//**********************************************************************
function b = flexibleDurbin(f1,f2,x,p,q)

  // Ensure that we have a column vector.
  x = x(:);

  if p < length(x)
    //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
    // Use the supplied correlation method to estimate the
    // denominator of an all-pole model for x(n).
    //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
    [a,epsilon] = f1(x,p);

    //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
    // Use the supplied correlation method to estimate the model
    // parameters of an all-pole model for a(n).  This will
    // be an estimate of b(n).
    //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
    [b,epsilon] = f2(length(x)*a/sqrt(epsilon),q);
    b = b*length(x)/sqrt(epsilon);
  else
    error('Model order too large');
  end

endfunction

//**********************************************************************
// Mainline code.
//**********************************************************************

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (a): The writing of our function,
// constructMaModel().
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (b)
// Test our function.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Generate 256 samples of white Gaussian noise with unit variance.
noisegen(1,256,1);
w = feval([1:256],Noise);

// Filter the noise.
x = filterBlock(w,[1 0 0.9],0);

// Construct MA(2) model.
b = constructMaModel(x);
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (c)
// Use Durbin's method.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
bd = durbin(x,8,2);
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (d)
// Use Durbin's method, but with the
// second autocorrelation replaced with the
// covariance.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
bdac = flexibleDurbin(acm,covm,x,8,2);
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (d)
// Use Durbin's method, but with both
// autocorrelations replaced with the
// covariance.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
bdcc = flexibleDurbin(covm,covm,x,8,2);
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Compute  the white noise response of the
// filters.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
y = filterBlock(w,b,0);
yd = filterBlock(w,bd,0);
ydac = filterBlock(w,bdac,0);
ydcc = filterBlock(w,bdcc,0);
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Plot the output so that comparasions can
// be made.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
subplot(511);
title('Original x(n)');
plot(x);

subplot(512);
title('Spectral Factorization');
plot(y);

subplot(513);
title('Durbin Method');
plot(yd);

subplot(514);
title('Durbin Method: Acm, Covm');
plot(ydac);

subplot(515);
title('Durbin Method: Covm, Covm');
plot(ydcc);



