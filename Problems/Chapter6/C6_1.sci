//**********************************************************************
// File Name: C6_1.sci
//**********************************************************************

exec('utils.sci',-1);

//**********************************************************************
//
//  Name: generateCosineWithNoise
//
//  Purpose: The purpose of this function is to generate a cosine
//  signal with additive Gaussian white noise.
//
//  Calling Sequence: y = generateCosineWithNoise(a,f0,phi,sigma,len)
//
//  Inputs:
//
//    a - The amplitude of the sinusoid.
//
//    f0 - The frequency of the sinusoid in cycles per sample.
//
//    phi - The phase of the sinusoide in radians.
//
//    sigma - The standard deviation of the Gaussian noise sequence
//    that is to be added to the sinusoid.
//
//    len - The length of the sequence.
//
//  Outputs:
//
//    y - The sinusoid plus noise sequence.
//
//**********************************************************************
function y = generateCosineWithNoise(a,f0,phi,sigma,len)

  // Construct time vector.
  n = 0:len-1;

  // Construct phase vector.
  theta = phi * ones(1,len);

  y = a * cos(2*%pi*f0*n + theta);

  if sigma > 0
    // Generate the Gaussian white noise sequence.
    noisegen(1,len,sigma);
    v = feval(1:len,Noise);

    // Add the noise to the signal.
    y = y + v;
  end

endfunction

//**********************************************************************
//
//  Name: estimateFrequency
//
//  Purpose: The purpose of this function is to estimate the frequency
//  of a sinusoid.
//
//  Calling Sequence: f0 = estimateFrequency(y,algorithmPtr)
//
//  Inputs:
//
//    y - A sinusoid combined with noise.
//
//    algorithmPtr - A pointer to the algorithm function that is to be
//    used for frequency estimation.
//
//  Outputs:
//
//    f0 - The frequency of the sinusoid in cycles per sample.
//
//**********************************************************************
function f0 = estimateFrequency(y,algorithmPtr)

  // Invoke the appropriate algorithm.
  f0 = algorithmPtr(y);

endfunction

//**********************************************************************
//
//  Name: constrainedLattic
//
//  Purpose: The purpose of this function is to estimate the frequency
//  of a sinusoid using the constrained lattice structure.  The
//  sinusoid may or may not be combined with noise. What's interesting
//  about this function is that it uses a second-order FIR lattice with
//  Gamma2 set equal to 1. This places the zeros of A(z) on the unit
//  circle.  To estimate the frequency of the sinusoid, all that is
//  needed is to construct a2 = [1 2*Gamma2 1]', solve for the roots,
//  and compute w0 = PI - Arg(z1)  - the first of of the two roots.
//  The frequency f0 is set equal to w0 / (2pi).
//
//  Calling Sequence: f0 = constrainedLattice(y)
//
//  Inputs:
//
//    y - A sinusoid combined with noise.
//
//  Outputs:
//
//    f0 - The frequency estimate of the sinusoid in cycles per
//    sample.
//
//**********************************************************************
function f0 = constrainedLattice(y)

  N = length(y);

  // Initialize sums.
  num = 0;
  den = 0;

  // Compute numerator.
  for n = 3:N
    num = num + (y(n) + y(n-2)) * conj(y(n-1));
  end

  // Construct denominator.
  for n = 3:n
    den = den + y(n-1) * conj(y(n-1));
  end

  // Scale denominator.
  den = den * 2;

  // Compute reflection coefficient.
  g1 = num ./ den;

  // Construct polynomial, and find the roots.
  a2 = [1 2*g1 1];
  r = roots(a2);

  // Compute estimated frequency in radians/sample.
  w0 = %pi - atan(imag(r(1)),real(r(1)));

  // Convert to cycles/sample.
  f0 = w0 / (2 * %pi);

endfunction

//**********************************************************************
//
//  Name: burgEstimator
//
//  Purpose: The purpose of this function is to estimate the frequency
//  of a sinusoid using the Burg algorithm.  The sinusoid may or may
//  not be combined with noise.  Once the reflection coefficients are
//  computed, they are converted to the coefficients of A2(z).  To
//  estimate the frequency of the sinusoid, all that is needed is to
//  solve for the roots of A2(z), and compute w0 = Arg(z1), where z1 is
//  the first of the two roots of A2(z). The frequency f0 is set equal
//  to w0 / (2pi).
//
//  Calling Sequence: f0 = burgEstimator(y)
//
//  Inputs:
//
//    y - A sinusoid combined with noise.
//
//  Outputs:
//
//    f0 - The frequency estimate of the sinusoid in cycles per
//    sample.
//
//**********************************************************************
function f0 = burgEstimator(y)

  // Compute reflection coefficients.
  [gamm,err] = burg(y,2);

  // Convert reflection coefficients to model parameters.
  a2 = gtoa(gamm);

  // Find the roots.
  r = roots(a2);

  // Compute estimated frequency in radians/sample.
  w0 = atan(imag(r(1)),real(r(1)));

  // Convert to cycles/sample.
  f0 = w0 / (2 * %pi);

endfunction

//**********************************************************************
//
//  Name: mcovEstimator
//
//  Purpose: The purpose of this function is to estimate the frequency
//  of a sinusoid using the modified covariance algorithm.  The sinusoid
//  may or may not be combined with noise.  Once the reflection
//  coefficients are computed, they are converted to the coefficients
//  of A2(z).  To estimate the frequency of the sinusoid, all that is
//  needed is to solve for the roots of A2(z), and compute w0 = Arg(z1),
//  where z1 is the first of the two roots of A2(z). The frequency f0
//  is set equal to w0 / (2pi).
//
//  Calling Sequence: f0 = mcovEstimator(y)
//
//  Inputs:
//
//    y - A sinusoid combined with noise.
//
//  Outputs:
//
//    f0 - The frequency estimate of the sinusoid in cycles per
//    sample.
//
//**********************************************************************
function f0 = mcovEstimator(y)

  // Compute reflection coefficients.
  [gamm,err] = mcov(y,2);

  // Convert reflection coefficients to model parameters.
  a2 = gtoa(gamm);

  // Find the roots.
  r = roots(a2);

  // Compute estimated frequency in radians/sample.
  w0 = atan(imag(r(1)),real(r(1)));

  // Convert to cycles/sample.
  f0 = w0 / (2 * %pi);

endfunction

//**********************************************************************
//
//  Name: acmEstimator
//
//  Purpose: The purpose of this function is to estimate the frequency
//  of a sinusoid using the modified autocorrelation method.  The sinusoid
//  may or may not be combined with noise.  Once the reflection
//  coefficients are computed, they are converted to the coefficients
//  of A2(z).  To estimate the frequency of the sinusoid, all that is
//  needed is to solve for the roots of A2(z), and compute w0 = Arg(z1),
//  where z1 is the first of the two roots of A2(z). The frequency f0
//  is set equal to w0 / (2pi).
//
//  Calling Sequence: f0 = mcovEstimator(y)
//
//  Inputs:
//
//    y - A sinusoid combined with noise.
//
//  Outputs:
//
//    f0 - The frequency estimate of the sinusoid in cycles per
//    sample.
//
//**********************************************************************
function f0 = acmEstimator(y)

  // Compute the model parameters.
  [a2,err] = acm(y,2)

  // Find the roots.
  r = roots(a2);

  // Compute estimated frequency in radians/sample.
  w0 = atan(imag(r(1)),real(r(1)));

  // Convert to cycles/sample.
  f0 = w0 / (2 * %pi);

endfunction

//**********************************************************************
// Mainline code.
//**********************************************************************

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Evaluate the performance and
// compare to other algorithms.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Generate frequency vector.
f0 = 0.05:0.05:0.45;
f0 = f0';

// Generate noise vector.
v = [0 0.25 0.5 0.9]';

// Generate buffer length vector.
len = [100 500 1000 5000]';

// Perform frequency sweep with various noise values.
for k = 1:length(v)
  for i = 1:length(f0)
    y = generateCosineWithNoise(1,f0(i),0,v(k),100);
    fhatnoisy(k,i) = estimateFrequency(y,constrainedLattice);
  end
end

// Perform a frequency sweep with various buffer sizes and fixed noise value.
for k = 1:length(len)
  for i = 1:length(f0)
    y = generateCosineWithNoise(1,f0(i),0,0,len(k));
    fhatlen(k,i) = estimateFrequency(y,constrainedLattice);
  end
end

// Perform a frequency sweep on the other algorithms.
for i = 1:length(f0)
  y = generateCosineWithNoise(1,f0(i),0,0,100);
  fhatburg(i) = estimateFrequency(y,burgEstimator);
  fhatmcov(i) = estimateFrequency(y,mcovEstimator);
  fhatacm(i) = estimateFrequency(y,acmEstimator);
end

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Plot constrained lattice results.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Display the first graphic.
scf(1);

subplot(421);
title('Constrained, Noise: 0, Len: 100');
plot(f0,fhatnoisy(1,:));
xgrid();

subplot(422);
s1 = sprintf("Constrained, Noise: %f, ",v(1));
s2 = s1 + 'Len: 100';
title(s2);
plot(f0,fhatnoisy(2,:));
xgrid();

subplot(423);
s1 = sprintf("Constrained, Noise: %f, ",v(2));
s2 = s1 + 'Len: 100';
title(s2);
plot(f0,fhatnoisy(3,:));
xgrid();

subplot(424);
s1 = sprintf("Constrained, Noise: %f, ",v(4));
s2 = s1 + 'Len: 100';
title(s2);
plot(f0,fhatnoisy(4,:));
xgrid();

subplot(425);
s1 = sprintf("Constrained, Len: %d, ",len(1));
s2 = s1 + 'Noise: 0.5';
title(s2);
plot(f0,fhatlen(1,:));
xgrid();

subplot(426);
s1 = sprintf("Constrained, Len: %d, ",len(2));
s2 = s1 + 'Noise: 0.5';
title(s2);
plot(f0,fhatlen(2,:));
xgrid();

subplot(427);
s1 = sprintf("Constrained, Len: %d, ",len(3));
s2 = s1 + 'Noise: 0.5';
title(s2);
plot(f0,fhatlen(3,:));
xgrid();

subplot(428);
s1 = sprintf("Constrained, Len: %d, ",len(4));
s2 = s1 + 'Noise: 0.5';
title(s2);
plot(f0,fhatlen(4,:));
xgrid();

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Plot Burg, modified covariance, and
// autocorrelation results.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Display the second graphic.
scf(2);

subplot(311);
title('Burg, Len: 100, Noise: 0');
plot(f0,fhatburg);
xgrid();

subplot(312);
title('Modified Covariance, Len: 100, Noise: 0');
plot(f0,fhatmcov);
xgrid();

subplot(313);
title('Autocorrelation Method, Len: 100, Noise: 0');
plot(f0,fhatacm);
xgrid();
