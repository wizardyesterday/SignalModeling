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
//  of a sinusoid.  The sinusoid may or may not be combined with noise.
//  What's interesting about this function is that it uses a second-
//  order FIR lattice with Gamma2 set equal to 1. This places the zeros
//  of A(z) on the unit circle.  To estimate the frequency of the
//  sinusoid, all that is needed is to construct a2 = [1 2*Gamma2 1]',
//  solve for the roots, and compute the Arg() of PI - the first of
//  of the two roots.  The frequency f0 is set equal to w0 / (2pi).  The
//  fact that Gamma2 has a value of unit ensures that all roots of the
//  predictor polynomial, A2(z), will lie on the unit circle.
//
//  Calling Sequence: f0 = constrainedLattice(y)
//
//  Inputs:
//
//    y - A sinusoid combined with noise.
//
//  Outputs:
//
//    f0 - The frequency of the sinusoid in cycles per sample.
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
  a = [1 2*g1 1];
  r = roots(a);

  // Compute estimated frequency in radians/sample.
  w0 = %pi - atan(imag(r(1)),real(r(1)));

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
f0 = 0:0.05:0.45

// Generate noise vector.
v = [0 0.25 0.5 0.9];

// Generate buffer length vector.
len = [100 500 1000 5000];

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
    y = generateCosineWithNoise(1,f0(i),0,0.5,len(k));
    fhatlen(k,i) = estimateFrequency(y,constrainedLattice);
  end
end

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Plot results.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
subplot(421);
title('Constrained, Noise: 0, Len: 100');
plot(f0,fhatnoisy(1,:));

subplot(422);
s1 = sprintf("Constrained, Noise: %f, ",v(1));
s2 = s1 + 'Len: 100';
title(s2);
plot(f0,fhatnoisy(2,:));

subplot(423);
s1 = sprintf("Constrained, Noise: %f, ",v(2));
s2 = s1 + 'Len: 100';
title(s2);
plot(f0,fhatnoisy(3,:));

subplot(424);
s1 = sprintf("Constrained, Noise: %f, ",v(4));
s2 = s1 + 'Len: 100';
title(s2);
plot(f0,fhatnoisy(4,:));

subplot(425);
s1 = sprintf("Constrained, Len: %d, ",len(1));
s2 = s1 + 'Noise: 0.5';
title(s2);
plot(f0,fhatlen(1,:));

subplot(426);
s1 = sprintf("Constrained, Len: %d, ",len(2));
s2 = s1 + 'Noise: 0.5';
title(s2);
plot(f0,fhatlen(2,:));

subplot(427);
s1 = sprintf("Constrained, Len: %d, ",len(3));
s2 = s1 + 'Noise: 0.5';
title(s2);
plot(f0,fhatlen(3,:));

subplot(428);
s1 = sprintf("Constrained, Len: %d, ",len(4));
s2 = s1 + 'Noise: 0.5';
title(s2);
plot(f0,fhatlen(4,:));



