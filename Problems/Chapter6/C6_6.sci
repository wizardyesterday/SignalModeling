//**********************************************************************
// File Name: C6_6.sci
//**********************************************************************

exec('utils.sci',-1);

//**********************************************************************
//
//  Name: generateSignals
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
//    len - The length of the sequences.
//
//  Outputs:
//
//    y - The sinusoid plus noise sequence.
//
//**********************************************************************
function A = generateSignals(len)

  // Construct time vector.
  n = 0:len-1;

  // Generate a Gaussian white noise sequence.
  noisegen(1,len,0.1);
  v = feval(1:len,Noise);

  // Generate noise-free signal.
  x = cos(2*%pi*0.35*n);

  // Add to the collection of signals.
  A = x';

  // Generate noisy signal.
  x = x + v;

  // Add to the collection of signals.
  A = [A x'];

  // Generate a decaying exponential.
  x = exp(-0.05*n);

  // Add to the collection of signals.
  A = [A x'];

endfunction

//**********************************************************************
// Mainline code.
//**********************************************************************
// Generate the collection of signals.
A = generateSignals(100);

// Generate Itkura models.
[g1Ita,e1Ita] = itakura(A(:,1),2);
[g2Ita,e2Ita] = itakura(A(:,2),2);
[g3Ita,e3Ita] = itakura(A(:,3),2);

// Generate forward covariance models.
[g1Fcov,e1Fcov] = fcov(A(:,1),2);
[g2Fcov,e2Fcov] = fcov(A(:,2),2);
[g3Fcov,e3Fcov] = fcov(A(:,3),2);

// Generate Burg models.
[g1Burg,e1Burg] = burg(A(:,1),2);
[g2Burg,e2Burg] = burg(A(:,2),2);
[g3Burg,e3Burg] = burg(A(:,3),2);


//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Print results.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
printf("\n--------------------------------------------------------\n");
printf("          x(n) = cos(2*pi*0.35*n)\n");
printf("--------------------------------------------------------\n");
printf("Ita e: %3.2e  Fcov e: %3.2e  Burg e: %3.2e\n", ...
        e1Ita($),e1Fcov($),e1Burg($));

printf("Ita g: %3.2e  Fcov g: %3.2e Burg g: %3.2e\n", ...
        g1Ita,g1Fcov,g1Burg);

printf("\n--------------------------------------------------------\n");
printf("          x(n) = cos(2*pi*0.35*n) + w(n)\n");
printf("--------------------------------------------------------\n");
printf("Ita e: %3.2e  Fcov e: %3.2e  Burg e: %3.2e\n", ...
        e2Ita($),e2Fcov($),e2Burg($));

printf("Ita g: %3.2e  Fcov g: %3.2e Burg g: %3.2e\n", ...
        g2Ita,g2Fcov,g2Burg);

printf("\n--------------------------------------------------------\n");
printf("          x(n) = exp(-0.05*n)\n");
printf("--------------------------------------------------------\n");
printf("Ita e: %3.2e  Fcov e: %3.2e  Burg e: %3.2e\n", ...
        e3Ita($),e3Fcov($),e3Burg($));

printf("Ita g: %3.2e  Fcov g: %3.2e Burg g: %3.2e\n", ...
        g3Ita,g3Fcov,g3Burg);


