//**********************************************************************
// File Name: C6_6.sci
//**********************************************************************

exec('utils.sci',-1);

//**********************************************************************
//
//  Name: generateSignals
//
//  Purpose: The purpose of this function is to generate a cosine
//  signal with additive Gaussian white noise.  Three signals are
//  generated:
//
//  1. A cosine signal.
//  2. A cosine signal combined with white Gaussian noise.
//  3. A decaying exponential signal.
//  4. A decaying exponential signal multiplied by a ramp.
//
//  Calling Sequence: A = generateSignals(sigma,len)
//
//  Inputs:
//
//    sigma - The standard deviation of the Gaussian noise sequence
//    that is to be added to the sinusoid.
//
//    len - The length of the sequences.
//
//  Outputs:
//
//    A - A set of signals for which each column contains a
//    particular signal.
//
//**********************************************************************
function A = generateSignals(sigma,len)

  // Construct time vector.
  n = 0:len-1;

  // Generate a Gaussian white noise sequence.
  noisegen(1,len,sigma);
  v = feval(1:len,Noise);

  // Generate noise-free signal.
  x = cos(2*%pi*0.1*n);

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

  // Generate a decaying exponential multiplied by a ramp.
  x = (n + 1) .* x;

  // Add to the collection of signals.
  A = [A x'];

endfunction

//**********************************************************************
//
//  Name: burgImproved
//
//  Purpose: The purpose of this function is to compute an all-pole
//  model for the input sequence, x(n).  The difference between this
//  algorithm and the Burg algorithm lies in the face that the
//  denominator used for the computation of the jth reflection
//  coefficient is performed recursively so that CPU cycles can be
//  saved.
//
//  Calling Sequence: [gamn,err] = burgImproved(x,p)
//
//  Inputs:
//
//    x - A vector of signal values that are to be modeled.
//
//    p - The order of the model.
//
//  Outputs:
//
//    gamm - The vector of reflection coefficients.
//
//   err - The vector of modeling errors.
//
//**********************************************************************
function [gamm,err] = burgImproved(x,p)

  // Force column vector.
  x = x(:);

  N = length(x);

  // Initialize forward prediction error.
  eplus  = x(2:N);

  // Initialize backward prediction error.
  eminus = x(1:N-1);

  // D1 = |eplus_0|^2 + |eminus_0|^2.
  D = eplus'*eplus + eminus'*eminus;

  N = N - 1;

  for j = 1:p;
    // Compute reflection coefficient.
    gamm(j) = -2*eminus'*eplus / D;

    // eplus{j}(n) = eplus{j-1}(n) + gamms{j}*eminus{j-1}(n-1).
    temp1 = eplus  + gamm(j)*eminus;

    // eminus{j}(n) = eminus{j-1}(n-1) + gamma*{j}*eplus{j-1}(n).
    temp2 = eminus + conj(gamm(j))*eplus;

    // Update error.
    err(j) = temp1'*temp1 + temp2'*temp2;

    // Update forward prediction error.
    eplus = temp1(2:N);

    // Update backward prediction error.
    eminus = temp2(1:N-1);

    // D_j+1 = D_j(1 - |gamma_j|^2) - |eplus_j(j)|^2 - |eminus_j(N)|^2.
    D = D * (1 - gamm(j) * conj(gamm(j))) - eplus(j) * conj(eplus(j)) ...
        - eminus($) * conj(eminus($));

    N = N - 1;
  end

endfunction

//**********************************************************************
// Mainline code.
//**********************************************************************
// Generate the collection of signals.
A = generateSignals(0.1,60);

// Generate Burg models.
[g1Burg,e1Burg] = burg(A(:,1),2);
[g2Burg,e2Burg] = burg(A(:,2),2);
[g3Burg,e3Burg] = burg(A(:,3),2);
[g4Burg,e4Burg] = burg(A(:,4),2);

// Generate more efficient Burg models.
[g1Burgi,e1Burgi] = burgImproved(A(:,1),2);
[g2Burgi,e2Burgi] = burgImproved(A(:,2),2);
[g3Burgi,e3Burgi] = burgImproved(A(:,3),2);
[g4Burgi,e4Burgi] = burgImproved(A(:,4),2);

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Print results.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
printf("\n--------------------------------------------------------\n");
printf("          x(n) = cos(2*pi*0.35*n)\n");
printf("--------------------------------------------------------\n");
printf("Burg e: %3.2e  Burgi e: %3.2e\n", ...
        e1Burg($),e1Burgi($));

printf("Burg g: %3.2e  Burgi g: %3.2e\n", ...
        g1Burg,g1Burgi);

printf("\n--------------------------------------------------------\n");
printf("          x(n) = cos(2*pi*0.35*n) + w(n), sigma = 0.1\n");
printf("--------------------------------------------------------\n");
printf("Burg e: %3.2e  Burg Improved e: %3.2e\n", ...
        e2Burg($),e2Burgi($));

printf("Burg g: %3.2e  Burg Improved g: %3.2e\n", ...
        g2Burg,g2Burgi);


printf("\n--------------------------------------------------------\n");
printf("          x(n) = exp(-0.05*n)\n");
printf("--------------------------------------------------------\n");
printf("Burg e: %3.2e  Burg Improved e: %3.2e\n", ...
        e3Burg($),e3Burgi($));

printf("Burg g: %3.2e  Burg Improved g: %3.2e\n", ...
        g3Burg,g3Burgi);


printf("\n--------------------------------------------------------\n");
printf("          x(n) = (n+1)*exp(-0.05*n)\n");
printf("--------------------------------------------------------\n");
printf("Burg e: %3.2e  Burg Improved e: %3.2e\n", ...
        e4Burg($),e4Burgi($));

printf("Burg g: %3.2e  Burg Improved g: %3.2e\n", ...
        g4Burg,g4Burgi);
