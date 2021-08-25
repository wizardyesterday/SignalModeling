//**********************************************************************
// File Name: C4_2.sci
//**********************************************************************

// Bring in all of the utility functions.
exec('utils.sci',-1);

//**********************************************************************
//
//  Name: performIpfProc
//
//  Purpose: The purpose of this function is to create models using
//  the iterative prefiltering method.
//
//  Calling Sequence: performIpfProc(h,sigma2,n,displayNumber)
//
//  Inputs:
//
//    h - The unit sample response to be modeled.
//
//    sigma2 - The variance of the noise to be added to h.
//
//    n - The number of iterations that the iterative prefilter
//    method should execute.
//
//    displayNumber - The identifier of the display for which plots
//    are to be rendered.
//
//  Outputs:
//
//    None.
//
//**********************************************************************
function performIpfProc(h,sigma2,n,displayNumber)

  N = length(h);

  if sigma2 == 0
    v = zeros(1,N);
  else
    noisegen(1,N,sqrt(sigma2));

    // Generate Gaussian noise sequence.
    v = feval([1:N],Noise);
  end

  // Force column vectors.
  h = h(:);
  v = v(:);

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Two-zero, four-pole model.
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Generate the model.
  [a4,b2,err42] = ipf(h+v,4,2,n);

  // Compute unit sample response.
  h42 = filterBlock(u,b2,a4(2:$));
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Four-zero, four-pole model.
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Generate the model.
  [a4,b4,err44] = ipf(h+v,4,4,n);

  // Compute unit sample response.
  h44 = filterBlock(u,b4,a4(2:$));
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Five-pole, five-zero model.
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Generate the model.
  [a5,b5,err55] = ipf(h+v,5,5,n);

  // Compute unit sample response.
  h55 = filterBlock(u,b5,a5(2:$));
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Plot the results.
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  scf(displayNumber);

  subplot(421);
  s1 = "Original Filter, h(n), ";
  s2 = msprintf("sigma^2: %1.2e, err: 0.00",sigma2);  
  title(s1+s2);
  plot(h+v);

  subplot(423);
  s1 = "IPF Model, p: 4, q: 2, ";
  s2 = msprintf("sigma^2: %1.2e, err42: %1.2e",sigma2,err42);  
  title(s1+s2);
  plot(h42);

  subplot(425);
  s1 = "IPF Model, p: 4, q: 4, ";
  s2 = msprintf("sigma^2: %1.2e, err44: %1.2e",sigma2,err44);  
  title(s1+s2);
  plot(h44);

  subplot(427);
  s1 = "IPF Model, p: 5, q: 5, ";
  s2 = msprintf("sigma^2: %1.2e, err55: %1.2e",sigma2,err55);  
  title(s1+s2);
  plot(h55);
 //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

endfunction

//**********************************************************************
//
//  Name: performPronyProc
//
//  Purpose: The purpose of this function is to create models using
//  the Prony method.
//
//  Calling Sequence: performPronyProc(h,sigma2,displayNumber)
//
//  Inputs:
//
//    h - The unit sample response to be modeled.
//
//    sigma2 - The variance of the noise to be added to h.
//
///    displayNumber - The identifier of the display for which plots
//    are to be rendered.
///
//  Outputs:
//
//    None.
//
//**********************************************************************
function performPronyProc(h,sigma2,displayNumber)

  N = length(h);

  if sigma2 == 0
    v = zeros(1,N);
  else
    noisegen(1,N,sqrt(sigma2));

    // Generate Gaussian noise sequence.
    v = feval([1:N],Noise);
  end

  // Force column vectors.
  h = h(:);
  v = v(:);

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Two-zero, four-pole model.
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Generate the model.
  [a4,b2,err42] = prony(h+v,4,2);

  // Compute unit sample response.
  h42 = filterBlock(u,b2,a4(2:$));
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Four-zero, four-pole model.
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Generate the model.
  [a4,b4,err44] = prony(h+v,4,4);

  // Compute unit sample response.
  h44 = filterBlock(u,b4,a4(2:$));
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Five-pole, five-zero model.
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Generate the model.
  [a5,b5,err55] = prony(h+v,5,5);

  // Compute unit sample response.
  h55 = filterBlock(u,b5,a5(2:$));
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Plot the results.
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  scf(displayNumber);

  subplot(422);
  s1 = "Original Filter, h(n), ";
  s2 = msprintf("sigma^2: %1.2e, err: 0.00",sigma2);  
  title(s1+s2);
  plot(h+v);

  subplot(424);
  s1 = "Prony Model, p: 4, q: 2, ";
  s2 = msprintf("sigma^2: %1.2e, err42: %1.2e",sigma2,err42);
  title(s1+s2);
  plot(h42);

  subplot(426)
  s1 = "Prony Model, p: 4, q: 4, ";
  s2 = msprintf("sigma^2: %1.2e, err44: %1.2e",sigma2,err44);
  title(s1+s2);
  plot(h44);

  subplot(428);
  s1 = "Prony Model, p: 5, q: 5, ";
  s2 = msprintf("sigma^2: %1.2e, err55: %1.2e",sigma2,err55);  
  title(s1+s2);
  plot(h55);
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

endfunction

//**********************************************************************
// Mainline code.
//**********************************************************************

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (a):
// Generate the first 100 samples of the
// unit sample response of a filter.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Generate a unit sample sequence.
u = [1 zeros(1,99)];

// Construct numerator.
b = [1 -0.9 0.81];

// Construct denominator.
a = [-1.978 2.853 -1.877 0.9036];

// Compute unit sample response.
h = filterBlock(u,b,a);
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (b):
// Use iterative prefiltering to generate
// some models.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
performIpfProc(h,0,6,1);
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (c):
// Use iterative prefiltering to generate
// some models with noise added to the signal.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
performIpfProc(h,0.0001,6,2);
performIpfProc(h,0.001,6,3);
performIpfProc(h,0.01,6,4);
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (d):
// Use Prony to generate some models with and
// without noise.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
performPronyProc(h,0,1);
performPronyProc(h,0.0001,2);
performPronyProc(h,0.001,3);
performPronyProc(h,0.01,4);
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
