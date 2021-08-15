//**********************************************************************
// File Name: C4_1.sci
//**********************************************************************

// Bring in all of the utility functions.
exec('utils.sci',-1);

//**********************************************************************
//
//  Name: performPadeProc
//
//  Purpose: The purpose of this function is to create lowpass filters
//  using the Pade approximation.  Plots are created for various values
//  of p and q.  The heavy lifting is performed by the pade() function.
//
//  Calling Sequence: performPadeProc(omegaC,n0,scaleFactor,
//                                    displayNumber)
//
//  Inputs:
//
//    omegaC - The cuttoff frequency in rad/sample.
//
//    n0 - The delay in samples.
//
//    scaleFactor - The factor for which the sinc() function is
//    dividec.
//
//    displayNumber - The identifier of the display for which plots
//    are to be rendered.
//
//  Outputs:
//
//    None.
//
//**********************************************************************
function performPadeProc(omegaC,n0,scaleFactor,displayNumber)

  // Generate time vector.
  n = 0:19;

  // Delay the time by n0 samples.
  nShifted = n - n0;

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Generate ideal impulse responses.
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  h = sinc(nShifted * omegaC) / scaleFactor;

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Generate filters.
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // p = 0, q = 19.
  [a_0,b_19] = pade(h,0,19);

  // p = 4, q = 15.
  [a_4,b_15] = pade(h,4,15);

  // p = 8, q = 11.
  [a_8,b_11] = pade(h,8,11);

  // p = 12, q = 7.
  [a_12,b_7] = pade(h,12,7);

  // p = 16, q = 3.
  [a_16,b_3] = pade(h,16,3);

  // p = 19, q = 0.
  [a_19,b_0] = pade(h,19,0);
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Construct frequency responses.
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  [m0_19,fr] = frmag(b_19,a_0,200);
  [m4_15,fr] = frmag(b_15,a_4,200);
  [m8_11,fr] = frmag(b_11,a_8,200);
  [m12_7,fr] = frmag(b_7,a_12,200);
  [m16_3,fr] = frmag(b_3,a_16,200);
  [m19_0,fr] = frmag(b_0,a_19,200);
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

  // Scale frequency to f/%pi.
  fr = fr * 2;

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Plot the results.
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  scf(displayNumber);

  // Build up string of parameters.
  s2 = msprintf("wc: %1.2f, n0: %d, ",omegaC,n0);

  subplot(321);
  s1 = "Pade LP Filter, p: 0, q: 19, ";
  title(s1+s2);
  plot(fr,20*log10(m0_19));

  subplot(322);
  s1 = "Pade LP Filter, p: 4, q: 15, ";
  title(s1+s2);
  plot(fr,20*log10(m4_15));

  subplot(323);
  s1 = "Pade LP Filter, p: 8, q: 11, ";
  title(s1+s2);
  plot(fr,20*log10(m8_11));

  subplot(324);
  s1 = "Pade LP Filter, p: 12, q: 7, ";
  title(s1+s2);
  plot(fr,20*log10(m12_7));

  subplot(325);
  s1 = "Pade LP Filter, p: 16, q: 3, ";
  title(s1+s2);
  plot(fr,20*log10(m16_3));

  subplot(326);
  s1 = "Pade LP Filter, p: 19, q: 0, ";
  title(s1+s2);
  plot(fr,20*log10(m19_0));
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

endfunction

//**********************************************************************
//
//  Name: performPronyProc
//
//  Purpose: The purpose of this function is to create lowpass filters
//  using the Prony approximation.  Plots are created for various values
//  of p and q.  The heavy lifting is performed by the prony() function.
//
//  Calling Sequence: performPronyProc(omegaC,n0,scaleFactor,
//                                     displayNumber)
//
//  Inputs:
//
//    omegaC - The cuttoff frequency in rad/sample.
//
//    n0 - The delay in samples.
//
//    scaleFactor - The factor for which the sinc() function is
//    dividec.
//
//    displayNumber - The identifier of the display for which plots
//    are to be rendered.
//
//  Outputs:
//
//    None.
//
//**********************************************************************
function performPronyProc(omegaC,n0,scaleFactor,displayNumber)

  // Generate time vector.
  n = 0:19;

  // Delay the time by n0 samples.
  nShifted = n - n0;

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Generate ideal impulse responses.
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  h = sinc(nShifted * omegaC) / scaleFactor;

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Generate filters.
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // p = 0, q = 19.
  [a_0,b_19,eps0_19] = prony(h,0,19);

  // p = 4, q = 15.
  [a_4,b_15,eps4_15] = prony(h,4,15);

  // p = 8, q = 11.
  [a_8,b_11,eps8_11] = prony(h,8,11);

  // p = 12, q = 7.
  [a_12,b_7,eps12_7] = prony(h,12,7);

  // p = 16, q = 3.
  [a_16,b_3,eps16_3] = prony(h,16,3);

  // p = 19, q = 0.
  [a_19,b_0,eps19_0] = prony(h,19,0);
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Construct frequency responses.
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  [m0_19,fr] = frmag(b_19,a_0,200);
  [m4_15,fr] = frmag(b_15,a_4,200);
  [m8_11,fr] = frmag(b_11,a_8,200);
  [m12_7,fr] = frmag(b_7,a_12,200);
  [m16_3,fr] = frmag(b_3,a_16,200);
  [m19_0,fr] = frmag(b_0,a_19,200);
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

  // Scale frequency to f/%pi.
  fr = fr * 2;

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Plot the results.
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  scf(displayNumber);

  // Build up string of parameters.
  s2 = msprintf("wc: %1.2f, n0: %d, ",omegaC,n0);

  subplot(321);
  s1 = "Prony LP Filter, p: 0, q: 19, ";
  s3 = msprintf("e0_19: %0.4f",eps0_19);
  title(s1+s2+s3);
  plot(fr,20*log10(m0_19));

  subplot(322);
  s1 = "Prony LP Filter, p: 4, q: 15, ";
  s3 = msprintf("e4_15: %0.4f",eps4_15);
  title(s1+s2+s3);
  plot(fr,20*log10(m4_15));

  subplot(323);
  s1 = "Prony LP Filter, p: 8, q: 11, ";
  s3 = msprintf("e8_11: %0.4f",eps8_11);
  title(s1+s2+s3);
  plot(fr,20*log10(m8_11));

  subplot(324);
  s1 = "Prony LP Filter, p: 12, q: 7, ";
  s3 = msprintf("e12_7: %0.4f",eps12_7);
  title(s1+s2+s3);
  plot(fr,20*log10(m12_7));

  subplot(325);
  s1 = "Prony LP Filter, p: 16, q: 3, ";
  s3 = msprintf("e16_3: %0.4f",eps16_3);
  title(s1+s2+s3);
  plot(fr,20*log10(m16_3));

  subplot(326);
  s1 = "Prony LP Filter, p: 19, q: 0, ";
  s3 = msprintf("e19_0: %0.5f",eps19_0);
  title(s1+s2+s3);
  plot(fr,20*log10(m19_0));
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

endfunction

//**********************************************************************
//
//  Name: performEllipticProc
//
//  Purpose: The purpose of this function is to create a lowpass filter
//  using an elliptic approximation.  A plot is created to illustrate
//  the frequency response.
//
//  Calling Sequence: performEllipticProc(n,omegaC,ripple,displayNumber)
//
//  Inputs:
//
//    n - The filter order.
//
//    omegaC - The cuttoff frequency in rad/sample.
//
//    ripple - The passband and stopband ripple.
//
//    displayNumber - The identifier of the display for which plots
//    are to be rendered.
//
//  Outputs:
//
//    None.
//
//**********************************************************************
function performEllipticProc(n,omegaC,ripple,displayNumber);

  // Construct elliptic filter.  
  hs = analpf(n,'ellip',[ripple,ripple],omegaC);

  // Construct independent variable.
  fr = 0:1/200:1 - (1/200);
  fr = fr * %pi;

  // Generate frequency response vector.
  hf = freq(hs(2),hs(3),%i*fr);

  // Convert to magnitude response.
  hm = abs(hf);

  // Compute number of elements in passband.
  l = omegaC / %pi;
  l = l * length(hm);

  // Compute number of elements in stopband.
  u = length(hm) - l;

  // Create ideal magnitude response vector.
  ideal = [ones(1,l) zeros(1,u)];

  // Compute differences.
  difference = ideal - hm;

  // Square the differences.
  difference = difference .* difference;

  // Compute mean square value.
  eps = sum(difference);
  
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Plot the results.
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  scf(displayNumber);

  // Plot frequency response.
  subplot(211)
  s1 = "Elliptic LP Filter, ";
  s2 = msprintf("wc: %1.2f, ",omegaC);
  title(s1+s2);
  plot(fr,hm);

  // Plot squared error.
  subplot(212)
  s1 = "Elliptic LP Filter, ";
  s2 = msprintf("epsilon: %f, ",eps);
  title(s1+s2);
  plot(difference);

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

endfunction

//**********************************************************************
// Mainline code.
//**********************************************************************

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (a):
// Design a filter with a cutoff frequency
// of PI/2 using the Pade' method.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
performPadeProc(%pi/2,9,2,1);
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (b):
// Design a filter with a cutoff frequency
// of PI/16 using the Pade' method.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
performPadeProc(%pi/16,9,16,2);
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (c):
// Design a filter with a cutoff frequency
// of PI/2 using the Prony method.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

performPronyProc(%pi/2,9,2,3);
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (d):
// Design a filter using the builtin
// analpf() command.  This will be a 10th
// order elliptic filter with a cutoff
// frequency of PI/2.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
performEllipticProc(10,%pi/2,10^(-1.25),4);
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
