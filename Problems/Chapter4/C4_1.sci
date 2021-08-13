//**********************************************************************
// File Name: C4_1.sci
//**********************************************************************

// Bring in all of the utility functions.
exec('utils.sci',-1);

//**********************************************************************
//
//  Name: performPadeProcessing
//
//  Purpose: The purpose of this function is to create lowpass filters
//  using the Pade approximation.  Plots are created for various values
//  of p and q.
//
//  Calling Sequence: performPadeProcessing()
//
//  Inputs:
//
//    None.
//
//  Outputs:
//
//    None.
//
//**********************************************************************
function performPadeProcessing(omegaC,n0,displayNumber)

  // Generate arguments to the sinc() function.
  n = 0:19;

  // Generate n - n0 (n0 = 9).
  nShifted = n - n0;

  // Generate ideal impulse responses.
  i9 = sinc(nShifted * omegaC) / 2;

  // Generate truncated impulse response.
  h9 = [i9(1:10) zeros(1,10)];

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Generate filters.
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // p = 0, q = 19.
  [a_0,b_19] = pade(h9,0,19);

  // p = 4, q = 15.
  [a_4,b_15] = pade(h9,4,15);

  // p = 8, q = 11.
  [a_8,b_11] = pade(h9,8,11);

  // p = 12, q = 7.
  [a_12,b_7] = pade(h9,12,7);

  // p = 16, q = 3.
  [a_16,b_3] = pade(h9,16,3);

  // p = 19, q = 0.
  [a_19,b_0] = pade(h9,19,0);

  // Construct frequency responses.
  [m0_19,fr] = frmag(b_19,a_0,200);
  [m4_15,fr] = frmag(b_15,a_4,200);
  [m8_11,fr] = frmag(b_11,a_8,200);
  [m12_7,fr] = frmag(b_7,a_12,200);
  [m16_3,fr] = frmag(b_3,a_16,200);
  [m19_0,fr] = frmag(b_0,a_19,200);

  // Scale frequency to f/%pi.
  fr = fr * 2;

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Plot the results.
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  scf(displayNumber);

  subplot(321);
  title('Pade-designed lowpass filter (p = 0, q = 19');
  plot(fr,20*log10(m0_19));

  subplot(322);
  title('Pade-designed lowpass filter (p = 4, q = 15');
  plot(fr,20*log10(m4_15));

  subplot(323);
  title('Pade-designed lowpass filter (p = 8, q = 11');
  plot(fr,20*log10(m8_11));

  subplot(324);
  title('Pade-designed lowpass filter (p = 12, q = 7');
  plot(fr,20*log10(m12_7));

  subplot(325);
  title('Pade-designed lowpass filter (p = 16, q = 3');
  plot(fr,20*log10(m16_3));

  subplot(326);
  title('Pade-designed lowpass filter (p = 19, q = 0');
  plot(fr,20*log10(m19_0));

endfunction

//**********************************************************************
// Mainline code.
//**********************************************************************

// Parform part (a) processing.
performPadeProcessing(%pi/2,1,1);

// Perform part (b) processing.
performPadeProcessing(%pi/16,1,2);

// Perform part (c) processing.

// Perform part (d) processing.





