//**********************************************************************
// File Name: C5_3.sci
//**********************************************************************

// Bring in all of the utility functions.
exec('utils.sci',-1);

//**********************************************************************
//
//  Name: inverseShur
//
//  Purpose: The purpose of this function is to compute the
//  autocorrelation sequency given the reflection coefficients and the
//  modeling error.  A backwards recursion is performed to recover the
//  autocorrelation sequence.
//  Presently, reflection coefficients lengths of 3 are supported.
//  This is still adequate for benchmark purposes.
//  This function will be enhanced in the fugure to handle any length
//  of reflection coefficient sequence.
//
//  Calling Sequence: r = inverseShur(gamm,epsilon)
//
//  Inputs:
//
//    a - The model coefficients.
//
//    e - The modeing error.
//
//  Outputs:
//
//    r - The output autocorrelation sequence.
//
//**********************************************************************
function r = inverseShur(gamm,epsilon)

  // Upper loop index.
  p = length(gamm);
  p1 = p + 1;

  if p <> 3
    error("Only sequence lengths of 3 are handled.");
  end

  r = zeros(1,p1)';

 // Compute |gamm|^2.
  gammSquared = gamm .* conj(gamm);

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // compute first autocorrelation value.
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  r(1) = e / prod(1 - gamm .* gamm);

  // Compute remaining autocorrelation values.
  r(2) = -gamm(1)*r(1);

  r(3) = -gamm(2)*r(1) - (gamm(1) + gamm(2)*gamm(1))*r(2);

  r(4) = -g(3)*r(1) ...
         - (g(2) + g(3)*g(1) + g(3)*g(2)*g(1))*r(2) ...
         - (g(1) + g(2)*g(1) + g(3)*g(2))*r(3);
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
 
 
endfunction


//**********************************************************************
// Mainline code.
//**********************************************************************

r = [1 .9 .8 .7]';

[g,e] = shur(r);

r1 = inverseShur(g,e);
disp([r r1]);




