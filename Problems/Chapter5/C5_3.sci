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
//  Presently something needs to be figured out how to determine
//  gR.j-1(p).  This function still has some use for benchmark purposes.
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

  // Create initial vectors.
  g = zeros(1:p+1)';
  gR = zeros(1,p+1)';
  gR(p+1) = epsilon;

 // Compute |gamm|^2.
  gammSquared = gamm .* conj(gamm);

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Compute g.j-1(k) and gR.j-1(k).
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  for j = p-1:-1:0

    //-------------------------------
    // Double buffering is needed so
    // that the past incarnations
    // of g(k) and gR(k) can be
    // utilized.  This is not so
    // obvious when looking at the
    // algorithm.
    //-------------------------------
    gPrev = g;
    gRPrev = gR;
    //-------------------------------

    // Compute g.j+1(k).
    g = (gPrev - gamm(j+1)*gRPrev) / (1 - gammSquared(j+1));

    //----------------------------------------------------
    // Compute all values of gR(k), k = 1, 2,...,p.  The
    // (p+1)st value cannot be computed since future
    // values of gRPrev(p+2) and gPrev(p+2) would be
    // needed.  These values do not exist in our finite
    // sequences.  Perhaps, this can be reconciled in the
    // future.
    //----------------------------------------------------
    for k = p:-1:j+1
      gR(k) = ...
        (gRPrev(k+1) - conj(gamm(j+1))*gPrev(k+1)) / (1 - gammSquared(j+1));
    end  
    //----------------------------------------------------

    // This cannot be computed since we don't know future values of the
    // right hand side.
    //gR = (gRPrev - conj(gamm(j+1))*gPrev) / (1 - gammSquared(j+1));

    //----------------------------------------------------
  end
 //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

  // Initialize the running product.
  d = 1;

  // Compute the denominator for computation of r(0).
  for j = 1:p
    d = d * (1 - gammSquared(j));
  end

  // Fill in what could not be computed in previous processing. 
  g(1) = epsilon / d;

  // Set returned value as the autocorrelation sequence.
  r = g;
 
endfunction


//**********************************************************************
// Mainline code.
//**********************************************************************

r = [1 .9 .8 .7]';

[g,e] = shur(r);

r1 = inverseShur(g,e);
disp([r r1]);




