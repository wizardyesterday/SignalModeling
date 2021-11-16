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

  a = [0.7 0.17 0.1789474 0.1888889]';

  disp([g gR]);

  // Compute |gamm|^2.
  gammSquared = gamm .* conj(gamm);

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Construct reflection coefficient
  // sequence.  Note that the indices
  // differ from the algorithm in the
  // book since the book assumes
  // zero-based arrays, whereas
  // Scilab assumes one-based arrays.
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

    g = (gPrev - gamm(j+1)*gRPrev) / (1 - gammSquared(j+1));

//    for k = j+1:p
//      gR(k+1) = gRPrev(k) + conj(gamm(j+1)) * gPrev(k+1);
//    end

    for k = p:-1:j+1
      gR(k) = ...
        (gRPrev(k+1) - conj(gamma(j+1))*gPrev(k+1)) / (1 - gammSquared(j+1));
    end  


//    gR = (gRPrev - conj(gamma(j+1))*gPrev) / (1 - gammSquared(j+1));

    printf("\nj: %d",j);
    disp([g gR]);

    //----------------------------------------------------
  end
 //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

  d = 1;
  for j = 1:p
    d = d * (1 - gammSquared(j));
  end

  g(1) = epsilon / d;
//  disp(g);

  r = g;
 
endfunction



//**********************************************************************
// Mainline code.
//**********************************************************************

r = [1 .9 .8 .7]';

[g,e] = shur(r);

r1 = inverseShur(g,e);



