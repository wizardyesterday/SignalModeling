//**********************************************************************
// File Name: C5_3.sci
//**********************************************************************

// Bring in all of the utility functions.
exec('utils.sci',-1);

//**********************************************************************
//
//  Name: newAtoR
//
//  Purpose: The purpose of this function is to provide an alternate
//  way of retrieving an autocorrelation sequence from the model
//  coefficients and the modeling error.  Matrix, Ap, is first
//  constructed, and once constructed, the autocorrelation sequence
//  is computed from Ap and the modeling error, e.
//
//  Calling Sequence: r = newAtoR(a,e)
//
//  Inputs:
//
//    a - The model coefficients.
//
//    e - The modeing error.
//
//  Outputs:
//
//    r - The input autocorrelation sequence.
//
//**********************************************************************
function [gamm,epsilon] = shurInv(r)

  // Ensure that we have a column vector.
  r = r(:);

  // Upper loop index.
  p = length(r) - 1;

  // g0(k) = r(k).
  g = r;

  // gR(k) = r(k).
  gR = r;

  disp(0);
  disp([g gR]);

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Construct reflection coefficient
  // sequence.  Note that the indices
  // differ from the algorithm in the
  // book since the book assumes
  // zero-based arrays, whereas
  // Scilab assumes one-based arrays.
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  for j = 0:p-1
    // Update gamma
    gamm(j+1) = -g(j+2) / gR(j+1);

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

    //----------------------------------------------------
    // Let's vectorize this stuff for speed.  I still
    // want to keep the for loops around as comments so
    // that the spirit of the algorithm is not lost.
    // Vectorizing without these comments would make the
    // code unmaintainable since the loops are implicit
    // in the vectorized statements.
    //----------------------------------------------------
//    for k = j+2:p
//      g(k+1) = gPrev(k+1) + gamm(j+1) * gRPrev(k);
//    end

    if (j == 0)
      g(2:$) = 0;
      g(1) = g(1) * (1 - gamm(j+1)^2);
    else
      g(1:$) = 0;
    end
    // Vectorized version of above loop.
    g(j+3:p+1) = gPrev(j+3:p+1) + gamm(j+1) .* gRPrev(j+2:p);

//    for k = j+1:p
//      gR(k+1) = gRPrev(k) + conj(gamm(j+1)) * gPrev(k+1);
//    end

gR(1:$) = 0;
    // Vectorized version of above loop.
    gR(j+2:p+1) = gRPrev(j+1:p) + conj(gamm(j+1)) .* gPrev(j+2:p+1);

v = (g - gamm(j+1)*gR) / (1 - gamm(j+1)^2);

    disp(j+1);
    disp([v g gR]);

    //----------------------------------------------------
  end
 //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

  // eP = gP^(R)(p).
  epsilon= gR(p+1);
 
endfunction



//**********************************************************************
// Mainline code.
//**********************************************************************

r = [1 .9 .8 .7];

[g,e] = shurInv(r);



