//**********************************************************************
// File Name: LevinsonRecursion.sci
//**********************************************************************

//**********************************************************************
//
//  Name: rtoa
//
///  The purpose of this function is to solve the Toeplitz normal equations,
//
//  R * a = epsilon [1 0 ... 0]'
//
//  where R=toeplitz(r) is a Toeplitz matrix that contains the
//  autocorrelation sequence r(k).  The ouputs are the coefficients of
//   the all-pole model a(k) and the constantepsilon.

//  Calling Sequence: [a,epsilon] = rtoa(r)
//
//  Inputs:
//
//    r - The autocorrelation sequence for which the model parameters
//    are to be constructed.
//
//  Outputs:
//
//    a - The denominator coefficients in the model.
//
//    epsilon - The mean square error eP.
//
//**********************************************************************
function [a,epsilon] = rtoa(r)

  // Ensure that we have a column vector.
  r = r(:);

  p = length(r) - 1;

  // a0(0) = 1.
  a = 1;

  // e0 = rx(0).
  epsilon = r(1);

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Construct coefficient sequence.
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  for j = 2:p+1
    gamm = -r(2:j)' * flipud(a) / epsilon;
    a=[a;0] + gamm * [0; conj(flipud(a))];
    epsilon = epsilon * (1 - abs(gamm)^2);
  end
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

endfunction

//**********************************************************************
//
//  Name: rtog
//
//  The purpose of this function is to construct the reflection
//  coefficients from the model coefficients.

//  Calling Sequence: [gamm,epsilon] = rtog(r)
//
//  Inputs:
//
//    r - The autocorrelation sequence for which the reflection
//    coefficients are to be constructed.
//
//  Outputs:
//
//    gamm - The reflection coefficients.
//
//    epsilon - The mean square error.
//
//**********************************************************************
function [gamm,epsilon] = rtog(r)

  // Construct model coefficients.
  [a,epsilon] = rtoa(r);

  // Construct reflection coefficients.
  gamm = atog(a);

endfunction

//**********************************************************************
//
//  Name: gtoa
//
//  The purpose of this function is to construct the direct-form
//  filter coefficients,
//
//  a=[1 a(1) ... a(p)],
//
//  from the reflection coefficients using the step-up recursion.
//

//  Calling Sequence: a = gtoa(gamm)
//
//  Inputs:
//
//    gamm - The reflection coefficients.
//
//  Outputs:
//
//    a - The model coefficients.
//
//**********************************************************************
function a=gtoa(gamm)

  // Ensure that we have a column vector.
  gamm = gamm(:);

  // a0(0) = 1.
  a = 1;

  p = length(gamm);

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Construct coefficient sequence.
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  for j = 2:p+1
     a = [a;0] + gamm(j-1) * [0; conj(flipud(a))];
 end
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

endfunction

//**********************************************************************
//
//  Name: atog
//
//  The purpose of this function is to construct the reflection
//  coefficients from the direct-form filter coefficients, a(k).
//  filter coefficients, using the step-down recursion.
//
//  Calling Sequence: a = gamm = atog(a)
//
//  Inputs:
//
//    a - The model coefficients.
//
//  Outputs:
//
//    gamm - The reflection coefficients.
//
//**********************************************************************
function gamm = atog(a)

  // Ensure that we have a column vector.
  a = a(:);

  p = length(a);

  a = a(2:p) / a(1);

  gamm(p-1) = a(p-1);
  
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Construct reflection coefficients.
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  for j = p-1:-1:2
    a= (a(1:j-1) - gamm(j) * flipud(conj(a(1:j-1)))) ./ (1 - abs(gamm(j))^2);
    gamm(j-1) = a(j-1);
  end
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

endfunction

//**********************************************************************
//
//  Name: gtor
//
//  The purpose of this function is to construct the autocorrelation
//  sequence from the reflection coefficients.  If the optional input
//  epsilon is omitted, then the autocorrelation sequence is
//  normalized so that r(0)=1.
//
//  Calling Sequence: a = gamm = atog(a)
//
//  Inputs:
//
//    gamm - The reflection coefficients.
//
//  Outputs:
//
//    r - The autocorrelation sequence.
//
//**********************************************************************
function r = gtor(gamm,epsilon)

  p = length(gamm);

  aa = gamm(1);

  r = [1 -gamm(1)];

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Construct autocorrelation sequence.
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  for j = 2:p
    aa = [aa; 0] + gamm(j) * [conj(flipud(aa));1];
    r = [r -fliplr(r) * aa];
  end
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // argn(2) returns is the number of arguments passed to the
  // function.
  // If 2 arguments were  passed to the function, it is
  // implied that the last parameter was passed.
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  if argn(2) == 2
    // Normalize.
    r = r * epsilon/prod(1-abs(gamm).^2);
  end

endfunction
      
//**********************************************************************
//
//  Name: ator
//
//  The purpose of this function is to construct the autocorrelation
//  sequence from the model coefficients.  If the optional input
//  epsilon is omitted, then the autocorrelation sequence is
//  normalized so that r(0) = 1.
//
//  Calling Sequence: r = ator(a,b)
//
//  Inputs:
//
//    a - The denominator coefficients of the model.
//
//    b - The numerator coefficient of the model.
//
//  Outputs:
//
//    r - The autocorrelation sequence.
//
//**********************************************************************
function r = ator(a,b)

  p = length(a) - 1;

  // Construct reflection coefficients.
  gamm = atog(a);

  // Construction autocorrelation sequence.
  r = gtor(gamm);

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // argn(2) returns is the number of arguments passed to the
  // function.
  // If 2 arguments were  passed to the function, it is
  // implied that the last parameter was passed.
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  if argn(2) == 2
    // Normalize.
    r = r * sqrt(b) / prod(1 - abs(gamm).^2);
  end

endfunction

//**********************************************************************
//
//  Name: glev
//
//  The purpose of this function is to solve the Toeplitz equations,
//  R*x = b, where R = toeplitz(r), and b is an arbitrary vector.  This
//  function performs the general Levinson recursion.
//
//  Calling Sequence: x = glev(r,b)
//
//  Inputs:
//
//    r - An autocorrelation sequence.
//
//    b - An arbitrary vector.
//
//  Outputs:
//
//    x - The computed sequence.
//
//**********************************************************************
function x = glev(r,b)

  // Force a column vector.
  r = r(:);

  p = length(b);

  // a0(0) = 1.
  a = 1;

  // x0(0) = b(0) / rx(0).
  x = b(1) /  r(1);

  // e0 = rx(0).
  epsilon = r(1);

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Construct data sequence.
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  for j = 2:p
    g = r(2:j)' * flipud(a)
    gamm(j-1) = -g / epsilon
    a = [a; 0] + gamm(j-1) * [0; conj(flipud(a))]
    epsilon = epsilon * (1 - abs(gamm(j-1))^2)
    delta =r(2:j)' * flipud(x)
    q = (b(j) - delta) / epsilon
    x = [x; 0] + q * [conj(flipud(a))]
  end
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

endfunction

//**********************************************************************
//
//  Name: shur
//
//  The purpose of this function is to apply the Shur recursion to a
//  sequence of autocorrelation values and the generate a sequence of
//  reflection coefficients.
//
//  Calling Sequence: [gamm,epsilon] = shur(r)
//
//  Inputs:
//
//    r - The autocorrelation sequence for which the reflection
//    coefficients are to be constructed.
//
//  Outputs:
//
//    gamm - The reflection coefficients.
//
//    epsilon - The mean square error eP.
//
//**********************************************************************
function [gamm,epsilon] = shur(r)

  // Ensure that we have a column vector.
  r = r(:);

  // Upper loop index.
  p = length(r) - 1;

  // g0(k) = r(k).
  g = r;

  // gR(k) = r(k).
  gR = r;

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Construct reflection coefficient
  // sequence.  Note that the indices
  // differ from the algorithm in the
  // book since the book assumes
  // zero-based arrays, whereas
  // Scilab assumes one-based arrays.
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  for j = 0:p-1
    // Update gamma.
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
    //for k = j+2:p
    //  g(k+1) = gPrev(k+1) + gamm(j+1) * gRPrev(k);
    //end

    // Vectorized version of above loop.
    g(j+3:p+1) = gPrev(j+3:p+1) + gamm(j+1) .* gRPrev(j+2:p);

    //for k = j+1:p
    //  gR(k+1) = gRPrev(k) + conj(gamm(j+1)) * gPrev(k+1);
    //end

    // Vectorized version of above loop.
    gR(j+2:p+1) = gRPrev(j+1:p) + conj(gamm(j+1)) .* gPrev(j+2:p+1);
    //----------------------------------------------------
  end
 //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

  // eP = gP^(R)(p).
  epsilon = gR(p+1);
 
endfunction

//**********************************************************************
//
//  Name: gtod
//
//  Purpose: The purpose of this function is to compute the split
//  Levinson coefficient sequence given the reflection coefficient
//  sequence.
//
//  Calling Sequence: d = gtod(g)
//
//  Inputs:
//
//    g - The reflection coefficient sequence.
//
//
//  Outputs:
//
//    d - The split Levinson coefficient sequence.
//
//**********************************************************************
function d = gtod(g)

  // Force a column vector.
  g = g(:);

  // Compute the number of coefficients
  p = length(g);

  // Set initial d(j) assuming gamma(0) = 0.
  d(1) = 1 - g(1);

  // Construct split Levenson sequence.
  for j = 2:p
    d(j) = (1 - g(j))*(1 + g(j-1));
  end
 
endfunction

//**********************************************************************
//
//  Name: dtog
//
//  Purpose: The purpose of this function is to compute the reflection
//  coefficient sequence given the split Levinson coefficient sequence.
///
//  Calling Sequence: g = dtog(d)
//
//  Inputs:
//
//    d - The split Levinson coefficient sequence.
//
//  Outputs:
//
//    g - The reflection coefficient sequence.
//
//**********************************************************************
function g = dtog(d)

  // Force a column vector.
  d = d(:);

  // Compute the number of coefficients
  p = length(d);

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Search for '-1' values in g, and adjust the values,
  // This avoids division by zero later on.
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  n = find(g == -1);
  g(n) = g(n) + 1e-6;
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

  // Extend the sequence to insert gamma(0) = 0.
  g = [0; g];

  // Construct reflection coefficient sequence.
  for j = 2:p+1
    g(j) = 1 - d(j-1) / (1 + g(j-1));
  end

  // Remove the initial condition.
  g(1) = [];
 
endfunction

//**********************************************************************
//
//  Name: slev
//
//  The purpose of this function is to solve for the linear predictor
//  coefficients, given the autocorrelation sequence, using the
//  split Levinson recursion algorithm.
//
//  Calling Sequence: a = slev(r)
//
//  Inputs:
//
//    r - An autocorrelation sequence.
//
//
//  Outputs:
//
//    a - The linear predictor coefficients.
//
//**********************************************************************
function a = slev(r)

  // Force a column vector.
  r = r(:);

  p = length(r) - 1;

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Preallocate vectors.
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  tau = zeros(1,p+1);

  sjm1 = zeros(1,p+1)';
  sj = zeros(1,p+1)';
  sjp1 = zeros(1,p+1)';

  a = zeros(1,p)';
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Set initial conditions.
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  //-------------------------------------
  // Initialize pipeline.
  //-------------------------------------
  // s0 = 2.
  sjm1(1) = 2;

  // s1 = [1 1]';
  sj(1) = 1;
  sj(2) = 1;
  //-------------------------------------

  // Note that tauJ = epsilonJ / (1 + gammaJ).
  tau(1) = r(1);

  // Initialize the first linear predictor coefficient.
  a(1) = 1;
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

  for j = 1:p

    // Generate tau
    for i = 0:j
      tau(j+1) = tau(j+1) + r(i+1)*sj(i+1);
    end

    delta = tau(j+1) / tau(j);

    // Initialize.
    sjp1(1) = 1;
    sjp1(j+2) = 1;

    // Compute the next singular predictor coefficients.    
    for k = 1:j
      sjp1(k+1) = sj(k+1) + sj(k) - delta*sjm1(k);
    end

    // Update the pipeline for the next recursion.
    sjm1 = sj;
    sj = sjp1;
  end

  onePlusGammaP = sum(sj) / sum(sjm1);

  // Compute the linear predictor coefficients.
  for j = 1:p
    a(j+1) = a(j) + sjp1(j+1) - onePlusGammaP*sjm1(j);
  end

endfunction
