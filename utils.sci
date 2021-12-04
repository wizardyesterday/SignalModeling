//**********************************************************************
// File Name: utils.sci
//**********************************************************************

//**********************************************************************
//
//  Name: shiftSampleInPipeline
//
//  Purpose: The purpose of this function is to place a new sample
//  into the pipeline referenced by the state parameter.  The format
//  of the pipeline appears below.
//
//      {x(n) x(n-1) x(n-2)...}
//
//  Calling Sequence: state = shiftSampleInPipeline(x,state)
//
//  Inputs:
//
//    x - The input sample to be shifted into the pipeline.
//
//    state - The pipeline that is typically used to store the state
//    of a digital filter.  The order of storage appears belos.
//
//  Outputs:
//
//    state - The updated pipeline.
//
//**********************************************************************
function state = shiftSampleInPipeline(x,state)

  // Ensure we have a column vector.
  state = state(:);

  // Shift the contents of the pipeline.
  state(2:$) = state(1:$-1);

  // Store sample value.
  state(1) = x;

endfunction

//**********************************************************************
//
//  Name: iirFilter
//
//  Purpose: The purpose of this function is to filter a signal using
//  a recursive filter (in this case an infinite impulse response).
//  Advantage has been taken of specific features of Scilab for the
//  computation of the convolution sums.  Specifically, the for loop
//  has been replaced with the inner product of a pipeline and the
//  corresponding filter coefficients.
//
//  Calling Sequence: [y,bState,aState = iirFilter(x,b,a,bState,aState)
//
//  Inputs:
//
//    x - The input sample.
//
//    b - The numerator (feedforward) coefficients of the filter.
//
//    a - The denominator (feedback) coefficients of the filter.
//
//    bState - The pipeline for the nonrecursive portion of the filter.
//
//    aState - The pipeline for the recursive portion of the filter.
//
//  Outputs:
//
//    y - The output of the filter.
//
//    bState - The pipeline for the nonrecursive portion of
//    the filter updated with the latest sample value.
//
//    aState - The pipeline for the recursive portion of the
//    filter updated with the latest output value.
//
//**********************************************************************
function [y,bState,aState] = iirFilter(x,b,a,bState,aState)

  // Ensure that we have column vectors.
  a = a(:);
  b = b(:);
  aState = aState(:);
  bState = bState(:);

  // Shift input sample into the forward pipeline.
  bState = shiftSampleInPipeline(x,bState);

 // Perform convolution sum for the nonrecursive portion of the filter.
  y = b'*bState;

  // Perform convolution sum for the recursive portion of the filter.
  y = y - a'*aState;

  // Shift output sample into the recursive pipeline.
  aState = shiftSampleInPipeline(y,aState);

endfunction

//**********************************************************************
//
//  Name: filterBlock
//
//  Purpose: The purpose of this function is to filter a signal using
//  a recursive filter.  The function accepts a block of input data
//  and generates a block of output data.  The pipelines are created in
//  this function since it only processes one block of data.
//
//  Calling Sequence: y = filterBlock(x,b,a)
//
//  Inputs:
//
//    x - The input vector.
//
//    b - The numerator coefficients of the filter.
//
//    a - The denominator coefficients of the filter.
//
//  Outputs:
//
//    y - The output vector.
//
//**********************************************************************
function y = filterBlock(x,b,a)

  // Instantiate the pipelines.
  aState = zeros(1,length(a));
  bState = zeros(1,length(b));

  // Compute the filter response on a sample by sample basis.
  for i = 1:length(x)
    [y(i),bState,aState] = iirFilter(x(i),b,a,bState,aState);
  end

endfunction

//**********************************************************************
//
//  Name: convm
//
//  Purpose: The purpose of this function is to construct a convolution
//  matrix.
//
//  Calling Sequence: X = convm(x,p)
//
//  Inputs:
//
//    x - The input vector to be processed.
//
//    p - The desired number of columns in the matrix.
//
//  Outputs:
//
//    X - The convolution matrix.
//
//**********************************************************************
function X = convm(x,p)

  N = length(x) + 2 * p - 2;

  // Ensure that we have a column vector.
  x = x(:);

  // Construct the appropriate zero padding.
  xpad = [zeros(p-1,1); x; zeros(p-1,1)];

  // Construct the matrix.
  for i = 1:p
    X(:,i) = xpad(p-i+1:N-i+1);
  end

endfunction

//**********************************************************************
//
//  Name: flipud
//
//  Purpose: The purpose of this function is to reverse the order of
//  the entries of a vector.
//
//  Calling Sequence: b = flipud(a)
//
//  Inputs:
//
//    a - The input vector.
//
//  Outputs:
//
//    b - The vector that contains reversed entries of a.
//
//**********************************************************************
function b = flipud(a)

  // Reverse the order of a.
  b = a($:-1:1);

endfunction

//**********************************************************************
//
//  Name: fliplr
//
//  Purpose: The purpose of this function is to reverse the order of
//  the entries of a vector.
//
//  Calling Sequence: b = fliplr(a)
//
//  Inputs:
//
//    a - The input vector.
//
//  Outputs:
//
//    b - The vector that contains reversed entries of a.
//
//**********************************************************************
function b = fliplr(a)

  // Reverse the order of a.
  b = a($:-1:1);

endfunction

//**********************************************************************
//
//  Name: quantize
//
//  Purpose: The purpose of this function is to quantize a fractional
//  number to a finite number of integer bits.
//
//  Calling Sequence: aQ = quantize(a,b)
//
//  Inputs:
//
//    a - The input vector.
//
//    b - The number integer bits of precision.
//
//  Outputs:
//
//    aQ - The quantized output vector, if the quantization was
//    successful, otherwise, the original vector is returned.
//
//    success - An indicator of the outcome of the operation.  A value
//    of 1 implies that quantization had occurred, and a value of 0
//    implies that the quantization has not occurred.
//
//**********************************************************************
function [aQ,success] = quantize(a,b)

  // Force a column vector.
  a = a(:)

  // Default to failure.
  success = 0;

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Search for any numbers whose magnitude matches or
  // exceeds unity.  Due to the fact that this function
  // processes fractional numbers, if numbers greater
  // than or equal to 1 in magnitude, the vector, a', is
  // not a candidate for conversion.
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  indexOfBigNumbers = find(abs(a) >= 1);
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

  if indexOfBigNumbers == [];
    positiveLimit = 2^(b - 1) - 1;
    negativeLimit = 2^(b - 1);

    // Locate negative values.
    n = find(a < 0);

    // Compute the quantized value.
    aQ = a * positiveLimit;
    aQ(n) = a(n) * negativeLimit;

    // Remove any fractional part.
    aQ = round(aQ);

    // Indicate a successful quantization.
    success = 1;
  else
    // Handle the zero vector case.
    aQ = a;
  end

endfunction


exec('SignalModeling.sci',-1);
exec('LevinsonRecursion.sci',-1);
exec('Lsp.sci',-1);

