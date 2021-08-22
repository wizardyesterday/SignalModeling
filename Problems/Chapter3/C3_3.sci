//**********************************************************************

//  Name: autocorrelate
//
//  Purpose: The purpose of this function is to compute the
//  estimated autocorrelation function for an input sequence.
//
//  Calling Sequence: r = autocorrelate(x,lag)
//
//  Inputs:
//
//    x - The input vector to be processed.
//
//    lag - The autocorrelation lag.
//
//  Outputs:
//
//    r - The estimated autocorrelation function.
//
//**********************************************************************
function r = autocorrelate(x,lag);

// Start with a clean slate.
accumulator = 0;

// Avoid constant reevaluation during loop execution.
upperLimit = length(x) - lag;

for i = 1:upperLimit
  accumulator = accumulator + x(i) * x(i + lag);
end

// Scale as appropriate.
r = accumulator / length(x);

endfunction

//**********************************************************************

//  Name: segmentedAutocorrelate
//
//  Purpose: The purpose of this function is to compute the
//  estimated autocorrelation function for an input sequence by
//  breaking up a 1000 element input vector into 100 element segments,
//  running an autocorrelation on each segment, and averaging the ten
//  segments.  It should be noted that advantage of the colon operator
//  of Scilab has been utilized to simply processing.
//
//  Calling Sequence: r = segmentedAutocorrelate(x)
//
//  Inputs:
//
//    x - The input vector to be processed.
//
//  Outputs:
//
//    r - The estimated autocorrelation function.
//
//**********************************************************************
function [r] = segmentedAutocorrelate(x)

  // Preallocate the temporary matrix.
  a = zeros(10,199);

  // Reference the first segment of the input vector.
  startIndex = 1;

  for i = 1:10

    // Grab the current segment of the input vector.
    y = x(startIndex:startIndex + 99);

    // Perform the autocorrelation.
    y = convol(y,y($:-1:1));

    // Fill in the ith row of the temporary matrix.
    a(i,1:$) = y;

    // Reference the next segment of the input vector.
    startIndex = startIndex + 100;

  end

  // Compute the sum of the columns of a.
  r = sum(a,'m');

  // Divide by the number of rows, since we're computing an average.
  r = r / 1000;

endfunction


//**********************************************************************
// Mainline code.
//**********************************************************************

// Generate 10000 samples of Gaussian noise with unity variance.
noisegen(1,10000,1);
a1000 = 1:1000;
a10000 = 1:10000;
y1000 = feval(a1000,Noise);
y10000 = feval(a10000,Noise);

// Preallocate output vectors.
r1000 = zeros(1,100);
r10000 = zeros(1,100);

// Compute autocorrelations for the small and large data sets.
for i = 0:99
  r1000(i+1) = autocorrelate(y1000,i);
  r10000(i+1) = autocorrelate(y10000,i);
end

rSegmented = segmentedAutocorrelate(y1000);


