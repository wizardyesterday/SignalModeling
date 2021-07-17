//**********************************************************************

//  Name: autocorrelate
//
//  Purpose: The purpose of this function is to compute the
//  estimated autocorrelation function for an input sequence,
//  Rather than performing the correlation in the time domain, the
//  estimate is computed by the following steps:
//
//    1. Zero pad the input vector to the desired length.
//    2. Compute the Fourier transform of the input vector.
//    3. Compute the power spectral density of the frequency data.
//    4. Compute the inverse Fourier transform of the power spectral
//    density.  This is the autocorrelation function.
//
//  Calling Sequence: r = autocorrelate(x,vectorLength)
//
//  Inputs:
//
//    x - The input vector to be processed.
//
//    vectorLength - The adjusted length after zeros are appended to
//    the input vector before processing occurs.
//
//  Outputs:
//
//    r - The estimated autocorrelation function.
//
//**********************************************************************
function r = autocorrelate(x,vectorLength)

// Compute number of zeros to pad to the input vector.
paddingLength = vectorLength - length(x);

if paddingLength < 0
  // Handle this corner case.
  paddingLength = 0;
end

// Perform the zero padding.
y = [x zeros(1,paddingLength)];

// Transform to the frequency domain.
Y = fft(y,-1);

// Compute power spectral density
P = Y .* conj(Y);

// Compute autocorrelation estimate via inverse FFT.
r = fft(P,1);

// Deal with the scrambling due to the FFT.
r = fftshift(r);

endfunction

//**********************************************************************
// Mainline code.
//**********************************************************************

// Allocate vector.
x = ones(1,8);

// Compute actual autocorrelation function.
r = convol(x,x($:-1:1));

// Compute autocorrelation estimates.
r8 = autocorrelate(x,8);
r16 = autocorrelate(x,16);
r32 = autocorrelate(x,32);

