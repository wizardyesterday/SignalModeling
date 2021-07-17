//**********************************************************************
// File Name: C3_4.sci
//**********************************************************************

//**********************************************************************
// Mainline code.
//**********************************************************************

//****************************************
// Part (a)
//****************************************
// Generate 10000 samples of Gaussian noise with unity variance.
noisegen(1,10000,1);

// Construct noise source.
v = feval([1:24],Noise);

// Compute autocorrelation of noise source.
rv = convol(v(1:24),v(24:-1:1));
indexMax = find(rv == max(rv));
rv = rv(indexMax:$);

// Set initial values of x(n).
x(1) = v(1);
x(2) = v(2);

// Compute output of filter when driven by a noise source.
for i = 3:24
  x(i) = v(i) - 0.81*x(i-2);
end

//****************************************
// Part (b)
//****************************************
// Construct impulse response of the filter.
n = 0:23;
h = (0.81)^(n/2) .* cos(%pi*n/2);

// Compute sample autocorrelation.
rxhat = convol(x,x($:-1:1));
indexMax = find(rxhat == max(rxhat));
rxhat = rxhat(indexMax:$);

// Compute autocorrelation of filter
rh = convol(h,h($:-1:1));
indexMax = find(rh == max(rh));
rh = rh(indexMax:$);

// Since the true value of rv is delta(n), the true autocorrelation is trivial.
rx = rh;

//****************************************
// Part (c)
//****************************************
// Compute power spectral density.
Pxhat = fft(rxhat,-1);

//****************************************
// Part (d)
//****************************************
r0 = rxhat(1);
r1 = rxhat(2);
r2 = rxhat(3);

// Construct autocorrelation matrix
R = toeplitz([r0 r1]);

// Construct constant vector.
r = -[r1 r2]';

// Compute estimated filter coefficients.
a = inv(R) * r;

a1 = a(1);
a2 = a(2);
b0 = sqrt(r0 + r1*a1 + r2*a2);

//****************************************
// Part (e)
//****************************************
// Set initial values of xe(n).
xe(1) = b0;
xe(2) = a1*b0;

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Compute impulse response of estimated filter.
// Note that the Yule-Walker equations are derived from
// the form (second order example),
// x(k) + a(1)x(n-1) + a(2)x(n-2) = b(0)delta(k),
// therefore, when we generate an impulse response,
// we write x(k) = -a(1)x(n-1) - a(2)x(n-2) + b(0)delta(k).
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
for i = 3:24
  xe(i) = -a1*xe(i-1) - a2*xe(i-2);
end

// Compute Fourier transform.
Xe = fft(xe,-1);

// Compute power spectral density.
Pxe = Xe .* conj(Xe);

//****************************************
// Part (f)
//****************************************
// Compute transfer function.
H = fft(h,-1);

// Compute output power spectral density.
Px = H .* conj(H);

//****************************************
// Create plots.
//****************************************
// Display first figure.
scf(1);

subplot(311);
title("White Noise Autocorrelation, rv(k)");
plot(rv);

subplot(312);
title("Sample Autocorrelation, rxhat(k)");
plot(rxhat);

subplot(313);
title("True Autocorrelation, rx(k)");
plot(rx);

// Display second figure.
scf(2);

subplot(311);
title("Sample Power Spectral Density, Pxhat");
plot(Pxhat);

subplot(312);
title("Power Spectral Density of Estimated Filter");
plot(Pxe);

subplot(313);
title("True Power Spectral Density, Px");
plot(Px);

