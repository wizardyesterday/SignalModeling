//**********************************************************************
// File Name: C6_5.sci
//**********************************************************************

exec('utils.sci',-1);

//**********************************************************************
// Mainline code.
//**********************************************************************
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (a), generate a second order model
// for x(n) = (0.8)^n + (1.25)^n for
// n = 0,1,...,15 using the modified
// covariance method.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Generate the time vector.
n = 0:15;

// Set the parameters.
alpha1 = 0.8;
beta1 = 1.25;

// Generate the data sequence.
x1 = alpha1.^n + beta1.^n;

// Generate the model coefficients and reflection coefficients.
[aaMcov,eaMcov] = mcov(x1,2);
gaMcov = atog(aaMcov);

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (b), generate second orders model
// for x(n) = 0.8^n + 1.25^n for
// n = 0,1,...,15 using the Burg's method,
// the covariance method, and the
// autocorrelation method.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Generate the model coefficients and reflection coefficients.
[gbBurg,ebBurg] = burg(x1,2);
abBurg = gtoa(gbBurg);

// Generate the model coefficients and reflection coefficients.
[abCovm,ebCovm] = covm(x1,2);
gbCovm = atog(abCovm);

// Generate the model coefficients and reflection coefficients.
[abAcm,ebAcm] = acm(x1,2);
gbAcm = atog(abAcm);

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Part (c), repeat parts (a) and (b) with
// x(n) = 0.75^n + 2^n for n = 0,1,...,15.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Set the parameters.
alpha2 = 0.75;
beta2 = 2;

// Generate the data sequence.
x2 = alpha2.^n + beta2.^n;

// Generate the model coefficients and reflection coefficients.
[acaMcov,ecaMcov] = mcov(x2,2);
gcaMcov = atog(acaMcov);

// Generate the model coefficients and reflection coefficients.
[gcbBurg,ecbBurg] = burg(x2,2);
acbBurg = gtoa(gcbBurg);

// Generate the model coefficients and reflection coefficients.
[acbCovm,ecbCovm] = covm(x2,2);
gcbCovm = atog(acbCovm);

// Generate the model coefficients and reflection coefficients.
[acbAcm,ecbAcm] = acm(x2,2);
gcbAcm = atog(acbAcm);

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Print results.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
printf("\n--------------------------------------------------------\n");
printf("          x(n) = %f^n + %f^n\n",alpha1,beta1);
printf("--------------------------------------------------------\n");
printf("Mcov e: %3.2e  Burg e: %3.2e  Covm e: %3.2e  Acm E:%3.2e\n", ...
        eaMcov($),ebBurg($),ebCovm($),ebAcm($));

printf("Mcov g: %3.2e  Burg g: %3.2e Covm g: %3.2e  Acm g:%3.2e\n", ...
        gaMcov,gbBurg,gbCovm,gbAcm);

printf("\n--------------------------------------------------------\n");
printf("          x(n) = %f^n + %f^n\n",alpha1,beta1);
printf("--------------------------------------------------------\n");
printf("Mcov e: %3.2e  Burg e: %3.2e  Covm e: %3.2e  Acm E:%3.2e\n", ...
        ecaMcov($),ecbBurg($),ecbCovm($),ecbAcm($));

printf("Mcov g: %3.2e  Burg g: %3.2e Covm g: %3.2e  Acm g:%3.2e\n", ...
        gcaMcov,gcbBurg,gcbCovm,gcbAcm);

