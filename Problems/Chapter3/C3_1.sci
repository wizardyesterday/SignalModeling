//**********************************************************************

//  Name: exponential
//
//  Purpose: The purpose of this function is to generate an
//  exponentially distributed random variable.  An exponentially
//  distributed random variable is described as follows:
//
//  Probability Density Function: f(x) a * exp(-a*x),
//  where 1/a is the mean value.
//
//  Cumulative Distribution Function: F(x) = 1 - exp(-a*x) 
//
//  The cumulative distribution function is generated from a
//  unformly distributed random variable, u, as follows:
//
//  x = -(1/a) * ln(1 - u)
//
//  Calling Sequence: x = exponential(mu)
//
//  Inputs:
//
//    mu - The mean value of the exponential distribution.
//
//  Outputs:
//
//    x - A uniformly distributed random variable with mean mu.
//
//**********************************************************************
function x = exponential(mu)

// Generate uniformly distributed random variable.
u = rand();

// Generate argument.
a = 1 - u;

// Generate exponentially distributed ransom number.
x = -mu * log(a);

endfunction

//**********************************************************************
// Mainline code.
//**********************************************************************

// Preallocate the vectors.
x = zeros(1,100000);
a = zeros(1,30);

// Perform 30 runs.
for i = 1:30
  for j = 1:102977
    x(j) = exponential(2);
  end
  // Save the mean value for the run.
  a(i) = mean(x);
end

