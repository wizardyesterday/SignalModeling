//**********************************************************************
// File Name: C5_4.sci
//**********************************************************************

// Bring in all of the utility functions.
exec('utils.sci',-1);

//**********************************************************************
// Mainline code.
//**********************************************************************

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Simple valid data case.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Start with a valid autocorrelation sequenc.
r = [4 1 1 -2]';

// Generate the associated reflection coefficient sequence.
[g,e] = rtog(r);

// Compute the split Levinson coefficient sequence.
d1 = gtod(g);

// Recover the reflection coefficient sequence.
g1 = dtog(d1);

disp("[g d1 g1] =");
disp([g d1 g1]);

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Generate boundary cases for |g(j) = 1.
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Construct a reflection coefficient sequence with g(p) = 1.
g2 = [.1 .5 1]';

// See what happens to the split Levinson coefficient sequence.
d2 = gtod(g2);

disp("[g2 d2] =");
disp([g2 d2]);

// Construct a reflection coefficience with g(p) = -1.
g3 = [.1 .5 -1]';

// See what happens to the split Levinson coefficient sequence.
d3 = gtod(g3);

disp("[g3 d3] =");
disp([g3 d3])

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Generate boundary cases for |g(j) = 1, -1,...
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
g4 = [1 -1 1 -1]';

d4 = gtod(g4);

disp("[g4 d4] =");
disp([g4 d4]);
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
// Generate boundary cases for |g(j) = 1, 0, 1` ,0,...
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
g5 = [1 0 1 0]';

d5 = gtod(g5);

disp("[g5 d5] =");
disp([g5 d5]);
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/








