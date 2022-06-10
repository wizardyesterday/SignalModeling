//**********************************************************************
// File Name: C8_5.sci
//**********************************************************************

// Bring in all of the utility functions.
exec('utils.sci',-1);

//**********************************************************************
// Mainline code.
//**********************************************************************
// Generate white noise sequence with unit variance.
noisegen(1,256,1);
w = feval([1:256],Noise);
w = w(:);

// Create the filter for the wideband process.
// Create the filter.
a1 = [1 -0.5 0.5];
a2 = [1 0 0.5];
awb = convol(a1,a2);
awb = awb(:);

// Create the filter for the narrowband process.
// Create the filter.
a1 = [1 -1.585 0.96];
a2 = [1 -1.152 0.96];
anb = convol(a1,a2);
anb = anb(:);

// Generate wideband process.
xwb = filterBlock(w,1,awb(2:$));

// Generate narrowband process.
xnb = filterBlock(w,1,anb(2:$));



