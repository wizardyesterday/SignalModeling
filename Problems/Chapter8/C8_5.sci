//**********************************************************************
// File Name: C8_5.sci
//**********************************************************************

// Bring in all of the utility functions.
exec('utils.sci',-1);

//**********************************************************************
//
//  Name: evaluateFilterOrders
//
//  The purpose of this function is to perform the drudgery
//  evaluating various error criteria for determining the filter
//  order to model a given random process.
//
//  Calling Sequence: evaluateFilterOrders(x,N)
//
//  Inputs:
//
//    x - The input sequency that represents a random process.
//
//    N - The data record length of the input sequence being
//    evaluated.
//
//  Outputs:
//
//    None.
//
//**********************************************************************
function evaluateFilterOrders(x,N,functionPtr);

  // Force a column vector.
  x = x(:);

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Evaluate using Akaike Information Criterion.
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Start with a first-order system.
  pAic = 1;

  // Pick something huge.
  previousC = 1e6;

  // Set up for loop entry.
  done = 0;

  while done == 0
    // Invoke model.
    [model,ep] = functionPtr(x,pAic);

    // Evaluate criteria.
    c = Aic(ep($),pAic,N);

    if previousC > c
      // Increase the model order.
     pAic = pAic + 1;

      // Update the previous
      previousC = c;
    else
      // We hit a minimum.

      if pAic > 1
        // Use previous filter order.
        pAic = pAic - 1;
      end

      // Bail out.
      done = 1;
    end
  end
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Evaluate using minimum description length.
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Start with a first-order system.
  pMdl = 1;

  // Pick something huge.
  previousC = 1e6;

  // Set up for loop entry.
  done = 0;

  // Pick something huge.
  previousC = 1e6;

  // Set up for loop entry.
  done = 0;

  while done == 0
    // Invoke model.
    [model,ep] = functionPtr(x,pMdl);

    // Evaluate criteria.
    c = Mdl(ep($),pMdl,N);

    if previousC > c
      // Increase the model order.
     pMdl = pMdl + 1;

      // Update the previous
      previousC = c;
    else
      // We hit a minimum.

      if pMdl > 1
        // Use previous filter order.
        pMdl = pMdl - 1;
      end

      // Bail out.
      done = 1;
    end
  end
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Evaluate using Akaike's Final Prediction Error.
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Start with a first-order system.
  pFpe = 1;

  // Pick something huge.
  previousC = 1e6;

  // Set up for loop entry.
  done = 0;

  while done == 0
    // Invoke model.
    [model,ep] = functionPtr(x,pFpe);

    // Evaluate criteria.
    c = Fpe(ep($),pFpe,N);

    if previousC > c
      // Increase the model order.
     pFpe = pFpe + 1;

      // Update the previous
      previousC = c;
    else
      // We hit a minimum.

      if pFpe > 1
        // Use previous filter order.
        pFpe = pFpe - 1;
      end

      // Bail out.
      done = 1;
    end
  end
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Evaluate using Parzen's Criterion Autoregressive
  // Transfer.
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // Start with a first-order system.
  pCat = 1;

  // Pick something huge.
  previousC = 1e6;

  // Set up for loop entry.
  done = 0;

  while done == 0
    // Invoke model.
    [model,e] = functionPtr(x,pCat);

    if length(e) == 1
      ep(pCat) = e;
    else
      ep = e;
    end
 
    // Evaluate criteria.
    c = Cat(ep,pCat,N);

    if previousC > c
      // Increase the model order.
     pCat = pCat + 1;

      // Update the previous
      previousC = c;
    else
      // We hit a minimum.

      if pCat > 1
        // Use previous filter order.
        pCat = pCat - 1;
      end

      // Bail out.
      done = 1;
    end
  end
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

// Print results.
printf("pAic: %f, pMdl: %f, pFpe: %f, pCat: %f\n", ...
       pAic,pMdl,pFpe,pCat);

endfunction

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

printf("Evaluation of model order for wideband process\n");
printf("Autocorrelation method\n");
evaluateFilterOrders(xwb,256,acm);
printf("Burg method\n");
evaluateFilterOrders(xwb,256,burg);

printf("\nEvaluation of model order for narrowband process\n");
printf("Autocorrelation method\n");
evaluateFilterOrders(xnb,256,acm);
printf("Burg method\n");
evaluateFilterOrders(xnb,256,burg);



