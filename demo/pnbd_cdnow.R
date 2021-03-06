## Authors: Lukasz Dziurzynski, Daniel McCarthy, Edward Wadsworth

data(cdnowSummary)

## Get the calibration period customer-by-sufficient-statistic matrix from the cdnow data:
cbs <- cdnowSummary$cbs

## Estimate parameters for the Pareto/NBD model from the CBS:
par.start <- c(0.5, 1, 0.5, 1)
params <- pnbd.EstimateParameters(cbs, 
                                  par.start, 
                                  hardie = TRUE)
params

## Check log-likelihood of the params:
pnbd.cbs.LL(params, 
            cbs, 
            hardie = TRUE)

## Plot the comparison of actual and expected calibration period frequencies:
pnbd.PlotFrequencyInCalibration(params = params, 
                                cal.cbs = cbs, 
                                censor = 7, 
                                plotZero = TRUE, 
                                hardie = TRUE)

T.star <- 39                # Length of holdout period
x.star <- cbs[,"x.star"]    # Transactions made by each customer in the holdout period

## Plot the comparison of actual and conditional expected holdout period frequencies,
## binned according to calibration period frequencies:
pnbd.PlotFreqVsConditionalExpectedFrequency(params = params, 
                                            T.star = T.star, 
                                            cal.cbs = cbs, 
                                            x.star = x.star, 
                                            censor = 7, 
                                            hardie = TRUE)

## Plot the comparison of actual and conditional expected holdout period frequencies,
## binned according to calibration period recencies:
pnbd.PlotRecVsConditionalExpectedFrequency(params, 
                                           cbs, 
                                           T.star, 
                                           x.star, 
                                           hardie = TRUE)
