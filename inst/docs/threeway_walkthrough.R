# Let's check BTYD vs. BTYD2 vs. Theo's fix
library('dplyr')
library('lubridate')
library('tidyr')
library('purrr')
library('ggplot2')

# BTYD canonical data, CD-Now set: people who made
# their first purchase during q1 of 1997, observed 
# through mid 1998.
cdnowElog <- system.file("data/cdnowElog.csv", package = "BTYD")
elog      <- BTYD::dc.ReadLines(cdnowElog, 
                                cust.idx = 2, date.idx = 3, sales.idx = 5)
elog$date <- as.Date(elog$date, "%Y%m%d")

# Here's how I know they made 1st purchase in q1 1997:
foo <- tbl_df(elog) %>% 
  group_by(cust) %>% 
  filter(date == min(date))
table(year(foo$date), month(foo$date))
rm(foo)

# And here's how I know they were seen through mid-1998:
max(elog$date)

# Consolidate transactions that happened on the same date. 
elog <- BTYD::dc.MergeTransactionsOnSameDate(elog)

# pick date that falls in the middle of the 
# time span covered by the data (this won't
# be the same as half the observations, but 
# this is what the vignette means)
getHalfTimeDate <- function(rfm, floor = TRUE) {
  sort(unique(rfm$date))[floor(length(unique(rfm$date))/2)]
  if(floor == FALSE) {
    sort(unique(rfm$date))[ceiling(length(unique(rfm$date))/2)]
  }
}
elog.vig <- getHalfTimeDate(elog, floor = FALSE)

# Common core of stuff that needs to be done the same way.
# These functions are unchanged in BTYD2.
replicateBTYDhelper <- function(df, enddate) {
  end.of.cal.period <- enddate
  elog              <- df
  elog.cal          <- elog[which(elog$date <= end.of.cal.period), ]
  
  split.data        <- BTYD::dc.SplitUpElogForRepeatTrans(elog.cal)
  clean.elog        <- split.data$repeat.trans.elog
  freq.cbt          <- BTYD::dc.CreateFreqCBT(clean.elog)
  tot.cbt           <- BTYD::dc.CreateFreqCBT(elog)
  cal.cbt           <- BTYD::dc.MergeCustomers(tot.cbt, freq.cbt)
  birth.periods     <- split.data$cust.data$birth.per
  last.dates        <- split.data$cust.data$last.date
  cal.cbs.dates     <- data.frame(birth.periods, last.dates, 
                                  end.of.cal.period)
  cal.cbs           <- BTYD::dc.BuildCBSFromCBTAndDates(cal.cbt, 
                                                        cal.cbs.dates, 
                                                        per="week")
  cal.cbs
}

# Make use of dplyr facilities to do transformations
# (quicker, easier to read, works). 
# You have:
# -- an order log, df, with purchases consolidated by day
# -- an end of window date, enddate.
# You want:
# -- in the calibration period (cal = TRUE):
#    -- x, the number of transactions after the first (may be 0)
#    -- t.x, the time of the last transaction after the first (ditto)
#    -- T.cal, the total observation time btw 1st trans and enddate.
#    -- if compress = TRUE, also a column 'custs' that counts customers 
#       who have this specific combo of (x, t.x, T.cal). in other words, 
#       this gets you fewer rows, one extra column. compress = FALSE 
#       would give you 3 columns and as many rows as unique customers in df.
# -- in the holdout period (cal = FALSE):
#    -- x.star, number of expected transactions after the holdout period starts
#    -- T.star, duration of the calibration period
dplyrBTYDhelper <- function(df, 
                            enddate, 
                            compress = FALSE, 
                            rounding = 3, 
                            cal = TRUE) {
  # dplyr idiom for replicateBTYDhelper()
  # without any calls to dc. functions.
  out <- tbl_df(df) %>% 
    arrange(cust, date) %>% 
    mutate(calper = date <= enddate) 
  latest <- max(out$date)
  
  # _repeat_ transactions of all customers  
  # in the calibration period. for customers 
  # who only bought once, that number is 0.
  cal.cbs <- filter(out, calper == TRUE) %>%
    mutate(calper = NULL) %>% 
    group_by(cust) %>%
    summarise(x = n() - 1, 
              t.x = as.numeric(difftime(max(date), 
                                        min(date), 
                                        units = 'weeks')), 
              T.cal = as.numeric(difftime(enddate, 
                                          min(date), 
                                          units = 'weeks'))) %>%
    ungroup()
  
  # transactions of any customers seen in the holdout period
  holdout.cbs <- filter(out, calper == FALSE) %>%
    mutate(calper = NULL) %>% 
    group_by(cust) %>% 
    summarise(x.star = n()) %>%
    mutate(T.star = as.numeric(difftime(latest, 
                                        enddate, 
                                        units = 'weeks'))) %>%
    ungroup()
  
  # add zeroes for customers who did not make transactions
  # during the holdout period
  zeroes <- anti_join(dplyr::select(cal.cbs, cust), 
                      dplyr::select(holdout.cbs, cust)) %>%
    mutate(x.star = 0, 
           T.star = as.numeric(difftime(latest, 
                                        enddate, 
                                        units = 'weeks')))
  holdout.cbs <- bind_rows(zeroes, holdout.cbs)
  
  if(cal == TRUE) {
    out <- mutate(cal.cbs, cust = NULL)
    # dplyr idiom for BTYD::pnbd.compress.cbs(), optional:
    if(compress == TRUE) {
      out <- group_by(out, T.cal, t.x, x) %>% 
        summarise(custs = n()) %>%
        ungroup() %>%
        select(x, t.x, T.cal, custs)
    }
  } else {
    out <- mutate(holdout.cbs, cust = NULL)
  }
  as.matrix(round(out, rounding))
}

# To replicate results in section 2.2 of BTYD Walkthrough using
# the CRAN package BTYD, my patched version BTYD2, or Theo's script.
# if verbose is TRUE, print the body of the pnbd.LL() function.
# if theo is TRUE, use Theo's script, assumed to be at
# '../GitHub/theofilos_BTYD/pnbd.R'
btyd.cal.cbs <- replicateBTYDhelper(df = elog, enddate = elog.vig)
dplyr.cal.cbs <- dplyrBTYDhelper(df = elog, enddate = elog.vig)

# Set hardie to FALSE if you want Re(hypergeo) 
replicateBTYD <- function(cal.cbs, 
                          mylib = 'BTYD', 
                          verbose = FALSE,
                          theo = FALSE, 
                          hardie = TRUE) {
  # make sooper sure you use the package you 
  # mean to use, declared in the mylib argument
  try(detach("package:BTYD", unload = TRUE), silent = TRUE)
  try(detach("package:BTYD2", unload = TRUE), silent = TRUE)
  library(mylib, character.only = TRUE)
  myenv <- paste('namespace', mylib, sep = ':')
  
  if(theo == TRUE) {
    source('../theofilos_BTYD/pnbd.R')
    myenv <- 'R_GlobalEnv'
    params  <- pnbd.EstimateParameters(cal.cbs)
    LL <- pnbd.cbs.LL(params, cal.cbs)
  } else if(mylib != 'BTYD') {
    params  <- pnbd.EstimateParameters(cal.cbs, hardie = hardie)
    LL <- pnbd.cbs.LL(params, cal.cbs, hardie)
  } else {
    params  <- pnbd.EstimateParameters(cal.cbs)
    LL <- pnbd.cbs.LL(params, cal.cbs)
  }
  # Which likelihood function will be used?
  print(paste('LL function should be in the', myenv, 'environment', sep=' '))
  if(verbose == TRUE) print(pnbd.LL)
  
  
  # check the environment: in both cases should be the mylib namespace,
  # unless you're using Theo's script, in which case it should be the
  # global environment.
  print('Environment of the calling function pnbd.cbs.LL:')
  print(environment(pnbd.cbs.LL))
  print('Environment of the called function pnbd.LL:') 
  print(environment(pnbd.LL))
  
  # now retrieve your results.
  out <- list(cal.cbs, params, LL)
  names(out) <- c('CBS matrix','Parameter estimates','Log likelihood')
  out
}

# Proof that dplyrBTYDhelper() == replicateBTYDhelper():
# foo <- replicateBTYD(cal.cbs = btyd.cal.cbs)
# bar <- replicateBTYD(cal.cbs = dplyr.cal.cbs)
# setdiff(foo$`CBS matrix`, bar$`CBS matrix`)
# round(foo$`Parameter estimates` - bar$`Parameter estimates`, digits = 8)
# round(foo$`Log likelihood` - bar$`Log likelihood`, digits = 8)

# Check convergence: 
# -- x is an object returned by replicateBTYD()
checkConvergence <- function(x, reps = 2, mylib = 'BTYD') {
  try(detach("package:BTYD", unload=TRUE))
  try(detach("package:BTYD2", unload=TRUE))
  library(mylib, character.only = TRUE)
  print(pnbd.LL)  
  
  cal.cbs <- x[['CBS matrix']] 
  params  <- x[['Parameter estimates']]
  LL      <- x[['Log likelihood']]
  p.matrix <- c(params, LL)
  for (i in 1:reps) {
    params <- BTYD::pnbd.EstimateParameters(cal.cbs, params)
    LL <- BTYD::pnbd.cbs.LL(params, cal.cbs)
    p.matrix.row <- c(params, LL)
    p.matrix <- rbind(p.matrix, p.matrix.row)
  }
  colnames(p.matrix) <- c("r", "alpha", "s", "beta", "LL")
  rownames(p.matrix) <- 1:(reps + 1)
  p.matrix
}

# Plot to compare rate heterogeneity (transaction or dropout)
# between two sets of gamma parameter estimates
plotRH <- function(pars1, pars2, trans = TRUE) {
  xl  <- 'Transaction Rate'
  if(trans == FALSE) xl <- 'Dropout Rate'
  mvl <- function(pars) {
    shape <- pars[1]
    rate  <- pars[2]
    r.mean <- shape/rate
    r.var  <- shape/rate^2
    lim    <- qgamma(0.99, shape = shape, rate = rate)
    return(c(r.mean, r.var, lim))
  }
  mvl1  <- mvl(pars1)
  mvl2  <- mvl(pars2)
  lim1 <- mvl1[3]
  lim2 <- mvl2[3]
  x.axis.ticks <- seq(0, max(lim1, lim2), length.out = 100)
  h1 <- dgamma(x.axis.ticks, shape = pars1[1], rate = pars1[2])
  h2 <- dgamma(x.axis.ticks, shape = pars2[1], rate = pars2[2])
  rate.mean <- paste('(', 
                     paste(round(c(mvl1[1], mvl2[1]), digits=4), collapse = ', '),
                     ')', sep = '')
  rate.var <- paste('(', 
                     paste(round(c(mvl1[2], mvl2[2]), digits=4), collapse = ', '),
                     ')', sep = '')  
  plot(x.axis.ticks, h1, type = "l", xlab = xl, 
       ylab = "Density", main = paste("Heterogeneity in", xl, sep = " "))
  mean.var.label <- paste("Mean: ", rate.mean, "    Var:", rate.var)
  mtext(mean.var.label, side = 3)
  lines(x.axis.ticks, h2, col = 'red')
  return(cbind(x.axis.ticks, h1, h2))
}

# pars is a matrix of two columns, r and alpha or s and beta, depending 
# on whether you want transaction rate or dropout rate heteroeneity.
# it has as many rows as there are estimates of each.
compareRH <- function(pars, rate = TRUE) {
  stopifnot(ncol(pars) == 2)
  if(colnames(pars)[1] == 'r' & colnames(pars)[2] == 'alpha') {
    xl  <- 'Transaction Rate'
    paran <- '(alpha is'
  } else if(colnames(pars)[1] == 's' & colnames(pars)[2] == 'beta') {
    xl <- 'Dropout Rate'
    paran <- '(beta is'
  } else {
    stop('Bad parameters.')
  }
  # if you wan to treat the second parameter as a scale
  # parameter, as opposed to rate, invert it here. that
  # way rest of the code (which assumes rate) still works.
  if(rate == FALSE) {
    pars[,2] = 1 / pars[,2]
    paran <- paste(paran, 'scale)', sep = ' ')
  } else {
    paran <- paste(paran, 'rate)', sep = ' ')
  }  
  mvl <- function(pars) {
    shape <- pars[,1]
    rate  <- pars[,2]
    r.mean <- shape/rate
    r.var  <- shape/rate^2
    lim    <- qgamma(0.99, shape = shape, rate = rate)
    return(cbind(pars, r.mean, r.var, lim))
  }

  myplots <- mvl(pars)
  x.axis.ticks <- seq(0, max(myplots[,'lim']), length.out = 100)
  h <- apply(pars, 1, function(x) dgamma(x.axis.ticks, shape = x[1], rate = x[2]))
  h <- cbind(x.axis.ticks, h)
  h <- h[rowSums(!apply(h, 2, is.infinite)) == ncol(h),] %>%
    as.data.frame()
  rate.mean <- paste('(', 
                     paste(round(myplots[,'r.mean'], digits=3), collapse = ', '),
                     ')', sep = '')
  rate.var <- paste('(', 
                    paste(round(myplots[,'r.var'], digits=3), collapse = ', '),
                    ')', sep = '')   
  lh <- tidyr::gather(h, key = q, value = pdens, -x.axis.ticks)
  ggplot(data = lh, aes(x = x.axis.ticks, y = pdens, group = q, colour = q)) + 
    geom_line() + 
    ylab('Density') + 
    xlab(xl) + 
    ggtitle(paste(paste("Heterogeneity in", xl, paran, sep = " "), 
                  paste("Mean: ", rate.mean, "\nVar:", rate.var), 
                  sep = "\n"))
}

# Store the goods here
vignettes <- list()
clocks <- list()

# Results using CRAN package
rm(list = ls()[grep('pnbd',ls())])
rm(list = ls()[grep('h2f1',ls())])
clocks$CRAN <- system.time(vignettes$CRAN <- replicateBTYD(cal.cbs = dplyr.cal.cbs))

# Results using BTYD2 package
rm(list = ls()[grep('pnbd',ls())])
rm(list = ls()[grep('h2f1',ls())])
clocks$BTYD2 <- system.time(vignettes$BTYD2 <- replicateBTYD(cal.cbs = dplyr.cal.cbs, mylib = 'BTYD2'))

# And with Re(hypergeo()) instead of h2f1()
clocks$hyper <- system.time(vignettes$hyper <- replicateBTYD(cal.cbs = dplyr.cal.cbs, mylib = 'BTYD2', hardie = FALSE))

# Results using Theo's script
rm(list = ls()[grep('pnbd',ls())])
rm(list = ls()[grep('h2f1',ls())])
clocks$Theo <- system.time(vignettes$Theo <- replicateBTYD(cal.cbs = dplyr.cal.cbs, theo = TRUE))

# Compare transactions and dropout rate heterogeneity between the three models
mypars <- vignettes %>% 
  map('Parameter estimates') %>% 
  bind_rows() %>% 
  t()
colnames(mypars) <- c('r', 'alpha', 's', 'beta')
mypars <- mypars %>% 
  cbind(elapsed_time = clocks %>% map_dbl('elapsed'))

# Compare transactions and dropout rate heterogeneity with ggplot()

# r and alpha
ralpha.rate <- compareRH(mypars[, c('r', 'alpha')], rate = TRUE)
ralpha.scale <- compareRH(mypars[, c('r', 'alpha')], rate = FALSE)

# s and beta
sbeta.rate <- compareRH(mypars[, c('s', 'beta')], rate = TRUE)
sbeta.scale <- compareRH(mypars[, c('s', 'beta')], rate = FALSE)

# s and beta no theo
sbeta.rate <- compareRH(mypars[c('CRAN', 'BTYD2'), c('s', 'beta')], rate = TRUE)
sbeta.scale <- compareRH(mypars[c('CRAN', 'BTYD2'), c('s', 'beta')], rate = FALSE)

# Compare transactions and dropout rate heterogeneity
# with plot() like the vignette does
comp.trans.hetero <- try(plotRH(vignettes$CRAN$`Parameter estimates`[1:2], 
                                vignettes$BTYD2$`Parameter estimates`[1:2]))
comp.drop.hetero <- try(plotRH(vignettes$CRAN$`Parameter estimates`[3:4], 
                               vignettes$BTYD2$`Parameter estimates`[3:4], 
                               trans = FALSE))
