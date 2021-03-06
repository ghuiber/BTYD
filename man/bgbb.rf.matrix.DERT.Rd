% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bgbb.R
\name{bgbb.rf.matrix.DERT}
\alias{bgbb.rf.matrix.DERT}
\title{BG/BB Discounted Expected Residual Transactions using a recency-frequency matrix}
\usage{
bgbb.rf.matrix.DERT(params, rf.matrix, d)
}
\arguments{
\item{params}{BG/BB parameters - a vector with alpha, beta, gamma, and
delta, in that order. Alpha and beta are unobserved parameters for the
beta-Bernoulli transaction process. Gamma and delta are unobserved
parameters for the beta-geometric dropout process.}

\item{rf.matrix}{recency-frequency matrix. It must contain columns for
frequency ("x"), recency ("t.x"), number of transaction opportunities in
the calibration period ("n.cal"), and the number of customers with this
combination of recency, frequency and transaction opportunities in the
calibration period ("custs"). Note that recency must be the time between
the start of the calibration period and the customer's last transaction,
not the time between the customer's last transaction and the end of the
calibration period.}

\item{d}{discount rate.}
}
\value{
The present value of the expected future transaction stream for a particular customer.
}
\description{
Computes the number of discounted expected residual transactions by a
customer, conditional on their behavior in the calibration period.
}
\examples{
data(donationsSummary)

rf.matrix <- donationsSummary$rf.matrix
# donationsSummary$rf.matrix already has appropriate column names

# starting-point parameters
startingparams <- c(1, 1, 0.5, 3)
# estimated parameters
est.params <- bgbb.EstimateParameters(rf.matrix, startingparams)

# compute DERT for a customer from every row in rf.matrix,
# discounted at 10\%.
bgbb.rf.matrix.DERT(est.params, rf.matrix, d = 0.1)
}
\references{
Fader, Peter S., Bruce G.S. Hardie, and Jen Shang.
"Customer-Base Analysis in a Discrete-Time Noncontractual Setting."
\emph{Marketing Science} 29(6), pp. 1086-1108. 2010. INFORMS.
\href{http://www.brucehardie.com/papers/020/}{Web.}
See equation 14.
}
\seealso{
\code{\link{bgbb.DERT}}
}
