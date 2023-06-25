\name{semi_rk_est}
\alias{semi_rk_est}
\title{Calculate the Subsample Estimator}
\usage{
semi_rk_est(x, y, delta, pi, n, 
            control = list(init = rep(0, ncol(x)), 
            tol = 1e-5, maxit = 200))
}
\arguments{
\item{x}{the input covariate matrix.}

\item{y}{the input response vector.}

\item{delta}{the input censoring indicator.}

\item{pi}{the subsampling probabilities.}

\item{n}{the size of the full sample.}

\item{control}{the list of controling variables, including the initial value
               ,tolerance and the maximum number of iterations.}
}

\value{
\item{coefficient}{the subample estimator.}

\item{converge}{the converging message: 0 for successful converging;
                1 for failing to converge;
                2 for the iteration limit `maxit` has been reached}

\item{ite}{the number of iterations if it successfully converges.}
}

\description{
This function calculate the subsampling estiamtor for the subsample
covariate matrix x, response vector y, censoring indicator delta and the 
subsampling probabilities. It returns a list with three elements: 
coefficient, the subample estimator; 
coverge, the converging message; ite, the number of itterations used.
}
\examples{
library(semiaftosmac)
# Load the full data
data(lymp)
x <- lymp$covariates
y <- log(lymp$response)
delta <- lymp$delta
n <- nrow(x)

# Sampling with replacement based on uniform subsampling method
# Subsample size 3000
set.seed(101)
ind <- sample(n, 3000, replace = TRUE)

# Get the subsample estimator
semi_ls_est(x[ind, ], y[ind], delta[ind], rep(1/n, 3000), n)
}
\keyword{semi_ls_est}
