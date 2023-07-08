\name{semi_rk_ssp}
\alias{semi_rk_ssp}
\title{Calculate the Optimal Subsample Probabilities}
\usage{
semi_rk_ssp(x, y, delta, r0, ssp_type, alpha)
}
\arguments{
\item{x}{the input covariate matrix.}

\item{y}{the input response vector.}

\item{delta}{the input censoring indicator.}

\item{r0}{the input pilot subsample size.}

\item{ssp_type}{the type of optimal subsampling probabilities (SSPs). `optA` for the
                `optA` SSPs; `optL` for the `optL` SSPs; `uniform` for the
                uniform SSPs.}
         
\item{alpha}{propotion of uniform SSPs mixed with the optimal SSPs.}
}

\value{
\item{ssp}{the subsampling probabilities.}

\item{converge}{the converging message of the pilot estimator: 
                0 for successful converging;
                1 for failing to converge;
                2 for the iteration limit `maxit` has been reached}

\item{index.pilot}{the index of pilot subsample in the full sample.}
}

\description{
This function calculate the optimal subsampling probabilities for the full sample.
The inputs are full sample covariate matrix \code{x}, response vector \code{y} and censoring
indicator delta. Moreover, the pilot subsample size \code{r0} needs to be specified. For
more arguments, please read the details section.
It returns a list with three elements: ssp, the subsampling probabilities;
index.pilot, the index of pilot subsample in the full sample; 
coverge, the converging message of the pilot estimator.
}
\examples{
library(semirkosmac)
# Load the full data
data(lymp)
x <- lymp$covariates
y <- log(lymp$response)
delta <- lymp$delta
n <- nrow(x)

# Get the optimal SSPs based on A-optimality
set.seed(101)
ssps <- semi_rk_ssp(x, y, delta, r0 = 500, 
                    ssp_type = "optA", alpha = 0.2)
}
\keyword{semi_rk_ssp}
