\name{semi_rk_fit}
\alias{semi_rk_fit}
\title{The Two-step Procedure}
\usage{
semi_rk_fit(x, y, delta, r0, r, ssp_type, 
            se = TRUE, alpha = 0.2) 
}
\arguments{
\item{x}{the input covariate matrix.}

\item{y}{the input response vector.}

\item{delta}{the input censoring indicator.}

\item{r0}{the input pilot subsample size.}

\item{r}{the input second-step subsample size.}

\item{ssp_type}{the type of optimal subsampling probabilities (SSPs). `optA` for the
                `optA` SSPs; `optL` for the `optL` SSPs; `uniform` for the
                uniform SSPs.}
              
\item{se}{the indicator of calculating standard errors. The default is `TRUE`.}
           
\item{alpha}{propotion of uniform SSPs mixed with the optimal SSPs. The default is 0.2.}
}

\value{
\item{coefficient}{the subample estimator.}

\item{converge}{the converging message: 
                0 for successful converging;
                1 for failing to converge;
                2 for the iteration limit `maxit` has been reached.
                }

\item{ite}{the number of iterations if it successfully converges.}

\item{std}{the standard error estimator}
}

\description{
This function apply the two-step method for the full sample
covariate matrix x, response vector y and censoring indicator delta with user
specified pilot subsample size and second-step subsample size.
For more arguments, please read the details section.
It returns a list with four elements: coefficient, the subample estimator; 
coverge, the converging message; ite, the number of iterations used; std, the
standard error estimator.
}
\examples{
library(semirkosmac)
# Load the full data
data(lymp)
x <- lymp$covariates
y <- log(lymp$response)
delta <- lymp$delta
n <- nrow(x)

# Get the subsample estimator using A-optimal SSPs when \code{r0 = 500} and \code{r = 2000}.
set.seed(101)
semi_rk_fit(x, y, delta, r0 = 500, r = 2000, 
            ssp_type = "optA", se = "TRUE", 
            alpha = 0.2)
}
\keyword{semi_rk_fit}
