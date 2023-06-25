\name{semi_rk_fit}
\alias{semi_rk_fit}
\title{The Two-step Procedure}
\usage{
semi_rk_fit(x, y, delta, r0, r, ssp_type, method, 
            se = TRUE, b = 20, itr = 5, alpha = 0.2) 
}
\arguments{
\item{x}{the input covariate matrix.}

\item{y}{the input response vector.}

\item{delta}{the input censoring indicator.}

\item{pi}{the subsampling probabilities.}

\item{r0}{the input pilot subsample size.}

\item{r}{the input second-step subsample size.}

\item{ssp_type}{the type of optimal subsampling probabilities (SSPs). `optA` for the
                `optA` SSPs; `optL` for the `optL` SSPs; `uniform` for the
                uniform SSPs.}
                
\item{method}{the method to calculate the subsample estimator and standard errors.
              `one` for using one second-step subsample to calculate the resultant
              estimator and the variance is calculated based on a sandwich form.
              `multiple` for using multiple second-step subsample to calculate the
              resultant estimator and the standard is calculated by the sample
              standard error.}
              
\item{se}{the indicator of calculating standard errors. The default is `TRUE`.}

\item{b}{the number of noises generated in the resampling method to estimate the
         slope matrix. This argument is necessary only when the `ssp_type` argument
         being `optA` or the `method` argument being `one`.}
         
\item{itr}{the number of subsample estimators are calculated when the `method`
           argument being `multiple`.}
           
\item{alpha}{propotion of uniform SSPs mixed with the optimal SSPs.}
}

\value{
\item{coefficient}{the subample estimator.}

\item{converge}{the converging message: 
                0 for successful converging;
                1 for failing to converge;
                2 for the iteration limit `maxit` has been reached.
                If the `method` argument is `multiple` and at least one 
                subsampling procedure does not converge, the output is a matrix
                with the first row representing the number of `1` and the second
                row representing the number of `2` among all non-converging subsampling
                procedures.
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
library(semiaftosmac)
# Load the full data
data(lymp)
x <- lymp$covariates
y <- log(lymp$response)
delta <- lymp$delta
n <- nrow(x)

# Get the subsample estimator using A-optimal SSPs when `r0 = 2000` and `r = 4000`.
set.seed(101)
semi_ls_fit(x, y, delta, r0 = 3000, r = 5000, ssp_type = "optA",
            method = "one", se = "TRUE", b = 100, alpha = 0.2)
}
\keyword{semi_ls_fit}
