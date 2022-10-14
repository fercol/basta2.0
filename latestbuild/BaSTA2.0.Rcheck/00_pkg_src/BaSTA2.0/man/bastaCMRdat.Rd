\name{bastaCMRdat}
\alias{bastaCMRdat}
\docType{data}
\title{
Example of capture-mark-recapture data for BaSTA analysis.
}
\description{
Simulated capture-mark-recapture (\code{CMR}) data for Bayesian Survival Trajectory Analysis (BaSTA).
}
\usage{data("bastaCMRdat")}
\format{
  A data frame with 500 observations on the following 25 variables (see \code{details}).
}
\details{
This dataset was created by stochastically simulating a hypothetical population with different mortality patterns between males and females and with proportional decreases in mortality as a function of a hypothetical continuous covariate (e.g. birth weight). The population was simulated for 100 years, at each one of which 100 individuals were born. The number of females per generation was randomly drawn from a binomial distribution using function \code{rbinom} with probability of 0.5 (i.e. 1:1 sex ratio). The individual continuous covariate was randomly drawn from a random normal distribution (with function \code{rnorm}) with mean parameter equal to 0 (e.g. anomaly of weights) and standard deviation equal to 1.  The time of death for each individual was inversed sampled from a Gompertz CDF of ages at death. The Gompertz parameters for females were: \eqn{b_0 = -4} and \eqn{b_1 = 0.15}; and for males at \eqn{b_0 = -3} and \eqn{b_1 = 0.15}. The gamma parameter for the continuous covariate was \eqn{\gamma = 0.2}.

The study was assumed to start at year 51 and to finish at year 70. Recapture probability was set to 0.6 and thus each observation per individual was randomly drawn from a Bernoulli trial with parameter \eqn{p = 0.6}. Captures at birth and recoveries were randomly drawn from a Bernoulli trial with parameters \eqn{p = 0.5} and \eqn{p = 0.2}, respectively. 

Therefore, the resulting dataset includes individuals where the data are left-truncated and/or right-censored. This is typical of capture mark recovery datasets.

Here is the description of each column:

\code{id}: ID for each individual.

\code{birth}: Integer times of birth or \code{NA} when unknown.

\code{death}: Integer times of death or \code{NA} when unknown.

\code{Y.51 - Y.70}: Columns of the recapture matrix

\code{Sex}: Character string indicating the sex for the individual record (i.e., \dQuote{\code{Female}} and \dQuote{\code{Male}})

\code{Weight}: Simumated body weights at birth.

}
\examples{
## Load data:
data("bastaCMRdat", package = "BaSTA2.0")

## Check data consistency:
checkedData  <- DataCheck(bastaCMRdat, dataType = "CMR", studyStart = 51, 
                          studyEnd = 70)
}
\keyword{datasets}
