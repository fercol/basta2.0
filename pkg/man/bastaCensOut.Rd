\name{bastaCensOut}
\alias{bastaCensOut}
\docType{data}
\title{
Output from a Bayesian Survival Trajectory Analysis (BaSTA) analysis on a simulated census dataset.
}
\description{
This dataset is the output of a BaSTA analysis on the simulated dataset \code{\link{bastaCensDat}}. The analysis consisted of four independent simulations run in parallel. The model chosen was Gompertz (\dQuote{\code{GO}}) with the \code{shape} argument set to \dQuote{\code{simple}} and \code{covarStruct} set to \dQuote{\code{fused}}.   
}
\usage{data("bastaCensOut")}
\format{
  A BaSTA output list (see details in \code{\link{basta}}). 
}
\examples{
## Load BaSTA output:
data("bastaCensOut", package = "BaSTA")

## Plot traces for survival parameters:
plot(bastaCensOut)

## Plot posterior densities of survival parameters:
plot(bastaCensOut, densities = TRUE)

## Plot survival and mortality curves:
plot(bastaCensOut, plot.type = "demorates")
}
\keyword{model output}
