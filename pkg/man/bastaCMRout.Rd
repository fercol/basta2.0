\name{bastaCMRout}
\alias{bastaCMRout}
\docType{data}
\title{
Output from a Bayesian Survival Trajectory Analysis (BaSTA) analysis on a simulated capture-mark-recapture (CMR) dataset.
}
\description{
This dataset is the output of a BaSTA analysis on the simulated dataset \code{\link{bastaCMRdat}}. The analysis consisted of four independent simulations run in parallel. Each simulation was run for 20,000 iterations. The model chosen was Gompertz (\dQuote{\code{GO}}) with the \code{shape} argument set to \dQuote{\code{simple}} and \code{covarStruct} set to \dQuote{\code{fused}}.   
}

\usage{data("bastaCMRout")}
\format{
  A BaSTA output list (see details in \code{\link{basta}}). 
}
\examples{
## Load BaSTA output:
data(bastaCMRout)

## Plot traces for survival parameters:
plot(bastaCMRout)

## Plot posterior densities of survival parameters:
plot(bastaCMRout, densities = TRUE)

## Plot traces for proportional hazards parameter:
plot(bastaCMRout, trace.name = "gamma")

## Plot survival and mortality curves:
plot(bastaCMRout, plot.type = "demorates")
}
\keyword{model output}
