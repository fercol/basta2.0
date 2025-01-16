\name{summary.basta}
\alias{summary.basta}
\alias{print.basta}
\alias{plot.basta}

\title{
Summarizing and plotting Bayesian Survival Trajectory Analysis (BaSTA) model outputs.
}
\description{
These functions are all generic methods for class \code{basta}.
}
\usage{
\method{summary}{basta}(object, \dots)
\method{print}{basta}(x, \dots)
\method{plot}{basta}(x, type = "traces", trace.name = "theta",
                       noCIs = FALSE, minSurv = NULL,
                       \dots)
}
\arguments{
  \item{object }{An object of class \code{basta}.}
  
  \item{x }{An object of class \code{basta}.}

  \item{type }{A character vector indicating the type of plot to be produced. Options are: \dQuote{\code{traces}} for the MCMC traces; \dQuote{\code{densities}} for posterior parameter densities; \dQuote{\code{demorates}} for the resulting survival and mortality curves; \dQuote{\code{gof}} for a comparison between the estimated parametric survival function and the life table \code{lx} variable; and \dQuote{\code{fancy}} for a combined densities and survival and mortality plots. }
  
  \item{trace.name }{Character string indicating the set of parameters or posteriors to be plotted. The options are: \dQuote{\code{theta}} to plot the survival model parameters; \dQuote{\code{gamma}} to plot the proportional hazards parameters (if it applies, else \code{plot.basta} returns an error); \dQuote{\code{pi}} to plot the recapture probabilities.}
  
  \item{noCIs }{Logical indicating whether the 95\% credible intervals should be included in the plot.}
  
  \item{minSurv }{Numerical value for the minimum survival level, acts as \code{xlim} using the survival as reference.}
  
  \item{\dots }{Additional arguments passed to functions \code{print}, \code{summary} and \code{plot} (see details).}
}
\details{
For objects of class \code{basta}, the \code{print} function returns three summary elements describing the model and its results, namely: \code{call}, \code{run}, \code{coefficients} and, if convergence was reached, the \code{DIC} values for model fit. \code{call} describes the basic model used (i.e. exponential, Gompertz, Weibull or logistic), the shape chosen, \dQuote{\code{simple}}, \dQuote{\code{Makeham}} or \dQuote{\code{bathtub}}, the covariate structure chosen, \dQuote{\code{fused}}, \dQuote{\code{prop.haz}} or \dQuote{\code{all.in.mort}} and which covariates (if any) were categorical and which continuous. Argument \code{digits} can be used for number formatting (see \code{summary}() or \code{signif}() for details).  

The summary element \code{coefficients} prints out the estimated coefficients for all parameters in the model, as well as their standard errors and 95\% upper and lower credible intervals. It also includes a measure of serial autocorrelation for each parameter calculated from the thinned parameter chains, an update rate per parameter, and the potential scale reduction factor for each parameter as a measure of convergence (Gelman \emph{et al.} 2004). 

Function \code{summary} includes all the previous elements, as well as a summary description of the priors and jump standard deviations for all survival parameters, a calibration of the Kullback-Leibler discrepancy as a measure of parameter similarities for those parameters associated to categorical covariates (McCulloch 1989), and a measure of model fit based on the deviance information criterion (DIC) (Spiegelhalter \emph{et al.} 2002). 

Function \code{plot} takes objects of class \code{basta} to create trace plots of the estimated parameters (default) or, if the argument for \code{type} is set to \dQuote{\code{demorates}}, it plots estimated survival probabilities and mortality rates with their 95\% predictive intervals. Argument \code{xlim} can be used to define a range of ages to visualize survival and mortality trends. Also, if logical argument \code{noCI} is set to \code{TRUE}, credible intervals around survival and mortality are not plotted, leaving only the mean trends. This can be handy when several categorical covariates have been evaluated and the plots become too crowded.  

Other arguments for \code{plot} include \code{names.legend} to indicate alternative names for the legend of vital rates plots. Also, argument \code{xlim} can be changed to display only a subest of the support. When argument \code{noCI} is \code{TRUE}, then the credible intervals around the vital rates are not displayed. In addition, argument \code{lwd} can be used to change the width of the survival and mortality mean lines.

For consistency with previous versions, arguments \code{plot.type} and \code{densities} can still be used instead of \code{type}. 
}

\value{ 
Function \code{summary}() outputs the folowing values:

\item{coefficients }{Matrix with estimated coefficients (i.e. mean values per parameter on the thinned sequences after burnin), which includes standard errors, upper and lower 95\% credible intervals, update rates per parameter (commonly the same for all survival and proportional hazards parameters), serial autocorrelation on the thinned sequences and the potential scale reduction factor for convergence (see \code{Convergence} value below). }
  
\item{DIC }{Numeric vector with basic deviance information criterion (DIC) calculations to be used for model selection (Spiegelhalter \emph{et al.} 2002).}
	
\item{KullbackLeibler }{List with Kullback-Leibler discrepancy matrices between pair of parameters for categorical covariates (McCulloch 1989, Burnham and Anderson 2001) and McCulloch's (1989) calibration measure. If only one simulation was ran or if no convergence was reached, then the returned value is \dQuote{\code{Not calculated}}.}

\item{convergence }{A matrix with convergence coefficients based on potential scale reduction as described by Gelman \emph{et al.} (2004). If only one simulation was ran, then the returned value is \dQuote{\code{Not calculated}}.}
    
\item{modelSpecs }{Model specifications inidicating the \code{model}, the \code{shape} and the covariate structure that were specified by the user.}

\item{settings }{A vector indicating the number of iterations for each MCMC, the burn in sequence, the thinning interval, and the number of simulations that were run.}
	
}

\references{
Gelman, A., Carlin, J.B., Stern, H.S. and Rubin, D.B. (2004) \emph{Bayesian data analysis}. 2nd edn. Chapman & Hall/CRC, Boca Raton, Florida, USA.

McCulloch, R.E. (1989) Local model influence. \emph{Journal of the American Statistical Association}, 84, 473-478.

Spiegelhalter, D.J., Best, N.G., Carlin, B.P. and Van Der Linde, A. (2002) Bayesian measures of model complexity and fit. \emph{Journal of the Royal Statistical Society: Series B} 64, 583-639.

See also:

Colchero, F. and J.S. Clark (2012) Bayesian inference on age-specific survival from capture-recapture data for censored and truncated data. \emph{Journal of Animal Ecology}. 81(1):139-149.

Colchero, F., O.R. Jones and M. Rebke. (2012) BaSTA: an R package for Bayesian estimation of age-specific survival from incomplete mark-recapture/recovery data with covariates. \emph{Method in Ecology and Evolution}. DOI: 10.1111/j.2041-210X.2012.00186.x
}
\author{
Fernando Colchero \email{fernando_colchero@eva.mpg.de}, Owen R. Jones \email{jones@biology.sdu.dk}
}

\seealso{
\code{\link{basta}}
}

\examples{
## Load BaSTA output:
data("bastaCMRout", package = "BaSTA")

## Print summary output:
summary(bastaCMRout)

## Plot traces for mortality parameters (theta):
plot(bastaCMRout)

## Plot traces for proportional hazards parameters (gamma):
plot(bastaCMRout, trace.name = "gamma")

## Plot traces for recapture probability(ies) (pi):
plot(bastaCMRout, trace.name = "pi")

## Plot predicted mortality and survival:
plot(bastaCMRout, type = "demorates")

## Change the color for each covariate on 
## the predicted vital rates:
plot(bastaCMRout, type = "demorates", 
     col = c("dark green", "dark blue"))

## Change the color and the legend text:
plot(bastaCMRout, type = "demorates", 
     col = c("dark green", "dark blue"),
     names.legend = c("Females", "Males"))

## Plot predicted mortality and survival 
## between 2 and 8 years of age:
plot(bastaCMRout, type = "demorates", xlim = c(2, 8))

## Plot predicted mortality and survival 
## between 2 and 8 years of age without
## credible intervals:
plot(bastaCMRout, type = "demorates", xlim = c(2, 8), 
     noCI = TRUE)

## Plot parameter densities and predicted vital  
## rates in the same plot (i.e. fancy):
plot(bastaCMRout, type = "fancy")

## Change colors and legend names for the 
## "fancy" plot:
plot(bastaCMRout, type = "fancy", col = c("dark green", "dark blue"),
     names.legend = c("Females", "Males"))
}

\keyword{methods}
