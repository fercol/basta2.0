\name{basta}
\alias{basta}
\alias{basta.default}

\title{
 Parametric Bayesian estimation of age-specific survival for left-truncated and 
right-censored capture-mark-recapture or census data.
}

\description{
 This function performs multiple Markov Chain Monte Carlo (MCMC) simulations for the Bayesian estimation of age-specific mortality and survival when a large proportion of records have unknown times of birth and/or death. Survival parameters and unknown (i.e. latent) birth (and death, for CMR data) times are estimated, allowing the user to test a range of mortality patterns, and to test the effect of continuous and/or discrete covariates following Colchero and Clark's (2012) general approach.
}

\usage{
basta(object, \dots)

\method{basta}{default} (object, dataType = "CMR", model = "GO", 
  shape = "simple",  studyStart = NULL, studyEnd = NULL, 
  minAge = 0, covarsStruct = "fused",  formulaMort = NULL, 
  formulaRecap = NULL, recaptTrans = studyStart,  niter = 11000, 
  burnin = 1001, thinning = 20, nsim = 1, parallel = FALSE,  
  ncpus = 2, updateJumps = TRUE, negSenescence = FALSE,  \dots)
}

\arguments{
	\item{ object }{A \code{data.frame} to be used as an input data file for BaSTA. Note: BaSTA can take two types of datasets, namely capture-mark-recapture (CMR) or census data. See \code{details} for a description of the data structures. }

	\item{dataType }{A \code{character string} indicating if the data are for capture-mark-recapture (CMR) or census. Options are \dQuote{\code{CMR}} (default) or \dQuote{\code{census}}.}

	\item{model }{The underlying mortality model to be used. \code{"EX"} = exponential,\code{"GO"} = Gompertz, \code{"WE"} = Weibull and \code{"LO"} = logistic (see \code{details}).}

	\item{shape }{The overall shape of the model. Values are: \code{simple} = no extra parameters added; \code{Makeham} = a constant parameter is added to the mortality; and \code{bathtub} = a Gompertz declining mortality for early ages and a constant parameter are added to the mortality model (see \code{details}).}

	\item{ studyStart }{Only required for \code{dataType =} \dQuote{\code{CMR}}, an integer indicating the first year of the study.}

	\item{studyEnd }{Only required for \code{dataType =} \dQuote{\code{CMR}}, an integer indicating the last year of the study.}

	\item{minAge }{Age at which the analysis should start (see \code{details}).}

	\item{covarsStruct }{Character string that indicates how covariates should be evaluated. The options are: \dQuote{\code{fused}}, which defines all categorical variables as covariates for each mortality parameter and all continuous covariates under a proportional hazards structure; \dQuote{\code{prop.haz}}, which puts all covariates under a proportional hazards structure; and \dQuote{\code{all.in.mort}} puts all covariates as a multilevel function of the mortality parameters (see \code{details}).}

	\item{ formulaMort }{An object of class \code{\link{formula}} specifying the covariates to be included on the mortality part of the model. Note that the syntax should not include a dependent variable. If \code{formulaMort = NULL}, then no covariates are included in the analysis. See \code{details}.}

	\item{formulaRecap }{\emph{Not yet implemented}, an object of class \code{\link{formula}} specifying the covariates to be included on the recapture part of the model. Note that the syntax should not include a dependent variable. If \code{formulaMort = NULL}, then no covariates are included in the analysis.}

	\item{recaptTrans }{A vector (of maximum length equal to the duration of the study) defining the recapture probability transition times (RPTP). These are points (years) where the recapture probability is thought to change.  The default setting is for the recapture probability to be constant throughout the study, so the \code{recaptTrans} is simply defined as a single element vector of the first year of the observation period (e.g. c(1985)). If recapture probabilities were known to change at year say, 1990, the RPTP should be defined as c(1985, 1990).}

\item{niter }{The total number of MCMC steps.}
  
\item{burnin }{The number of iterations for the burn in (see \code{details}).}
  
\item{thinning }{The number of skipped MCMC steps to minimize serial autocorrelation (see \code{details}).}

\item{nsim }{A numerical value for the number of simulations to be run.}
  
\item{parallel }{A logical argument indicating whether the multiple simulations should be run in parallel or not. If \code{TRUE}, package \pkg{snowfall} is called and multiple simulations are run in parallel. If \pkg{snowfall} is not installed, the model is ran in series.}
  
\item{ncpus }{a numerical value that indicates the number of cpus to be used if \code{parallel} is \code{TRUE} and package \pkg{snowfall} is installed. The default is 2 cpus. If package pkg{snowfall} is not installed, the simulations are run in series.}

	\item{updateJumps }{A logical argument indicating wheter to update jump standard deviations (adaptive independent Metropolis) until an update rate of 0.25 is achieved (see \code{details}). }

	\item{negSenescence }{Logical indicating if negative senescence should be allowed, only applicable for \code{model = \dQuote{GO}} and \code{shape = \dQuote{Makeham}}.}

	\item{\dots }{Additional arguments to be passed to function \code{basta} (see \code{details})}

}

\details{
\code{1) DATA TYPES:}

The input \code{object} required by BaSTA needs to be constructed differently whether the data are of capture-mark-recapture or census. 

\code{a) capture-mark-recapture (CMR):}

If \code{dataType = \dQuote{CMR}}, then the data frame requires the following structure. The first column is a vector of individual unique IDs, the second and third columns are birth and death years respectively. Columns \eqn{4, \dots, T+3} represent the observation window (i.e., recapture matrix) of \eqn{T} years. This is followed (optionally) by columns for categorical and continuous covariates (see \code{\link{bastaCMRdat}} for a CMR dataset example).

\code{b) Census:}

If \code{dataType = \dQuote{census}}, then the input data \code{object} requires at least five dates columns, namely \dQuote{Birth.Date}, \dQuote{Min.Birth.Date}, \dQuote{Max.Birth.Date}, \dQuote{Entry.Date}, and \dQuote{Depart.Date}. All dates need to be format as \dQuote{\%Y-\%m-\%d}. In addition, a \dQuote{Depart.Type} column is required with two types of departures \dQuote{C} for Censored and \dQuote{D} for dead (see \code{\link{bastaCensDat}} for a census dataset example).


\code{2) INFERENCE ON AGE-SPECIFIC MORTALITY:}

\code{basta} uses parametric mortality functions to estimate age-specific mortality (survival) from capture-recapture/recovery data. The mortality function describes how the risk of mortality changes with age, and is defined as \eqn{\mu(x | \theta)}, where \eqn{x} corresponds to age and \eqn{\theta} is a vector of parameters to be estimated. 

The \code{model} argument allows the user to choose between four basic mortality functions, namely (a) Exponential (\dQuote{\code{EX}}; Cox and Oakes 1974), with constant mortality with age specified as

\deqn{\mu_b(x | b) = b,}

where \eqn{b > 0}, (b) the Gompertz mortality function (\dQuote{\code{GO}}; Gompertz 1925, Pletcher 1999), calculated as 

\deqn{\mu_b(x | b) = exp(b_0 + b_1 x),}

where \eqn{-\infty < b_0, b_1 < \infty}, (c) the Weibull mortality model (\dQuote{\code{WE}}; Pinder III \emph{et al.} 1978) calculated as

\deqn{\mu_b(x | b) = b_0 b_1^{b_0} x^{b_0 -1},}

where \eqn{b_0, b_1 > 0}, and (d) the logistic mortality model (\dQuote{\code{LO}}; Pletcher 1999), calculated as 

\deqn{\mu_b(x | b) = exp(b_0 + b_1 x) / (1 + b_2 exp(b_0)/b_1 (exp(b_1 x)-1)),}
where \eqn{b_0, b_1, b_2 > 0}.

The \code{shape} argument allows the user to extend these models in order to explore more complex mortality shapes. The default value is \dQuote{\code{simple}} which leaves the model as defined above. With value \dQuote{\code{Makeham}}, a constant is added to the mortality, making the model equal to \eqn{\mu_0(x | \theta)= \mu_b(x | b) + c}, where \eqn{\theta = [c, b]}. With value \dQuote{\code{bathtub}}, concave shapes in mortality can be explored. This is achieved by adding a declining Gompertz term and a constant parameter to the basic mortality model, namely

\deqn{\mu_0(x | \theta) = exp(a_0 - a_1 x) + c + \mu_b(x | b)},

where \eqn{-\infty < a_0 < \infty}, \eqn{a_1 \leq 0} and \eqn{c \leq 0}.

\code{3) COVARIATES:}

Covariates are selected by means of the argument \code{formulaMort}, which requires an object of class formula, just as with other statistical inference functions such as \code{lm} or \code{glm}.

When covariates are included in the dataset, the \code{basta} function provides three different ways in which these can be evaluated by using argument \code{covarsStruct}: 

1. \dQuote{\code{fused}} will make the mortality parameters linear functions of all categorical covariates (analogous to a generalised linear model (GLM) structure) and will put all continuous covariates under a proportional hazards structure. Thus, for a simple exponential model with constant mortality of the form \eqn{\mu_0(x | b) = b}, the parameter is equal to \eqn{b = b_0 + b_1 z_1 + \dots, b_k z_k}, where \eqn{[b_0, \dots, b_k]} are paramters that link the mortality parameter \eqn{b} with the categorical covariates \eqn{[z_1,\dots,z_k]}.

2. \dQuote{\code{prop.haz}} will put all covariates under a proportional hazards structure irrespective of the type of variable. In this case, the mortality model is be further extended by including a proportional hazards structure, of the form 

\deqn{\mu(x | \theta, \Gamma, Z_a, Z_c) = \mu_0(x | \theta, Z_a) exp(\Gamma Z_c)},

where \eqn{\mu_0(x | \theta, Z_a)} represents the mortality section as defined above, while the second term \eqn{exp(\Gamma Z_c)} corresponds to the proportional hazards function.  \eqn{Z_a} and \eqn{Z_c} are covariate (design) matrices for categorical and continuous covariates, respectively.

3. \dQuote{\code{all.in.mort}} will put all covariates as linear functions of the survival parameters as explained above. Since most models require the lower bounds for the mortality parameters to be equal to 0, the only model that can be used for this test is Gompertz with \code{shape} set to \dQuote{\code{simple}}. In case these arguments are specified deferently, a warning message is printed noting that \code{model} will be forced to be \dQuote{\code{GO}} and \code{shape} will be set to \dQuote{\code{simple}}.

\code{4) MCMC SETTINGS:}

The \code{burnin} argument represents the number of steps at the begining of the MCMC run that is be discarded. This sequence commonly corresponds to the non-converged section of the MCMC sequence. Convergence and model selection measures are calculated from the remaining thinned parameter chains if multiple simulations are run, and all if all of them run to completion.

The \code{thinning} argument specifies the number of steps to be skipped in order to reduce serial autocorrelation. The thinned sequence, which only includes steps after burn in, is then used to calculate convergence statistics and model for selection.

The \code{updateJumps} argument specifies wether to run a simulation to find appropriate jump standard deviations for theta and gamma parameters. If argument \dQuote{\code{nsim}} is set to 1, then the simulation runs with the update jumps routine active. If \dQuote{\code{nsim}} is larger than 1, then an initial simulation is ran to find apropriate jumps before the main analysis is ran. 

\code{5) ADDITIONAL ARGUMENTS:}

Additional arguments for priors, jumps and start values can be passed on the \dots section. For instance, argument \code{thetaStart} can be specified as a vector defining the initial values for each parameter in the survival model. If this argument is not specified, a set of random parameters is generated for each simulation. Similarly, argument \code{gammaStart} can be specified for all parameters in the proportional hazards section of the model. Jump standard deviations (i.e. the standard error in the Metropolis step) can be specified with arguments \code{thetaJumps} and \code{gammaJumps}. As with \code{thetaStart}, default values are assigned if these arguments are not specified. 

To specify priors, arguments \code{thetaPriorMean}, \code{thetaPriorSd}, \code{gammaPriorMean} and \code{gammaPriorSd} can be used for prior means and standard errors for each survival and proportional hazards parameters. If not specified, default values are assigned.

The number of parameters in \code{thetaStart}, \code{thetaJumps}, \code{thetaPriorMean} and \code{thetaPriorSd} should be a vector or matrix for the parameters in the mortality function. The number of parameters will depend on the model chosen with \code{model} (see above). If the number of parameters specified does not match the number of parameters inherent to the model and shape selected, the function returns an error. 

As described above, the number of parameters for \code{gammaStart}, \code{gammaJumps}, \code{gammaPriorMean} and \code{gammaPriorSd} arguments (i.e. section b), namely the proportional hazards section, will be a function of the number of continuous covariates if argument \code{covarsStruct} is \dQuote{\code{fused}}, or to the total number of covariates when \code{covarsStruct} is \dQuote{\code{prop.haz}}.

%Another additional argument is \code{returnAges}, which outputs a value \code{estAges}, which is a matrix with all estimated ages after the burnin and thinning. This matrix also includes ages that were known, thus each row corresponds to the individual in the same row as the original dataset.

}

\value{
	\item{params }{If requested, a matrix with the thinned, converged parameter traces of all runs. This matrix is used to calculate quantiles for parameters, survival probability and mortality (see below).}

\item{theta }{If requested, a matrix with only the parameters of the mortality function after convergence and thinning.} 

\item{coefficients }{A matrix with estimated coefficients (i.e. mean values per parameter on the thinned sequences after burnin), which includes standard errors, upper and lower 95\% credible intervals, update rates per parameter, serial autocorrelation on the thinned sequences and the potential scale reduction factor for convergence (see \code{Convergence} value below). }

\item{names }{Names of all parameters}

\item{DIC }{Basic deviance information criterion (DIC) calculations to be used for model selection (Spiegelhalter \emph{et al.} 2002, Celeux \emph{et al.} 2006). Small differences between values should only be used a reference (see comments in Spiegelhalter \emph{et al.} 2002). If all or some of the simulations failed, then the returned value is \dQuote{\code{Not calculated}}.}

\item{KullbackLeibler }{If called by \code{summary}, list with Kullback-Leibler discrepancy matrices between pair of parameters for categorical covariates (McCulloch 1989, Burnham and Anderson 2001) and McCulloch's (1989) calibration measure. If only one simulation was run or if no convergence is reached, then the returned value is \dQuote{\code{Not calculated}}.}

\item{PS }{If requested, a list with summary statistics of the PDF of ages at death, including the life expectancy, lifespan inequality, lifespan equality, and Gini index. These are separated by categorical covariate and, if continuous covariates are provided, they are evaluated at the average value of each continous covariate. The list object provides a table with the mean and lower and upper 95\% credible intervals and vectors of the converged and thinned values for each variable.}

\item{mort }{If requested or called by functions \code{plot} or \code{\link{plot.basta}} median and 95\% predictive intervals for the estimated mortality rates separated by categorical covariate and calculated at the mean for each continuous covariate, if provided.}

\item{surv }{If requested or called by functions \code{plot} or \code{\link{plot.bastaZIMS}} median and 95\% predictive intervals for the estimated survival probability, separated by categorical covariate and calculated at the mean for each continuous covariate, if provided.}

\item{dens }{If requested, median and 95\% predictive intervals for the estimated probability density function of ages at death, separated by categorical covariate and calculated at the mean for each continuous covariate, if provided.}

\item{x }{If requested, a vector of the ages used to calculate \code{mort}, \code{surv}, and \code{dens}.}

\item{cuts }{An index vector per categorical covariate of the ages where the survival is larger than 0.05, used for display purposes when producing the plots with function \code{plot}.}

\item{convergence }{If requested, a matrix with convergence coefficients based on potential scale reduction as described by Gelman \emph{et al.} (2004). If only one simulation was ran, then the returned value is \dQuote{\code{Not calculated}}.}

\item{convmessage }{Only used with functions \code{summary} and \code{print} to indicate whether the parameter traces converged appropriately.}

\item{runs }{A list object with the outputs of each individual MCMC run. Used with function \code{plot}.}

\item{fullpar }{A list object with the input parameter information for the model, including starting values, priors, initial jumps, lower bound, among other. Used with functions \code{summary} and \code{print}.}

\item{simthe }{A list object with information on the basic mortality model. Used with function \code{plot}.}

\item{jumps }{A list object with the final jump standard deviations for each parameter.}

\item{covs }{A list object with general information on the type of covariates, i.e., \code{cat} and \code{cont}, and the type of \code{covarStruct} used.}

\item{settings }{If called by \code{summary}, this is a vector indicating the number of iterations for each MCMC, the burn in sequence, the thinning interval, and the number of simulations that were run.}

\item{modelSpecs }{Model specifications inidicating the \code{model}, the \code{shape}, the covariate structure and the minimum age that were specified by the user.}
	
\item{lifeTable }{A period life table calculated from the estimated times of birth (and death for \code{\dQuote{CMR}}) accounting for truncation and censoring.}
}

\references{
Burnham, K.P. and Anderson, D.R. (2001) Kullback-Leibler information as a basis for strong inference in ecological studies. \emph{Widlife Research}, 28, 111-119.

Celeux, G., Forbes, F., Robert, C. P., and Titterington, D. M. (2006) Deviance information criteria for missing data models. \emph{Bayesian Analysis}, 1(4), 651-673. 

Colchero, F. and J.S. Clark (2012) Bayesian inference on age-specific survival from capture-recapture data for censored and truncated data. \emph{Journal of Animal Ecology}. 81, 139-149.

Colchero, F., O.R. Jones and M. Rebke. (2012) BaSTA: an R package for Bayesian estimation of age-specific survival from incomplete mark-recapture/recovery data with covariates. \emph{Method in Ecology and Evolution}. 3, 466-470. 

Colchero, F., \emph{et al.}  (2021) The long lives of primates and the "invariant rate of aging" hypothesis. \emph{Nature Communications} 12:3666

Cox, D. R., and Oakes D. (1984) \emph{Analysis of Survival Data}. Chapman and Hall, London.

Gelman, A., Carlin, J.B., Stern, H.S. and Rubin, D.B. (2004) \emph{Bayesian data analysis}. 2nd edn. Chapman & Hall/CRC, Boca Raton, Florida, USA.

Gompertz, B. (1825) On the nature of the function expressive of the law of human mortality, and on a new mode of determining the value of life contingencies. \emph{Philosophical Transactions of the Royal Society of London}, 115, 513-583.

King, R. and Brooks, S.P. (2002) Bayesian model discrimination for multiple strata capture-recapture data. \emph{Biometrika}, 89, 785-806. 

McCulloch, R.E. (1989) Local model influence. \emph{Journal of the American Statistical Association}, 84, 473-478.

Pinder III, J.E., Wiener, J.G. and Smith, M.H. (1978) The Weibull distribution: a new method of summarizing survivorship data. \emph{Ecology}, 59, 175-179.

Spiegelhalter, D.J., Best, N.G., Carlin, B.P. and van der Linde, A. (2002) Bayesian measures of model complexity and fit. \emph{Journal of the Royal Statistical Society: Series B}, 64, 583-639.
}

\author{Fernando Colchero \email{colchero@imada.sdu.dk}}

\seealso{
\code{\link{summary.basta}}, \code{print.basta}, \code{plot.basta} to visualise summary outputs for objects of class \dQuote{\code{basta}}. 

\code{\link{bastaCMRdat}} and \code{\link{bastaCensDat}} for examples of input CMR and census datasets, respectively.

%%\code{\link{CensusToCaptHist}} for raw capture-recapture matrix formatting.
}

\examples{
## ---------- #
## CMR data:
## ---------- #
## Load data:
data("bastaCMRdat", package = "BaSTA")

## Check data consistency:
checkedData  <- DataCheck(bastaCMRdat, dataType = "CMR", studyStart = 51, 
                          studyEnd = 70)

## Run short version of BaSTA on the data:
out <- basta(bastaCMRdat, studyStart = 51, studyEnd = 70, niter = 200, 
             burnin = 11, thinning = 10, updateJumps = FALSE)

## ------------- #
## Census data:
## ------------- #
## Load data:
data("bastaCensDat", package = "BaSTA")

## Check data consistency:
checkedData  <- DataCheck(bastaCensDat, dataType = "census")

## Run short version of BaSTA on the data:
out <- basta(bastaCensDat, niter = 200, burnin = 11, thinning = 10, 
             updateJumps = FALSE)

## --------------------- #
## Check BaSTA outputs:
## --------------------- #
## Print results:
summary(out, digits = 3)

## Plot traces for survival parameters:
plot(out)

## Plot posterior densities of survival parameters:
plot(out, densities = TRUE)

## Plot traces for proportional hazards parameter:
plot(out, trace.name = "gamma")

## Plot survival and mortality curves:
plot(out, plot.type = "demorates")

}

\keyword{Methods}

