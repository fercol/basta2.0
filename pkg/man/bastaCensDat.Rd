\name{bastaCensDat}
\alias{bastaCensDat}
\docType{data}
\title{
Example of census data for BaSTA analysis.
}
\description{
Simulated census data (i.e., continuous observation of individuals) for Bayesian Survival Trajectory Analysis (BaSTA).
}
\usage{data("bastaCensDat")}
\format{
  A data frame with 500 observations on the following 8 variables.
  \describe{
    \item{\code{ID}}{ID for each individual.}
    \item{\code{Birth.Date}}{Dates of birth, formated as \dQuote{\code{YYYY-mm-dd}}.}
    \item{\code{Min.Birth.Date}}{Minimum estimated dates of birth, formated as \dQuote{\code{YYYY-mm-dd}}. If the date of birth is known, then \code{Min.Birth.Date} is equal to \code{Birth.Date}.}
    \item{\code{Max.Birth.Date}}{Maximum estimated dates of birth, formated as \dQuote{\code{YYYY-mm-dd}}. If the date of birth is known, then \code{Max.Birth.Date} is equal to \code{Birth.Date}}
    \item{\code{Entry.Date}}{Dates of entry to the study, formated as \dQuote{\code{YYYY-mm-dd}}.}
    \item{\code{Depart.Date}}{Dates of departure from the study, formated as \dQuote{\code{YYYY-mm-dd}}.}
    \item{\code{Depart.Type}}{a character vector indicating whether the \code{Depart.Date} is because of death (i.e., \dQuote{\code{D}}) or censored (i.e., \dQuote{\code{C}}).}
    \item{\code{Sex}}{a character vector indicating the sex covariate.}
  }
}
\details{
This dataset was created by stochastically simulating a hypothetical population with different mortality patterns between males and females and with proportional decreases in mortality as a function of a hypothetical continuous covariate (e.g. birth weight, average adult weight, etc.). The population was simulated for 40 years, with uniform times of birth within this period. Sex ratios were assumed to be 1:1.  The time of death for each individual was inversed sampled from a Gompertz CDF of ages at death. The Gompertz parameters for females were: \eqn{b_0 = -3} and \eqn{b_1 = 0.15}; and for males at \eqn{b_0 = -2} and \eqn{b_1 = 0.2}.

The resulting dataset includes individuals where the data are left-truncated and/or right-censored. This is typical of capture mark recovery datasets.
}
\examples{
## Load data:
data("bastaCensDat", package = "BaSTA")

## Check data consistency:
checkedData  <- DataCheck(bastaCensDat, dataType = "census")

}
\keyword{datasets}
