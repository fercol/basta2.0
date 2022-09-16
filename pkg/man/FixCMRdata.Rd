\name{FixCMRdata}
\alias{FixCMRdata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fix issues on CMR input data for BaSTA.
}
\description{
This function provides general fixes to common issues of capture-mark-recapture (CMR) data.
}
\usage{
FixCMRdata(object, studyStart, studyEnd, autofix = rep(0, 6), 
                       silent = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object }{A data.frame to be used as an input data file for BaSTA for \code{dataType =} \code{\dQuote{CMR}}. The first column is the individual's ID, the second and third columns are birth and death years respectively. Columns 4 to \code{nt+3} represent the observation window of nt years. This is followed (optionally) by columns for covariate. }
	\item{ studyStart }{An integer indicating the first year of the study.}

	\item{studyEnd }{An integer indicating the last year of the study.}

  \item{autofix }{A vector argument with a length of 6 indicating whether to automatically fix any errors (see \code{details}). This should be used with extreme caution. We recommend going back to the individual-based data and fixing each error \dQuote{by hand}.}
  \item{silent }{A logical argument indicating whether to print a detailed report to the screen or not.}
}
\details{
Argument \code{autofix} allows the user to fix the potential errors by specifying a code for each fix. Below are the descriptions of the actions that are taken depending on the error type and the fix code:

Type 1: 0 = do nothing; 1 = remove from dataframe.

Type 2: 0 = do nothing; 1 = remove from dataframe.

Type 3: 0 = do nothing; 1 = replace death records with 0; 2 = replace birth records with 0; 3 = replace both birth and death records with 0.

Type 4: 0 = do nothing; 1 = remove spurious post-death observations.

Type 5: 0 = do nothing; 1 = remove observations that pre-date year of birth.

Type 6: 0 = do nothing; 1 = replace birth year element of observation matrix with 0.
}
\value{
For \code{dataType = \dQuote{CMR}:}
  \item{newData }{A corrected data frame.}
  \item{type1 }{A vector of row numbers in the original data frame where there are deaths occurring before the study starts.}
  \item{type2 }{A vector of row numbers in the original data frame where there are no birth/death AND no obervations.}
  \item{type3 }{A vector of row numbers in the original data frame where there are births recorded after death.}
  \item{type4 }{A vector of row numbers in the original data frame where there are observations (i.e. recaptures) after death.}
  \item{type5 }{A vector of row numbers in the original data frame where there are observations (i.e. recaptures) before birth.}
  \item{type6 }{A vector of row numbers in the original data frame where the year of birth is not a zero in the recapture matrix.}
    \item{stopExec }{Logical that indicates if the data are free of errors or not. i.e. \code{TRUE} = the data have no apparent errors, and \code{FALSE} = there is at leat one error.}

}

\author{Fernando Colchero \email{colchero@imada.sdu.dk}}

\seealso{\code{\link{DataCheck}} for running a data check on the input data for function \code{\link{basta}}.}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x) 
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
