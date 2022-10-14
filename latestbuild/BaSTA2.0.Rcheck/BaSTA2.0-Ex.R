pkgname <- "BaSTA2.0"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('BaSTA2.0')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("CensusToCaptHist")
### * CensusToCaptHist

flush(stderr()); flush(stdout())

### Name: CensusToCaptHist
### Title: Constructs a capture-history matrix from repeated individual
###   observations to be used in Bayesian Survival Trajectory Analysis
###   (BaSTA).
### Aliases: CensusToCaptHist

### ** Examples

## Create a simulated vector of repeated IDs:
IDvec <- sort(sample(1:5, size = 15, replace = TRUE))

## Simulate dates (e.g., years) of observation per individual:
dVec <- rep(0, length(IDvec))
for(i in unique(IDvec)) {
  svec <- which(IDvec == i)
  dVec[svec] <- sort(sample(1990:1995, length(svec)))
}

## Construct the capture-recapture matrix:
Y <- CensusToCaptHist(ID = IDvec,  d = dVec)



cleanEx()
nameEx("DataCheck")
### * DataCheck

flush(stderr()); flush(stdout())

### Name: DataCheck
### Title: Error checking for BaSTA input data.
### Aliases: DataCheck
### Keywords: FILL UP

### ** Examples

## CMR data:
## --------- #
## Load data:
data("bastaCMRdat", package = "BaSTA2.0")

## Check data consistency:
checkedData  <- DataCheck(bastaCMRdat, dataType = "CMR", studyStart = 51, 
                          studyEnd = 70)

## census data:
## ------------ #
## Load data:
data("bastaCensDat", package = "BaSTA2.0")

## Check data consistency:
checkedData  <- DataCheck(object = bastaCensDat, dataType = "census")

## Printed output:
## --------------- #
## Print DataCheck results:
print(checkedData)




cleanEx()
nameEx("FixCMRdata")
### * FixCMRdata

flush(stderr()); flush(stdout())

### Name: FixCMRdata
### Title: Fix issues on CMR input data for BaSTA.
### Aliases: FixCMRdata

### ** Examples

## Load data:
data("bastaCMRdat", package = "BaSTA2.0")

## Fix data:
fixedData  <- FixCMRdata(bastaCMRdat, studyStart = 51, 
                          studyEnd = 70, autofix = rep(1, 6))



cleanEx()
nameEx("basta")
### * basta

flush(stderr()); flush(stdout())

### Name: basta
### Title: Parametric Bayesian estimation of age-specific survival for
###   left-truncated and right-censored capture-mark-recapture or census
###   data.
### Aliases: basta basta.default
### Keywords: Methods

### ** Examples

## ---------- #
## CMR data:
## ---------- #
## Load data:
data("bastaCMRdat", package = "BaSTA2.0")

## Check data consistency:
checkedData  <- DataCheck(bastaCMRdat, dataType = "CMR", studyStart = 51, 
                          studyEnd = 70)

## Run short version of BaSTA on the data:
out <- basta(bastaCMRdat, studyStart = 51, studyEnd = 70, niter = 100, 
             burnin = 11, thinning = 10, updateJumps = FALSE)

## ------------- #
## Census data:
## ------------- #
## Load data:
data("bastaCensDat", package = "BaSTA2.0")

## Check data consistency:
checkedData  <- DataCheck(bastaCensDat, dataType = "census")

## Run short version of BaSTA on the data:
out <- basta(bastaCensDat, dataType = "census", niter = 100, burnin = 11, 
             thinning = 10, updateJumps = FALSE)

## --------------------- #
## Check BaSTA outputs:
## --------------------- #
## Print results:
summary(out, digits = 3)

## Plot traces for survival parameters:
plot(out)

## Plot posterior densities of survival parameters:
plot(out, densities = TRUE)

## Plot survival and mortality curves:
plot(out, plot.type = "demorates")




cleanEx()
nameEx("bastaCMRdat")
### * bastaCMRdat

flush(stderr()); flush(stdout())

### Name: bastaCMRdat
### Title: Example of capture-mark-recapture data for BaSTA analysis.
### Aliases: bastaCMRdat
### Keywords: datasets

### ** Examples

## Load data:
data("bastaCMRdat", package = "BaSTA2.0")

## Check data consistency:
checkedData  <- DataCheck(bastaCMRdat, dataType = "CMR", studyStart = 51, 
                          studyEnd = 70)



cleanEx()
nameEx("bastaCMRout")
### * bastaCMRout

flush(stderr()); flush(stdout())

### Name: bastaCMRout
### Title: Output from a Bayesian Survival Trajectory Analysis (BaSTA)
###   analysis on a simulated capture-mark-recapture (CMR) dataset.
### Aliases: bastaCMRout
### Keywords: model output

### ** Examples

## Load BaSTA output:
data("bastaCMRout", package = "BaSTA2.0")

## Plot traces for survival parameters:
plot(bastaCMRout)

## Plot posterior densities of survival parameters:
plot(bastaCMRout, densities = TRUE)

## Plot traces for proportional hazards parameter:
plot(bastaCMRout, trace.name = "gamma")

## Plot survival and mortality curves:
plot(bastaCMRout, plot.type = "demorates")



cleanEx()
nameEx("bastaCensDat")
### * bastaCensDat

flush(stderr()); flush(stdout())

### Name: bastaCensDat
### Title: Example of census data for BaSTA analysis.
### Aliases: bastaCensDat
### Keywords: datasets

### ** Examples

## Load data:
data("bastaCensDat", package = "BaSTA2.0")

## Check data consistency:
checkedData  <- DataCheck(bastaCensDat, dataType = "census")




cleanEx()
nameEx("bastaCensOut")
### * bastaCensOut

flush(stderr()); flush(stdout())

### Name: bastaCensOut
### Title: Output from a Bayesian Survival Trajectory Analysis (BaSTA)
###   analysis on a simulated census dataset.
### Aliases: bastaCensOut
### Keywords: model output

### ** Examples

## Load BaSTA output:
data("bastaCensOut", package = "BaSTA2.0")

## Plot traces for survival parameters:
plot(bastaCensOut)

## Plot posterior densities of survival parameters:
plot(bastaCensOut, densities = TRUE)

## Plot survival and mortality curves:
plot(bastaCensOut, plot.type = "demorates")



cleanEx()
nameEx("summary.basta")
### * summary.basta

flush(stderr()); flush(stdout())

### Name: summary.basta
### Title: Summarizing and plotting Bayesian Survival Trajectory Analysis
###   (BaSTA) model outputs.
### Aliases: summary.basta print.basta plot.basta

### ** Examples

## Load BaSTA output:
data("bastaCMRout", package = "BaSTA2.0")

## Print summary output:
summary(bastaCMRout)

## Plot traces for mortality parameters (theta):
plot(bastaCMRout)

## Plot traces for proportional hazards parameters (gamma):
plot(bastaCMRout, trace.name = "gamma")

## Plot traces for recapture probability(ies) (pi):
plot(bastaCMRout, trace.name = "pi")

## Plot predicted mortality and survival:
plot(bastaCMRout, plot.type = "demorates")

## Change the color for each covariate on 
## the predicted vital rates:
plot(bastaCMRout, plot.type = "demorates", 
     col = c("dark green", "dark blue"))

## Change the color and the legend text:
plot(bastaCMRout, plot.type = "demorates", 
     col = c("dark green", "dark blue"),
     names.legend = c("Females", "Males"))

## Plot predicted mortality and survival 
## between 2 and 8 years of age:
plot(bastaCMRout, plot.trace = FALSE, xlim = c(2, 8))

## Plot predicted mortality and survival 
## between 2 and 8 years of age without
## credible intervals:
plot(bastaCMRout, plot.trace = FALSE, xlim = c(2, 8), 
     noCI = TRUE)

## Plot parameter densities and predicted vital  
## rates in the same plot (i.e. fancy):
plot(bastaCMRout, fancy = TRUE)

## Change colors and legend names for the 
## "fancy" plot:
plot(bastaCMRout, fancy = TRUE, col = c("dark green", "dark blue"),
     names.legend = c("Females", "Males"))



cleanEx()
nameEx("summary.bastaCheckCMR")
### * summary.bastaCheckCMR

flush(stderr()); flush(stdout())

### Name: summary.bastaCheckCMR
### Title: Summary of outputs from the data checking function in BaSTA.
### Aliases: summary.bastaCheckCMR summary.bastaCheckCens
###   print.bastaCheckCMR print.bastaCheckCens
### Keywords: FILL UP

### ** Examples

## CMR data:
## --------- #
## Load data:
data("bastaCMRdat", package = "BaSTA2.0")

## Check data consistency:
checkedData  <- DataCheck(bastaCMRdat, dataType = "CMR", studyStart = 51, 
                          studyEnd = 70)

## census data:
## ------------ #
## Load data:
data("bastaCensDat", package = "BaSTA2.0")

## Check data consistency:
checkedData  <- DataCheck(object = bastaCensDat, dataType = "census")

## Printed output:
## --------------- #
## Print DataCheck results:
print(checkedData)



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
