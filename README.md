# BaSTA 2.0.1: Bayesian Survival Trajectory Analysis.

BaSTA is an R package for parametric Bayesian estimation of age-specific survival and mortality for left-truncated and right-censored data.

## What's new in version 2.0.1?

This is the next installment of the R package BaSTA, with several new and improved features, which include:  

- Ability to analyze capture-mark-recapture (CMR) and census data;
- Improved life table estimation using product limit estimators on the mean estimated ages at death;
- Plots for goodness of fit on survival and mortality;
- Allows users to specify minimum and maximum birth and death to improve estimation;
- Includes summary statistics such as remaining life expectancy, lifespan inequality, lifespan equality;
- Allows to produce `fancy` plots as in version 1.9.5;
- Includes function `multibasta()` to test multiple mortality functions.

## How to install BaSTA version 2.0.1?
BaSTA version 2.0.1 is now available in CRAN. To install it from CRAN, simply type to the console
```R
install.packages("BaSTA")
```

## How to install development versions of BaSTA?
To install the latest development version of BaSTA from GitHub, type the following lines of code on the R console:

```R
# Install and load 'devtools':
install.packages("devtools")
library(devtools)

# Install BaSTA:
install_git("https://github.com/fercol/basta2.0", subdir = "pkg/")
```
