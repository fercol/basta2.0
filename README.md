# BaSTA 2.0: Bayesian Survival Trajectory Analysis (Beta).

BaSTA is an R package for parametric Bayesian estimation of age-specific survival and mortality for left-truncated and 
right-censored data.

## What's new in version 2.0?

This is the next installment of the R package BaSTA, with several new and improved features, which include:  

- Ability to analyze capture-mark-recapture (CMR) and census data;
- Improved life table estimation using product limit estimators on the mean estimated ages at death;
- Plots for goodness of fit on survival and mortality;
- Allows users to specify minimum and maximum birth and death to improve estimation;

## How to install BaSTA2.0?
To install BaSTA2.0 from GitHub, type the following lines of code on the R console:

>`install.packages("devtools")`
> 
> `library(devtools)`
> 
> `install_git("https://github.com/fercol/basta2.0", subdir = "pkg/")`
