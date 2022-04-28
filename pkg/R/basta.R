# ============================== CODE METADATA =============================== #
# Functions to analyse mortality from capture-mark-recapture (CMR)
# or from census data 
# Name:        basta.R
# Date:        13 JAN 2021
# Version:     2.0
# Created by:  Fernando Colchero
# Modified by: Fernando Colchero

# General Comments: 
# -----------------
# - Combine packages BaSTA and BaSTA.ZIMS 
# - Extend the CMR part of the package to allow for covariates in recapture
#   probabilities.
# - Allow the CMR part of the package to have Min. and Max. Birth dates
# - Allow the CMR part to have times of censoring before the end of the 
#   study.

# List of names to change:
# ------------------------
# 1.- dataObj to dataObj
# 2.- parObj$theta$priorMu to parObj$theta$priorMean
# 3.- parObj$gamma$priorMu to parObj$theta$priorMean
# 4.- parObj$eta$priorMu to parObj$lambda$priorMean

# ================================ CODE START ================================ #
# =========================== #
# === A) USER FUNCTIONS: ====
# =========================== #

# A.1) Data check function:
# ------------------------- #



# A.2) main basta function:
# ------------------------- #
basta <-
  function(object, ... ) UseMethod("basta")

basta.default <- function(object, dataType = "CMR", 
                          model = "GO", shape = "simple", 
                          studyStart = NULL, studyEnd = NULL, minAge = 0,
                          covarsStruct = "fused", formulaMort = NULL, 
                          formulaRecap = NULL, recaptTrans = studyStart, 
                          niter = 11000, burnin = 1001, thinning = 20, 
                          nsim = 1, parallel = FALSE, ncpus = 2, 
                          lifeTable = TRUE, updateJumps = TRUE, 
                          negSenescence = FALSE, ...) {
  # Create BaSTA environment:
  # bastaenv <- new.env()
  
  # Additional arguments:
  argList <- list(...)
  
  # Check dataset and stop if problems found:
  # dcheck <- DataCheck(object, silent = TRUE)
  dcheck <- list(stopExec = FALSE)
  if (dcheck$stopExec) {
    stop("\nProblems detected with the data.\n", 
         "Use DataCheck() to find problems.", call. = FALSE)
  }
  
  # Variables to pass to all cpus:
  intVars <- c("algObj", "dataObj", "ageObj", "covObj", "defTheta", 
               "fullParObj", "parObj", "parCovObj", 
               "updateJumps", "niter", "nsim", "burnin", "thinning", "minAge", 
               ".CalcMort", '.CalcSurv', ".CalcMort.numeric", 
               ".CalcMort.matrix", ".CalcSurv.numeric", ".CalcSurv.matrix", 
               ".CalcCumHaz", ".CalcCumHaz.matrix", ".CalcCumHaz.numeric")
  
  # ------------------- #
  # General data setup:
  # ------------------- #
  # Algorithm information:
  algObj <- .CreateAlgObj(dataType, model, shape, studyStart, studyEnd, 
                          minAge, covarsStruct, formulaMort, formulaRecap, 
                          recaptTrans, niter, burnin, thinning, updateJumps,
                          nsim, negSenescence)
  
  # Dataset object:
  dataObj <- .CreateDataObj(object, algObj)
  
  # Covariate object:
  covObj <- .CreateCovObj(object, dataObj, algObj)
  
  # Create covariate names object:
  covsNames <- list(cat = NA, con = NA, class = NA)
  if (inherits(covObj, c("cateCov", "bothCov"))) {
    covsNames$cat <- covObj$cat
  }
  if (inherits(covObj, c("contCov", "bothCov"))) {
    covsNames$con <- covObj$cont
  }
  covsNames$class <- class(covObj)[1]
  
  # Define theta parameters:
  defTheta <- .SetDefaultTheta(algObj)
  
  # Incorporate user parameters:
  userPars <- .CreateUserPar(covObj, argList)
  
  # Construct full parameter object:
  fullParObj <- .CreateFullParObj(covObj, defTheta, algObj, userPars, 
                                  dataObj)
  
  # Initial parameter object:
  parObj <- .CreateParObj(fullParObj)
  
  # Covariate-parameter matrices:
  parCovObj <- .CalcCovPars(parObj = parObj, parCov = NULL, 
                            covObj = covObj, dataObj = dataObj, 
                            type = "both")
  
  # Age object:
  ageObj <- .CreateAgeObj(dataObj, algObj)
  
  # Index of kept MCMC records:
  keep <- seq(burnin, niter, thinning)
  
  # ---------------------- #
  # Demographic functions:
  # ---------------------- #
  # a) Mortality function:
  .CalcMort <- function(theta, ...) UseMethod(".CalcMort")
  .CalcMort.matrix <- .DefineMortMatrix(algObj)
  .CalcMort.numeric <- .DefineMortNumeric(algObj)
  
  # b) Cummulative hazard:
  .CalcCumHaz <- function(theta, ...) UseMethod(".CalcCumHaz")
  .CalcCumHaz.matrix <- .DefineCumHazMatrix(algObj)
  .CalcCumHaz.numeric <- .DefineCumHazNumeric(algObj)
  
  # c) Survival:
  .CalcSurv <- function(theta, ...) UseMethod(".CalcSurv")
  .CalcSurv.matrix <- .DefineSurvMatrix(algObj)
  .CalcSurv.numeric <- .DefineSurvNumeric(algObj)
  
  # --------------------------------------- #
  # Assign variables to global environment:
  # --------------------------------------- #
  for (name in intVars) assign(name, get(name), envir = globalenv())
  
  # ------------------------------ #
  # Initial likelihood and priors:
  # ------------------------------ #
  # Likelihood:
  likeObj <- .CalcLike(dataObj, parCovObj, parObj, fullParObj, ageObj, 
                       likeObj = NULL, ind = 1:dataObj$n)
  
  # Priors age object:
  priorAgeObj <- .SetPriorAgeDist(fullParObj, dataObj, covObj)
  
  # include in intVars:
  intVars <- c(intVars, "priorAgeObj")
  assign("priorAgeObj", get("priorAgeObj"), envir = globalenv())
  
  # ----------------------- #
  # Run MCMC to find jumps:
  # ----------------------- #
  cat("\nRunning sequence to find jump SDs... ")
  Start <- Sys.time()
  jumpRun <- .RunMCMC(1, UpdJumps = TRUE, parJumps = NA)
  End <- Sys.time()
  cat("Done\n")
  compTime <- round(as.numeric(End-Start, units = units(End - Start)), 2)
  cat(sprintf("Total jump SDs computing time: %.2f %s.\n\n", compTime, 
              units(End - Start)))
  
  
  # -------------- #
  # Run main MCMC:
  # -------------- #
  Start <- Sys.time()
  if (nsim > 1) {
    cat("Multiple simulations started...\n\n") 
    if (parallel) {
      opp <- options()
      options(warn = -1)
      sfInit(parallel = TRUE, cpus = ncpus)
      sfExport(list = c(intVars, ".Random.seed"))
      sfSource(file = "/Users/fernando/FERNANDO/PROJECTS/4.PACKAGES/BaSTA2.0/tests/sourceBaSTA.R")
      # sfLibrary(BaSTA.ZIMS)
      bastaOut <- sfClusterApplyLB(1:nsim, .RunMCMC, UpdJumps = FALSE, 
                                   parJumps = jumpRun$jumps)
      sfRemoveAll(hidden = TRUE)
      sfStop()
      options(opp)
    } else {
      bastaOut <- lapply(1:nsim, .RunMCMC, UpdJumps = FALSE, 
                         parJumps = jumpRun$jumps)
    }
  } else {
    cat("Simulation started...\n\n")
    bastaOut <- lapply(1:nsim, .RunMCMC, UpdJumps = FALSE, 
                       parJumps = jumpRun$jumps)
  }
  End <- Sys.time()
  cat("Simulations finished.\n")
  compTime <- round(as.numeric(End-Start, units = units(End - Start)), 2)
  cat(sprintf("Total MCMC computing time: %.2f %s.\n\n", compTime, 
              units(End - Start)))
  names(bastaOut) <- paste("sim.", 1:nsim, sep = "")
  
  # Calculate summaries:
  bastaSumars <- .ExtractParalOut(bastaOut, keep, fullParObj, covsNames, nsim,
                                  dataObj, algObj, defTheta, .CalcMort, 
                                  .CalcMort.numeric, .CalcMort.matrix, 
                                  .CalcSurv, .CalcSurv.matrix, 
                                  .CalcSurv.numeric, covObj)
  # Create final BaSTA object:
  bastaFinal <- bastaSumars
  
  # Include individual runs in output:
  bastaFinal$runs <- bastaOut
  
  # Summary information for parameters:
  bastaFinal$fullpar <- fullParObj
  bastaFinal$simthe <- defTheta
  
  # Covariate information:
  bastaFinal$covs <- covsNames
  
  # MCMC settings:
  bastaFinal$settings <- c(niter, burnin, thinning, nsim)
  names(bastaFinal$settings) <- c("niter", "burnin", "thinning", "nsim")
  
  # Model specifications:
  bastaFinal$modelSpecs <- 
    c(model, shape, minAge, covarsStruct,
      paste(names(covObj$cat), collapse = ", "), 
      paste(names(covObj$cont), collapse = ", "), dataType)
  names(bastaFinal$modelSpecs) <- c("model", "shape", "min. age", 
                                    "Covar. structure", "Categorical", 
                                    "Continuous", "DataType")
  
  # Life table:
  bastaFinal$lifeTable <- .CalcLifeTable(bastaOut, lifeTable, covObj, algObj,
                                         dataObj)
  
  # Remove variables from global environment:
  rm(list = intVars)
  
  # Assign class:
  class(bastaFinal) <- "basta"
  return(bastaFinal)
}

# A.3) plotting BaSTA outputs:
# ---------------------------- #
plot.basta <- function(x, plot.type = "traces", trace.name = "theta",
                       densities = FALSE, noCIs = FALSE, 
                       addKM = TRUE, ...) {
  args <- list(...)
  nv <- ifelse(plot.type == "traces", x$settings['nsim'], length(x$surv))
  
  if ("col" %in% names(args)) {
    Palette <- args$col
    if (length(Palette) < nv) {
      ncwarn <- ifelse(plot.trace, "simulation", "covariates")
      warning(sprintf("Insufficient number of colors. Not all %s will be displayed.",
                      ncwarn), call. = FALSE)
    }
  } else {
    if (nv <= 9) {
      Palette <- c('#E41A1C', '#377EB8', '#4DAF4A', '#984EA3', '#FF7F00', 
                   '#FFFF33', '#A65628', '#F781BF', '#999999')
    } else {
      Palette <- rainbow(nv)
    }
  }
  if ("lwd" %in% names(args)) {
    lwd <- args$lwd
  } else {
    lwd <- 1
  }
  if ("lty" %in% names(args)) {
    lty <- args$lty
  } else {
    lty <- 1
  }
  op <- par(no.readonly = TRUE)
  # ======= #
  # TRACES:
  # ======= #
  if (plot.type == "traces") {
    nsim <- x$settings["nsim"]
    if (trace.name == "theta") {
      ncat <- length(x$covs$cat)
      pcol <- ifelse(is.null(ncat), 1, ncat)
      prow <- ceiling(x$fullpar$theta$len / pcol)
      npar <- x$simthe$length
      allnpar <- x$fullpar$theta$len
    } else if (trace.name == "gamma") {
      npar <- ncol(x$runs$sim.1$gamma)
      if (is.null(npar)) {
        stop("'gamma' parameters not calculated (not propHaz)", call. = FALSE)
      }
      pcol <- ceiling(npar / 2)
      prow <- ceiling(npar / pcol)
      allnpar <- x$fullpar$gamma$len
    } else if (trace.name == "lambda") {
      npar <- ifelse("lambda" %in% colnames(x$params), 1, NA)
      if (is.na(npar)) {
        stop("'lambda' parameters not calculated (no minAge)", call. = FALSE)
      }
      pcol <- 1
      prow <- 1
      allnpar <- 1
    } else if (trace.name == "pi") {
      if (x$modelSpecs["DataType"] == "CMR") {
        npar <- ncol(x$runs$sim.1$pi)
        pcol <- ceiling(npar / 2)
        prow <- ceiling(npar / pcol)
        allnpar <- x$fullpar$pi$len
      } else {
        stop("Argument 'trace.name' cannot be 'pi'.\n", 
        "No recapture probabilities estimated, dataType is not CMR.",
             call. = FALSE)
      }
    } else {
      stop("Wrong 'trace.name' specified.\n", 
           "Names are 'theta', 'gamma', 'lambda', or 'pi'", call. = FALSE)
    }
    keep <- seq(1, x$setting["niter"], x$settings["thinning"])
    par (mfrow = c(prow, pcol), mar = c(4, 4, 1, 1))
    for (pp in 1:allnpar) {
      ylim <- 
        range(sapply(1:nsim, function(ii) 
          range(x$runs[[ii]][[trace.name]][, pp]))) 
      xlim <- c(0, nrow(x$runs[[1]][[trace.name]]))
      main <- x$fullpar[[trace.name]]$names[pp]
      plot(xlim, ylim, col = NA, xlab = "", ylab = "", 
           main = main)
      for (tt in 1:nsim) {
        lines(keep, x$runs[[tt]][[trace.name]][keep, pp], 
              col = Palette[tt])
      }
    }
    # =========== #
    # DEMO RATES:
    # =========== #
  } else if (plot.type == "demorates") {
    par(mfrow = c(2, 1), mar = c(4, 4, 1, 1)) 
    demvname <- c("Mortality", "Survival")
    names(demvname) <- c("mort", "surv")
    for (demv in c("mort", "surv")) {
      ylim <- c(0, 0)
      if ("xlim" %in% names(args)) {
        xlim <- args$xlim
      } else {
        xlim <- c(0, 0)
      }
      vars <- names(x$mort)
      minAge <- as.numeric(x$modelSpecs["min. age"])
      for (nta in 1:length(x$mort)) {
        cuts <- x$cuts[[nta]]
        ylim <- range(c(ylim, x[[demv]][[nta]][, cuts]), 
                      na.rm = TRUE)
        if (! "xlim" %in% names(args)) {
          xlim <- range(c(xlim, x$x[cuts] + minAge), na.rm = TRUE)
        }
      }
      plot(xlim, ylim, col = NA, xlab = "", ylab = demvname[demv])
      if (minAge > 0) abline(v = minAge, lty = 2)
      nn <- 0
      for (nta in 1:length(x$mort)) {
        nn <- nn + 1
        cuts <- x$cuts[[nta]]
        yy <- x[[demv]][[nta]][, cuts]
        xx <- x$x[cuts]
        if (!noCIs) {
          polygon(c(xx, rev(xx)) + minAge, c(yy[2, ], rev(yy[3, ])), 
                  col = adjustcolor(Palette[nn], alpha.f = 0.25),
                  border = NA)
        }
        lwdd <- ifelse(length(lwd) > 1, lwd[nn], lwd)
        ltyy <- ifelse(length(lty) > 1, lty[nn], lty)
        lines(xx + minAge, yy[1, ], lwd = lwdd, col = Palette[nn], 
              lty = ltyy)
      }
      if (nn > 1 & demv == 'surv') {
        legend('topright', vars, col = Palette[1:nn], pch = 15, bty = 'n')
      }
    }
    # ================ #
    # GOODNESS OF FIT: 
    # ================ #
  } else {
    ncat <- length(x$surv)
    catname <- names(x$surv)
    pcol <- ceiling(ncat / 2)
    prow <- ceiling(ncat / pcol)
    
    par(mfrow = c(prow, pcol), mar = c(4, 4, 1, 1)) 
    for (nta in 1:ncat) {
      ylim <- c(0, 1)
      xlim <- c(0, 0)
      minAge <- as.numeric(x$modelSpecs["min. age"])
      cuts <- x$cuts[[nta]]
      if (! "xlim" %in% names(args)) {
        xlim <- range(c(xlim, x$x[cuts] + minAge), na.rm = TRUE)
      } else {
        xlim <- args$xlim
      }
      if (is.data.frame(x$lifeTable[[nta]])) {
        lifeTab <- x$lifeTable[[nta]]
      } else {
        lifeTab <- x$lifeTable[[nta]]$Mean
      }
      
      plot(xlim, ylim, col = NA, xlab = "", ylab = "Survival", 
           main = catname[nta])
      if (minAge > 0) abline(v = minAge, lty = 2)
      nn <- 0
      cuts <- x$cuts[[nta]]
      yy <- x$surv[[nta]][, cuts]
      xx <- x$x[cuts]
      if (!noCIs) {
        polygon(c(xx, rev(xx)) + minAge, c(yy[2, ], rev(yy[3, ])), 
                col = adjustcolor(Palette[1], alpha.f = 0.25),
                border = NA)
      }
      lwdd <- ifelse(length(lwd) > 1, lwd[1], lwd)
      ltyy <- ifelse(length(lty) > 1, lty[1], lty)
      lines(lifeTab$Ages, lifeTab$lx, type = "s")
      lines(xx + minAge, yy[1, ], lwd = lwdd, col = Palette[1], 
            lty = ltyy)
      if (nta == ncat) {
        legend('topright', c("Kapplan-Meier", "Estimated survival"), 
               col = c(1, Palette[1]), lwd = 2, bty = 'n')
      }
    }
  }
  par(op)
}

# A.4) Printing BaSTA outputs:
# ---------------------------- #
print.basta <- function(x, ...) {
  extraArgs <- list(...)
  if (length(extraArgs) > 0) {
    if (!is.element('digits', names(extraArgs))){
      digits <- 4
    } else {
      digits <- extraArgs$digits
    }
  } else {
    digits <- 4
  }
  if ("ModelSpecs" %in% names(x)) {
    x$modelSpecs <- x$ModelSpecs
  }
  cat("\nCall:\n")
  cat(paste("Model             \t\t: ", x$modelSpecs[1], "\n", sep = ""))
  cat(paste("Shape             \t\t: ", x$modelSpecs[2], "\n", sep = ""))
  cat(paste("Minimum age       \t\t: ", x$modelSpecs[3], "\n", sep = ""))
  cat(paste("Covars. structure \t\t: ", x$modelSpecs[4], "\n", sep = ""))
  cat(paste("Cat. covars.      \t\t: ", x$modelSpecs[5], "\n", sep = ""))
  cat(paste("Cont. covars.     \t\t: ", x$modelSpecs[6], "\n", 
            collapse = ""))
  
  cat("\nCoefficients:\n")
  print.default(x$coefficients, digits, ...)
  cat("\nConvergence:\n")
  cat(x$convmessage)
  if (is.na(x$DIC[1])){
    cat("\nDIC not calculated.")
  } else {
    cat(sprintf("\nDIC = %s", round(x$DIC["DIC"], 2)))
  }
}

# A.5) Summary for BaSTA outputs:
# ------------------------------- #
summary.basta <-
  function(object,...){
    extraArgs       <- list(...)
    if (length(extraArgs) > 0) {
      if (!is.element('digits', names(extraArgs))){
        digits <- 4
      } else {
        digits <- extraArgs$digits
      }
    } else {
      digits <- 4
    }
    if ("ModelSpecs" %in% names(object)) {
      object$modelSpecs <- object$ModelSpecs
    }
    if ("version" %in% names(object)) {
      cat(sprintf("\nOutput from BaSTA version %s\n", object$version))
    }
    cat("\nCall:\n")
    cat(paste("Model             \t\t: ", object$modelSpecs[1], "\n", sep = ""))
    cat(paste("Shape             \t\t: ", object$modelSpecs[2], "\n", sep = ""))
    cat(paste("Covars. structure \t\t: ", object$modelSpecs[3], "\n", sep = ""))
    cat(paste("Minimum age       \t\t: ", object$modelSpecs[4], "\n", sep = ""))
    cat(paste("Cat. covars.      \t\t: ", object$modelSpecs[5], "\n", sep = ""))
    cat(paste("Cont. covars.     \t\t: ", object$modelSpecs[6], "\n", 
              collapse = ""))
    
    cat("\nModel settings:\n")
    print(object$set)
    
    
    cat("\nMean Kullback-Leibler\ndiscrepancy calibration (KLDC):\n")
    if (object$K[1] != "Not calculated") {
      if ("qkl1" %in% names(object$K)) {
        meanKLcalib  <- t((object$K$qkl1 + object$K$qkl2) / 2)
      } else {
        meanKLcalib  <- (object$K$q12 + object$K$q21) / 2
      }
      print.default(meanKLcalib, digits = digits)
    } else {
      if (object$set['nsim'] == 1) {
        cat("KLDC was not calculated due to insufficient number\n",
            " of simulations to estimate convergence.\n")
      } else {
        cat("KLDC was not calculated due to lack of convergence,\n",
            "or because covariates were not included in the model.\n")
      }
    }
    
    
    cat("\nCoefficients:\n")
    print.default(object$coefficients, digits = digits)
    
    cat("\nConvergence:\n")
    if ("Convergence" %in% names(object)) {
      object$convergence <- object$Convergence
    }
    if (object$convergence[1] == "Not calculated") {
      if (object$set['nsim'] == 1) {
        cat("\nConvergence calculations require more than one run.",
            "\nTo estimate potential scale reduction run at least",
            "two simulations.\n")
      } else {
        cat("\nWarning: Convergence not reached for some parameters",
            " (i.e. 'PotScaleReduc' values larger than 1.1).",
            "\nThese estimates should not be used for inference.\n")
      }
    } else {
      if (all(object$convergence[, "Rhat"] < 1.1)) {
        cat("Appropriate convergence reached for all parameters.\n")
      } else {
        cat("\nWarning: Convergence not reached for some parameters",
            " (i.e. 'PotScaleReduc' values larger than 1.1).",
            "\nThese estimates should not be used for inference.\n")
      }
    } 
    cat("\nDIC:\n")
    if (!is.na(object$DIC[1])){
      cat(object$DIC["DIC"],"\n")
      if ("Convergence" %in% names(object)) {
        warning("Model fit in versions older than BaSTA 1.5 had a mistake in the",
                "\ncalculation of DIC values. In case you are interested in\n",
                "comparing the fit of different models, please try to run them\n",
                "with BaSTA versions 1.5 or above.",
                " We apologise for the inconvenience.", call. = FALSE)
      }
    } else {
      if (object$set['nsim'] == 1) {
        cat("DIC was not calculated due to insufficient number",
            "of simulations to estimate convergence.\n")
      } else {
        cat("DIC was not calculated due to lack of convergence.\n")
      }
    }
    ans <- c(list(coefficients = object$coef, DIC = object$DIC,
                  KullbackLeibler = object$KullbackLeibler, 
                  convergence = object$convergence,
                  modelSpecs = object$modelSpecs, settings = object$set))
    return(invisible(ans))
  }

# ============================================ #
# ==== B) FUNCTIONS FOR INTERNAL OBJECTS: ==== 
# ============================================ #
# Algorithm object function:
.CreateAlgObj <- function(dataType, model, shape, studyStart, studyEnd, 
                          minAge, covarsStruct, formulaMort, formulaRecap, 
                          recaptTrans, niter, burnin, thinning, updateJumps,
                          nsim, negSenescence) {
  return(list(dataType = dataType, model = model, shape = shape, 
              start = studyStart, end = studyEnd, minAge = minAge, 
              covStruc = covarsStruct, formulaMort = formulaMort, 
              formulaRecap = formulaRecap, recap = recaptTrans, 
              niter = niter, burnin = burnin, thinning = thinning, 
              updJump = updateJumps, nsim = nsim,
              negSenescence = negSenescence))
}

# Prepare data object:
.CreateDataObj <- function(object, algObj) {
  classDataObj <- c("bastacmr", "ageUpd")
  # Data Object for CMR data type:
  if (algObj$dataType == "CMR") {
    dataObj <- list()
    # Extract study year sequence and length:
    dataObj$study <- algObj$start:algObj$end
    dataObj$studyLen <- length(dataObj$study)
    
    # Number of observations:
    dataObj$n <- nrow(object)
    
    # Recapture matrix:
    Y <- as.matrix(object[, 1:dataObj$studyLen + 3])
    dataObj$Y <- Y
    dataObj$Y[Y > 1] <- 1
    colnames(dataObj$Y) <- dataObj$study
    
    # NOTE: addition of censTime for studies where individuals are censored
    #       before the end of the study (e.g. two studies of different 
    #       duration).
    # Find possible times of censoring before the study end:
    censTime <- rep(algObj$end, dataObj$n)
    for (tt in 1:dataObj$studyLen) {
      idCens <- which(Y[, tt] > 1)
      censTime[idCens] <- dataObj$study[tt]
    }
    dataObj$censTime <- censTime
    
    # Birth - death matrix:
    bd <- as.matrix(object[, 2:3])
    dataObj$bi <- bd[, 1]
    dataObj$di <- bd[, 2]
    bi0 <- which(dataObj$bi == 0 | is.na(dataObj$bi))
    if (length(bi0) > 0) {
      dataObj$idNoB <- bi0
      dataObj$updB <- TRUE
      dataObj$nUpdB <- length(bi0)
    } else {
      dataObj$updB <- FALSE
      dataObj$nUpdB <- 0
    }
    di0 <- which(dataObj$di == 0 | is.na(dataObj$di))
    if (length(di0) > 0) {
      dataObj$idNoD <- di0
      dataObj$updD <- TRUE
      dataObj$nUpdD <- length(di0)
    } else {
      dataObj$updD <- FALSE
      dataObj$nUpdD <- 0
    }
    
    if (!dataObj$updB & !dataObj$updD) {
      classDataObj[2] <- "noAgeUpd"
      dataObj$updA <- FALSE
    } else {
      dataObj$idNoA <- sort(unique(c(dataObj$idNoB, dataObj$idNoD)))
      dataObj$nUpdA <- length(dataObj$idNoA)
      dataObj$updA <- TRUE
      
      # 4.1.2 Calculate first and last time observed 
      #       and total number of times observed:
      ytemp <- t(t(dataObj$Y) * dataObj$study)
      dataObj$lastObs <- c(apply(ytemp, 1, max))
      ytemp[ytemp == 0] <- 10000
      dataObj$firstObs <- c(apply(ytemp, 1, min))
      dataObj$firstObs[dataObj$firstObs == 10000] <- 0
      dataObj$oi <- dataObj$Y %*% rep(1, dataObj$studyLen)
      
      # 4.1.3 Define study duration:
      dataObj$Tm <- matrix(dataObj$study, dataObj$n, 
                           dataObj$studyLen, byrow = TRUE)
      fii <- dataObj$firstObs
      id1 <- which(dataObj$bi > 0 & dataObj$bi >= algObj$start)
      fii[id1] <- dataObj$bi[id1] + 1
      fii[dataObj$bi > 0 & dataObj$bi < algObj$start]  <- algObj$start
      lii <- dataObj$lastObs
      id2 <- which(dataObj$di > 0 & dataObj$di <= dataObj$censTime)
      lii[id2] <- dataObj$di[id2] - 1
      idCens <- which(dataObj$di > 0 & dataObj$di > dataObj$censTime)
      lii[idCens] <- dataObj$censTime[idCens]
      dataObj$obsMat <- .BuildAliveMatrix(fii, lii, dataObj)
      dataObj$obsMat[lii == 0 | fii == 0, ] <- 0
      classDataObj[2] <- "ageUpd"
    }
    dataObj$Dx <- 1 
    classDataObj[1] <- "bastacmr"
  }
  
  # Data Object for Census data type:
  else {
    n <- nrow(object)
    # Calculate Julian times:
    bi <- round(as.numeric(as.Date(object$Birth.Date, 
                                   format = "%Y-%m-%d")) / 
                  365.25, 2) + 1970
    bil <- round(as.numeric(as.Date(object$Min.Birth.Date, 
                                    format = "%Y-%m-%d")) /
                   365.25, 2) + 1970
    biu <- round(as.numeric(as.Date(object$Max.Birth.Date, 
                                    format = "%Y-%m-%d")) /
                   365.25, 2) + 1970
    firstObs <- round(as.numeric(as.Date(object$Entry.Date, 
                                      format = "%Y-%m-%d")) /
                     365.25, 2) + 1970
    lastObs <- round(as.numeric(as.Date(object$Depart.Date, 
                                       format = "%Y-%m-%d")) /
                      365.25, 2) + 1970
    
    # Entry and departure types:
    entryType <- as.character(object$Entry.Type)
    departType <- as.character(object$Depart.Type)
    
    # Censored individuals:
    idCens <- which(departType %in% c("O", "C"))
    nCens <- length(idCens)
    
    # Birth times to be estimated:
    idNoBirth <- which(bi != bil | bi != biu)
    nNoBirth <- length(idNoBirth)
    if (nNoBirth == 0) {
      updB <- FALSE
      classDataObj[2] <- "noAgeUpd"
    } else {
      updB <- TRUE
      classDataObj[2] <- "ageUpd"
    }

    # Create data object:
    dataObj <- list(bi = bi, bil = bil, biu = biu, firstObs = firstObs, 
                    lastObs = lastObs, idCens = idCens,
                    idNoB = idNoBirth, nUpdB = nNoBirth, updB = updB, 
                    updD = FALSE, idNoD = NA, nUpdD = 0,
                    idNoA = idNoBirth, nUpdA = nNoBirth, updA = updB, n = n, 
                    nCens = nCens, studyLen = NA)
    classDataObj[1] <- "bastacensus"
  }
  class(dataObj) <- classDataObj

  # Output:
  return(dataObj)
}

# Create initial age object:
.CreateAgeObj <- function(dataObj, algObj) {
  ageObj <- list()
  bi <- dataObj$bi
  if (dataObj$updB & inherits(dataObj, "bastacmr")) {
    bi[dataObj$idNoB] <- dataObj$firstObs[dataObj$idNoB] - 
      sample(6:1, size = dataObj$nUpdB, replace = TRUE)
  }
  if (inherits(dataObj, "bastacmr")) {
    di <- dataObj$di
    if (dataObj$updD) {
      di[dataObj$idNoD] <- apply(cbind(bi[dataObj$idNoD], 
                                       dataObj$lastObs[dataObj$idNoD]),
                                 1, max) + 
        sample(6:1, size = dataObj$nUpdD, replace = TRUE)
    }
    age <- di - bi
    ageTr <- algObj$start - bi
    ageTr[ageTr < 0] <- 0
  } else {
    di <- dataObj$lastObs
    age <- dataObj$lastObs - bi
    ageTr <- dataObj$firstObs - bi
    ageTr[ageTr < 0] <- 0
  }
  
  # indicator for uncensored:
  indUncens <- rep(1, dataObj$n)
  if (inherits(dataObj, "bastacensus")) {
    indUncens <- rep(1, dataObj$n)
    indUncens[dataObj$idCens] <- 0
  }
  
  # Recalculate ages based on minAge:
  # ---------------------------------
  if (algObj$minAge > 0) {
    # Ages and indicator after min age:
    ageAftMa <- age - algObj$minAge
    ageAftMa[ageAftMa < 0] <- 0
    
    # Ages and indicator before min age:
    ageBefMa <- age
    indBefMa <- rep(0, dataObj$n)
    indBefMa[ageBefMa < algObj$minAge] <- 1
    ageBefMa[age >= algObj$minAge] <- algObj$minAge
    
    # Ages at trunction and indicator after min age:
    ageTrAftMa <- ageTr - algObj$minAge
    ageTrAftMa[ageTrAftMa < 0] <- 0
    
    # Ages at trunction and indicator before min age:
    ageTrBefMa <- ageTr
    ageTrBefMa[ageTr >= algObj$minAge] <- algObj$minAge
  } else {
    ageBefMa <- ageTrBefMa <- indBefMa <- rep(0, dataObj$n)
    ageAftMa <- age
    ageTrAftMa <- ageTr
  }
  
  
  # Create alive matrix for bastacmr:
  # ---------------------------------
  if (inherits(dataObj, "bastacmr")) {
    firstObs <- c(apply(cbind(algObj$start, bi + 1), 1, max))
    lastObs <- c(apply(cbind(algObj$end, dataObj$censTime, di), 1, min))
    alive <- .BuildAliveMatrix(firstObs, lastObs, dataObj)
  } else {
    alive <- NA
  }
  
  # Fill-in matrices for age object:
  # --------------------------------
  ageObj$ages <- data.frame(birth = bi, death = di, age = age, ageTr = ageTr,
                            ageAft = ageAftMa, ageBef = ageBefMa, 
                            truAft = ageTrAftMa, truBef = ageTrBefMa)
  ageObj$inds <- data.frame(ageBef = indBefMa, uncens = indUncens)
  ageObj$alive <- alive
  
  # Assign class to age object:
  # ---------------------------
  minAgeClass <- ifelse(algObj$minAge > 0, "minAge", "noMinAge")
  if (inherits(dataObj, "bastacmr")) {
    class(ageObj) <- c("agecmr", minAgeClass)
  } else {
    class(ageObj) <- c("agecensus", minAgeClass)
  }
  return(ageObj)
}

# Build a matrix with 1 when alive, 0 otherwise:
.BuildAliveMatrix <- function(f, l, dataObj) {
  Fm <- dataObj$Tm - f
  Fm[Fm >= 0] <- 1
  Fm[Fm < 0] <- 0
  Lm <- dataObj$Tm - l
  Lm[Lm <= 0] <- -1
  Lm[Lm > 0] <- 0
  return(Fm * (-Lm))
}

# Create covariates object:
.CreateCovObj <- function(object, dataObj, algObj) {
  covObj <- list()
  # ------------------------- #
  # Covariates for mortality:
  # ------------------------- #
  covClass <- c("noCov", "noCovType", "noCovRecap")
  # Find if mortality covariates are not needed:
  # -------------------------------------------- #
  if (is.null(algObj$formulaMort) & algObj$dataType == "CMR" & 
      ncol(object) == dataObj$studyLen + 3) {
    covObj$covs <- NULL
  } else if (is.null(algObj$formulaMort) & algObj$dataType == "census") {
    covObj$covs <- NULL
  } else {
    # When there should be covariates for mortality:
    # ---------------------------------------------- #
    # Construct covariate matrix:
    if (is.null(algObj$formulaMort) & ncol(object) > dataObj$studyLen + 3) {
      covMat <- as.matrix(object[, (dataObj$studyLen + 4):ncol(object)])
      mode(covMat) <- "numeric"
      colnames(covMat) <- colnames(object)[(dataObj$studyLen + 4):ncol(object)]
    } else {
      covMat <- .MakeCovMat(algObj$formulaMort, data = object)
    }
    # Find covariate types (categorical vs continuous):
    covType <- .FindCovType(covMat)
    
    # Separate covariates by their type and by covarsStruct:
    if (algObj$covStruc == "fused") {
      covClass[1] <- "fused"
      if (!is.null(covType$cat)) {
        covObj$inMort <- covMat[, covType$cat]
        covObj$imLen <- ncol(covObj$inMort)
      } else {
        covClass[1] <- "propHaz"
      }
      if (!is.null(covType$cont)) {
        covObj$propHaz <- matrix(covMat[, c(covType$int, covType$cont)], 
                                 ncol = length(c(covType$int, covType$cont)),
                                 dimnames = list(NULL, c(names(covType$int), 
                                                         names(covType$cont))))
        covObj$phLen <- ncol(covObj$propHaz)
      } else {
        covClass[1] <- "inMort"
      }
    } else if (algObj$covStruc == "all.in.mort") {
      if (is.null(covType$int) & is.null(covType$cat)) {
        covObj$inMort <- cbind(1, covMat)
        colnames(covObj$inMort) <- c("Intercept", colnames(covMat))
      } else {
        covObj$inMort <- covMat
      }
      covObj$imLen <- ncol(covObj$inMort)
      covClass[1] <- "inMort"
    } else {
      if (!is.null(covType$int)) {
        covObj$propHaz <- 
          matrix(covMat[, -covType$int], dataObj$n, ncol(covMat) -1, 
                 dimnames = list(NULL, colnames(covMat)[-covType$int]))
      } else if (!is.null(covType$cat)) {
        covObj$propHaz <- 
          matrix(covMat[, -covType$cat[1]], dataObj$n, ncol(covMat) -1, 
                 dimnames = list(NULL, colnames(covMat)[-covType$cat[1]]))
      } else {
        covObj$propHaz <- covMat
      }
      covObj$phLen <- ncol(covObj$propHaz)
      covClass[1] <- "propHaz"
    }
    if (!is.null(covType$cat) & !is.null(covType$cont)) {
      covClass[2] <- "bothCov"
      covObj$cat <- covType$cat
      covObj$cont <- covType$cont
    } else if (!is.null(covType$cat)) {
      covClass[2] <- "cateCov"
      covObj$cat <- covType$cat
    } else if (!is.null(covType$cont)) {
      covClass[2] <- "contCov"
      covObj$cont <- covType$cont
    }
  }
  
  # ------------------------------------- #
  # Covariates for recapture probability:
  # ------------------------------------- #
  
  # ------------------------ #
  # Return covariate object:
  # ------------------------ #
  class(covObj) <- covClass
  return(covObj)
}

# Find covariate type (i.e. continuous vs categorical):
.FindCovType <- function(covMat) {
  # This functions finds and returns if an intercecpt was included 
  # and which covariates are categorical or continuous.
  if (!is.null(covMat)) {
    lu <- apply(covMat, 2, function(x) length(unique(x)))
    ru <- apply(covMat, 2, range)
    idcat <- which(lu == 2 & apply(ru, 2, sum) == 1)
    if (length(idcat) == 0) {
      idcat <- NULL
    }
    idint <- which(lu == 1)
    if (length(idint) == 0) {
      idint <- NULL
    }
    idcon <- which(lu > 2)
    if (length(idcon) == 0) {
      idcon <- NULL
    }
  }
  else {
    idcat <- NULL
    idint <- NULL
    idcon <- NULL
  }
  return(list(int = idint, cat = idcat, cont = idcon))
}

# Construct a design matrix (i.e. covariate matrix) following a formula:
.MakeCovMat <- function(covform, data) {
  covs <- model.matrix(covform, data = data)
  return(covs)
}

# ------------------------------- #
# Functions to manage parameters:
# ------------------------------- #
# Extract parameters specified by the user:
.CreateUserPar <- function(covObj, argList) {
  userPars <- list()
  genParName <- c("theta", "gamma")
  parTypes <- c("Start", "Jumps", "PriorMean", "PriorSd")
  parTypesList <- c("start", "jump", "priorMean", "priorSd")
  for (genPp in 1:2) {
    if (all(genPp == 2 & inherits(covObj, c("fused", "propHaz"))) |
        genPp == 1) {
      userPars[[genParName[genPp]]] <- list()
      for (pp in 1:4) {
        usrPar <- sprintf("%s%s", genParName[genPp], parTypes[pp])
        if (usrPar %in% names(argList)) {
          userPars[[genParName[genPp]]][[parTypesList[pp]]] <- 
            argList[[usrPar]]
        } else {
          userPars[[genParName[genPp]]][[parTypesList[pp]]] <- NULL 
        }
      }
    }    
  }
  return(userPars)
}

# Set the default mortality parameters:
.SetDefaultTheta  <- function(algObj) {
  if (algObj$model == "EX") {
    nTh <- 1
    startTh <- 0.2 
    jumpTh <- 0.1
    priorMean <- 0.06
    priorSd <- 1
    nameTh <- "b0"
    lowTh <- 0
    jitter <- 0.5
  } else if (algObj$model == "GO") {
    nTh <- 2 
    startTh <- c(-2, 0.01) 
    jumpTh <- c(0.1, 0.1)
    priorMean <- c(-3, 0.01)
    priorSd <- c(1, 1)
    nameTh <- c("b0", "b1")
    lowTh <- c(-Inf, 0)
    if (algObj$negSenescence) lowTh[2] <- -Inf
    jitter <- c(0.5, 0.2) 
    if (algObj$shape == "bathtub") {
      lowTh <- c(-Inf, 0)
    }
  } else if (algObj$model == "WE") {
    nTh <- 2
    startTh <- c(1.5, 0.2) 
    jumpTh <- c(.01, 0.1)
    priorMean <- c(1.5, .05)
    priorSd <- c(1, 1)
    nameTh <- c("b0", "b1")
    lowTh <- c(0, 0)
    jitter <- c(0.5, 0.2) 
  } else if (algObj$model == "LO") {
    nTh <- 3 
    startTh <- c(-2, 0.01, 1e-04) 
    jumpTh <- c(0.1, 0.1, 0.1) 
    priorMean <- c(-3, 0.01, 1e-10)
    priorSd <- c(1, 1, 1)
    nameTh <- c("b0", "b1", "b2")
    lowTh <- c(-Inf, 0, 0)
    jitter <- c(0.5, 0.2, 0.5) 
  }
  if (algObj$shape == "Makeham") {
    nTh <- nTh + 1 
    startTh <- c(0, startTh) 
    jumpTh <- c(0.1, jumpTh) 
    priorMean <- c(0, priorMean)
    priorSd <- c(1, priorSd)
    nameTh <- c("c", nameTh)
    lowTh <- c(0, lowTh)
    jitter <- c(0.25, jitter) 
  } else if (algObj$shape == "bathtub") {
    nTh <- nTh + 3 
    startTh <- c(-0.1, 0.6, 0, startTh)
    jumpTh <- c(0.1, 0.1, 0.1, jumpTh) 
    priorMean <- c(-2, 0.01, 0, priorMean)
    priorSd <- c(1, 1, 1, priorSd)
    nameTh <- c("a0", "a1", "c", nameTh)
    lowTh <- c(-Inf, 0, 0, lowTh)
    jitter <- c(0.5, 0.2, 0.2, jitter) 
  }
  defaultTheta  <- list(length = nTh, start = startTh, jump = jumpTh, 
                        priorMean = priorMean, priorSd = priorSd, name = nameTh, 
                        low = lowTh, jitter = jitter)
  attr(defaultTheta, "model") = algObj$model
  attr(defaultTheta, "shape") = algObj$shape
  return(defaultTheta)
}

# Build the parameters object:
.CreateFullParObj <- function(covObj, defTheta, algObj, 
                             userPars, dataObj) {
  fullParObj <- list()
  fullParObj$theta <- list()
  statName <- c("start", "priorMean", "priorSd", "jump", "low", "jitter")
  nstat <- length(statName)
  for (st in 1:nstat) {
    if (inherits(covObj, c("inMort", "fused"))) {
      if (is.null(userPars$theta[[statName[st]]])) {
        thetaMat <- matrix(defTheta[[statName[st]]], covObj$imLen, 
                           defTheta$length, byrow = TRUE, 
                           dimnames = list(colnames(covObj$inMort), 
                                           defTheta$name))
        if (st %in% c(1, 2)) {
          if (inherits(covObj, "inMort") & 
              inherits(covObj, c("contCov", "bothCov"))) {
            thetaMat[names(covObj$cont), ] <- 0
          }
        }
        fullParObj$theta[[statName[[st]]]] <- thetaMat
      } else {
        if (is.element(length(userPars$theta[[statName[st]]]), 
                       c(defTheta$length, defTheta$length * covObj$imLen))) {
          if (length(userPars$theta[[statName[st]]]) == defTheta$length) {
            fullParObj$theta[[statName[[st]]]] <- 
              matrix(userPars$theta[[statName[st]]], covObj$imLen, 
                     defTheta$length, byrow = TRUE, 
                     dimnames = list(colnames(covObj$inMort), defTheta$name))
          } else {
            fullParObj$theta[[statName[[st]]]] <- userPars$theta[[statName[st]]]
            dimnames(fullParObj$theta[[statName[[st]]]]) <- 
              list(colnames(covObj$inMort), defTheta$name)
          }
        } else {
          stop(paste("\nDimensions of theta ", statName[st], 
                     " matrix are incorrect.\n",
                     "Provide a single vector of length ", defTheta$length,
                     "\nor a matrix of dimensions ", covObj$imLen ," times ", 
                     defTheta$length, 
                     ".\n(i.e. number of covariates times number", 
                     " of\n parameters for model ", 
                     algObj$model," with ", algObj$shape, " shape).", 
                     sep = ""), call. = FALSE)
        }
      }
      allstatName <- paste(rep(defTheta$name, 
                               each = ncol(covObj$inMort)), 
                           rep(colnames(covObj$inMort), defTheta$len), 
                           sep = ".")
    } else {
      if (is.null(userPars$theta[[statName[st]]])) {
        fullParObj$theta[[statName[[st]]]] <- 
          matrix(defTheta[[statName[st]]], 1, defTheta$length, 
                 dimnames = list(NULL, defTheta$name))
      } else {
        if (length(userPars$theta[[statName[st]]]) == defTheta$length) {
          fullParObj$theta[[statName[[st]]]] <- 
            matrix(userPars$theta[[statName[st]]], 1, defTheta$length, 
                   dimnames = list(NULL, defTheta$name))
        } else {
          stop(paste("\nLength of theta ", statName[st], " is incorrect.\n",
                     "Provide a single vector of length ", defTheta$length,
                     ".\n(i.e. number of parameters for model ", 
                     algObj$model," with ", algObj$shape, " shape).", 
                     sep = ""), call. = FALSE)
        }
      }
      allstatName <- defTheta$name
    }
  }
  fullParObj$theta$low <- t(t(fullParObj$theta$start) * 0 + defTheta$low)
  fullParObj$theta$len <- length(fullParObj$theta$start)
  fullParObj$theta$names <- allstatName
  if (inherits(covObj, c("propHaz", "fused"))) {
    fullParObj$gamma <- list()
    for (st in 1:nstat) {
      if (is.null(userPars$gamma[[statName[st]]])) {
        fullParObj$gamma[[statName[st]]] <- 
          rep(c(0.01, 0, 1, 0.1, -Inf, 0.2)[st], covObj$phLen)
        names(fullParObj$gamma[[statName[st]]]) <- colnames(covObj$propHaz)
      } else {
        if (length(userPars$gamma[[statName[st]]]) == covObj$phLen) {
          fullParObj$gamma[[statName[st]]] <- userPars$gamma[[statName[st]]]
          names(fullParObj$gamma[[statName[st]]]) <- colnames(covObj$propHaz)
        } else {
          stop(paste("\nLength of gamma parameters is incorrect.\n",
                     "Provide a single vector of length ", covObj$phLen,
                     ".\n(i.e. number of proportional hazards covariates).", 
                     sep = ""), call. = FALSE)
        }
      }
    }
    fullParObj$gamma$len <- length(fullParObj$gamma$start)
    allstatName <- c(allstatName, paste("gamma", colnames(covObj$propHaz), 
                                        sep = "."))
    fullParObj$gamma$names <- paste("gamma", colnames(covObj$propHaz), 
                                    sep = ".")
  }
  if (inherits(covObj, c("inMort", "noCov"))) {
    Classes <- "theta"
  } else {
    Classes <- "theGam"
  }
  # Minimum age's lambda:
  if (algObj$minAge > 0) {
    fullParObj$lambda <- list(start = 0.01, priorMean = 0.01, priorSd = 1, 
                              jump = 0.01, low = 0, jitter = 0.1,  len = 1)
    Classes <- c(Classes, "lambda")
  } else {
    Classes <- c(Classes, "noLambda")
  }
  
  # Detection probability:
  if (inherits(dataObj,  "bastacmr")) {
    if (inherits(dataObj, "ageUpd")) {
      fullParObj$pi <- list()
      study <- algObj$start:algObj$end
      if (length(algObj$recap) == length(study)) {
        if (all(algObj$recap %in% study)) {
          idpi <- rep(1, length(study))
          for (ipi in 2:length(idpi)) {
            if (algObj$recap[ipi] == algObj$recap[ipi - 1]) {
              idpi[ipi] <- idpi[ipi - 1]
            } else if (algObj$recap[ipi] %in% algObj$recap[1:(ipi - 1)]) {
              idpi[ipi] <- idpi[which(algObj$recap[1:(ipi - 1)] ==
                                        algObj$recap[ipi])[1]]
            } else {
              idpi[ipi] <- max(idpi[1:(ipi - 1)] + 1)
            }
          }
        } else if (all(algObj$recap %in% 1:length(study))) {
          idpi <- algObj$recap
        }
        namespi <- unique(algObj$recap)
      } else {
        idpi <- findInterval(study, algObj$recap)
        namespi <- algObj$recap
      }
      names(idpi) <- study
      npi <- length(unique(idpi))
      fullParObj$pi$start <- rep(0.5, npi)
      names(fullParObj$pi$start) <- namespi
      fullParObj$pi$idpi <- idpi
      fullParObj$pi$n <- npi
      fullParObj$pi$prior2 <- 1
      fullParObj$pi$Prior1 <- tapply(1 + t(t(dataObj$Y) %*% rep(1, dataObj$n)),
                                     idpi, sum)
      fullParObj$pi$len <- length(fullParObj$pi$start)
      if (length(namespi) > 1) {
        piNames <- paste("pi", namespi, sep = ".")
      } else {
        piNames <- "pi"
      }
      fullParObj$pi$names <- piNames
      fullParObj$allNames <- c(allstatName, piNames)
      Classes <- c(Classes, "pi", "noEta")
    } else {
      fullParObj$allNames <- allstatName
      Classes <- c(Classes, "noPi", "noEta")
    }
  } else {
    Classes <- c(Classes, "noPi", "noEta")
  }
  fullParObj$class <- Classes
  class(fullParObj) <- Classes
  return(fullParObj)
}

# Define initial set of parameters:
.CreateParObj <- function(fullParObj) {
  parObj <- list()
  parObj$theta <- fullParObj$theta$start
  if (inherits(fullParObj,  "theGam")) {
    parObj$gamma <- fullParObj$gamma$start
  } else {
    parObj$gamma <- 0
  }
  if (inherits(fullParObj,  "lambda")) {
    parObj$lambda <- fullParObj$lambda$start
  } else {
    parObj$lambda <- 0
  }
  if (inherits(fullParObj,  c("pi", "piEta"))) {
    parObj$pi <- fullParObj$pi$start
  }
  class(parObj) <- class(fullParObj)
  return(parObj)
}

# Create and manage design matrices:
.CalcCovPars <- function(parObj, parCov, covObj, dataObj, type = "theta") {
  if (is.null(parCov)) {
    parCovObjn <- list()
  } else {
    parCovObjn <- parCov
  }
  if (type %in% c("theta", "both")) {
    if (inherits(covObj, c("fused", "inMort"))) {
      parCovObjn$theta <- covObj$inMort %*% parObj$theta
    } else {
      parCovObjn$theta <- matrix(1, nrow = dataObj$n) %*% parObj$theta
    }
  }
  if (type %in% c("gamma", "both")) {
    if (inherits(parObj, "theGam")) {
      parCovObjn$gamma <- covObj$propHaz %*% parObj$gamma
    } else {
      parCovObjn$gamma <- rep(0, dataObj$n)
    }
  }
  return(parCovObjn)
}

# ================================ #
# ==== STATISTICAL FUNCTIONS: ====
# ================================ #
# --------------------------------- #
# Truncated distribution functions:
# --------------------------------- #
# Truncated normal:
.rtnorm <- function(n, mean, sd, lower = -Inf, upper = Inf) {
  Flow <- pnorm(lower, mean, sd)
  Fup <- pnorm(upper, mean, sd)
  ru <- runif(n, Flow, Fup)
  rx <- qnorm(ru, mean, sd)
  return(rx)
}

.dtnorm <- function(x, mean, sd, lower = -Inf, upper = Inf, log = FALSE) {
  Flow <- pnorm(lower, mean, sd)
  Fup <- pnorm(upper, mean, sd)
  densx <- dnorm(x, mean, sd) / (Fup - Flow)
  if (log) densx <- log(densx)
  return(densx)
}

.ptnorm <- function(q, mean, sd, lower = -Inf, upper = Inf, log = FALSE) {
  p <- (pnorm(q, mean, sd) - pnorm(lower, mean, sd)) / 
    (pnorm(upper, mean, sd) - pnorm(lower, mean, sd))
  if (log) {
    p <- log(p)
  }
  return(p)
}

.qtnorm <- function (p, mean = 0, sd = 1, lower = -Inf, upper = Inf) {
  p2 <- (p) * (pnorm(upper, mean, sd) - pnorm(lower, mean, sd)) + 
    pnorm(lower, mean, sd)
  q <- qnorm(p2, mean, sd)
  return(q)
}

# Truncated Gamma:
.rtgamma <- function(n, shape, rate = 1, scale = 1/rate, lower = -Inf, 
                     upper = Inf) {
  Flow <- pgamma(lower, shape, scale = scale)
  Fup <- pgamma(upper, shape, scale = scale)
  ru <- runif(n, Flow, Fup)
  rx <- qgamma(ru, shape, scale = scale)
  return(rx)
}

.dtgamma <- function(x, shape, rate = 1, scale = 1/rate, lower = -Inf, 
                     upper = Inf, log = FALSE) {
  Flow <- pgamma(lower, shape, scale = scale)
  Fup <- pgamma(upper, shape, scale = scale)
  densx <- dgamma(x, shape, scale = scale) / (Fup - Flow)
  if (log) densx <- log(densx)
  return(densx)
}

.ptgamma <- function(q, shape, rate = 1, scale = 1/rate, lower = -Inf, 
                     upper = Inf, log = FALSE) {
  p <- (pgamma(q, shape, scale = scale) - 
          pgamma(lower, shape, scale = scale)) / 
    (pgamma(upper, shape, scale = scale) - 
       pgamma(lower, shape, scale = scale))
  if (log) {
    p <- log(p)
  }
  return(p)
}

.qtgamma <- function (p, shape, rate = 1, scale = 1/rate, lower = -Inf, 
                      upper = Inf) {
  p2 <- (p) * (pgamma(upper, shape, scale = scale) - 
                 pgamma(lower, shape, scale = scale)) + 
    pgamma(lower, shape, scale = scale)
  q <- qgamma(p2, shape, scale = scale)
  return(q)
}

# ---------------------------------- #
# Basic survival analysis functions:
# ---------------------------------- #
# a) General mortality functions:
.DefineMortMatrix <- function(algObj) {
  if (algObj$model == "EX") {
    mortfun <- function(theta, x) rep(theta, length(x))
  } else if (algObj$model == "GO") {
    if (algObj$shape == "simple") {
      mortfun <- function(theta, x) {
        exp(theta[ ,"b0"] + theta[, "b1"] * x)
      }
    } else if (algObj$shape == "Makeham") {
      mortfun <- function(theta, x) {
        theta[, "c"] + exp(theta[, "b0"] + theta[, "b1"] * x)
      }
    } else {
      mortfun <- function(theta, x) {
        exp(theta[, "a0"] - theta[, "a1"] * x) + theta[, "c"] + 
          exp(theta[, "b0"] + theta[, "b1"] * x)
      }
    }
  } else if (algObj$model == "WE") {
    if (algObj$shape == "simple") {
      mortfun <- function(theta, x) {
        mu <- theta[, "b0"] * theta[, "b1"]^theta[, "b0"] * 
          x^(theta[, "b0"] - 1)
        id0 <- which(x == 0)
        if (length(id0) > 0) mu[id0] <- min(mu[-id0])
        return(mu)
      }
    } else if (algObj$shape == "Makeham") {
      mortfun <- function(theta, x) {
        mu <- theta[, "c"] + theta[, "b0"] * theta[, "b1"]^theta[, "b0"] * 
          x^(theta[, "b0"] - 1)
        id0 <- which(x == 0)
        if (length(id0) > 0) mu[id0] <- theta[id0, "c"]
        return(mu)      
      }
    } else {
      mortfun <- function(theta, x) {
        mu <- exp(theta[, "a0"] - theta[, "a1"] * x) + theta[, "c"] + 
          theta[, "b0"] * theta[, "b1"]^theta[, "b0"] * 
          x^(theta[, "b0"] - 1)
        id0 <- which(x == 0)
        if (length(id0) > 0) {
          mu[id0] <- exp(theta[id0, "a0"]) + theta[id0, "c"]
        }
        return(mu)
      }
    }
  } else if (algObj$model == "LO") {
    if (algObj$shape == "simple") {
      mortfun <- function(theta, x) {
        exp(theta[, "b0"] + theta[, "b1"] * x) / 
          (1 + theta[, "b2"] * exp(theta[, "b0"]) / 
             theta[, "b1"] * (exp(theta[, "b1"] * x) - 1))
      }
    } else if (algObj$shape == "Makeham") {
      mortfun <- function(theta, x) {
        theta[, "c"] + exp(theta[, "b0"] + theta[, "b1"] * x) / 
          (1 + theta[, "b2"] * exp(theta[, "b0"]) / 
             theta[, "b1"] * (exp(theta[, "b1"] * x) - 1))
      }
    } else {
      mortfun <- function(theta, x) {
        exp(theta[, "a0"] - theta[, "a1"] * x) + theta[, "c"] + 
          exp(theta[, "b0"] + theta[, "b1"] * x) / 
          (1 + theta[, "b2"] * exp(theta[, "b0"]) / 
             theta[, "b1"] * (exp(theta[, "b1"] * x) - 1))
      }
    }
  }
  return(mortfun)
}

.DefineMortNumeric <- function(algObj) {
  if (algObj$model == "EX") {
    mortfun <- function(theta, x) rep(theta, length(x))
  } else if (algObj$model == "GO") {
    if (algObj$shape == "simple") {
      mortfun <- function(theta, x) {
        exp(theta["b0"] + theta["b1"] * x)
      }
    } else if (algObj$shape == "Makeham") {
      mortfun <- function(theta, x) {
        theta["c"] + exp(theta["b0"] + theta["b1"] * x)
      }
    } else {
      mortfun <- function(theta, x) {
        exp(theta["a0"] - theta["a1"] * x) + theta["c"] + 
          exp(theta["b0"] + theta["b1"] * x)
      }
    }
  } else if (algObj$model == "WE") {
    if (algObj$shape == "simple") {
      mortfun <- function(theta, x) {
        mu <- theta["b0"] * theta["b1"]^theta["b0"] * 
          x^(theta["b0"] - 1)
        id0 <- which(x == 0)
        if (length(id0) > 0) mu[id0] <- min(mu[-id0])
        return(mu)
      }
    } else if (algObj$shape == "Makeham") {
      mortfun <- function(theta, x) {
        mu <- theta["c"] + theta["b0"] * theta["b1"]^theta["b0"] * 
          x^(theta["b0"] - 1)
        id0 <- which(x == 0)
        if (length(id0) > 0) mu[id0] <- theta["c"]
        return(mu)
      }
    } else {
      mortfun <- function(theta, x) {
        mu <- exp(theta["a0"] - theta["a1"] * x) + theta["c"] + 
          theta["b0"] * theta["b1"]^theta["b0"] * 
          x^(theta["b0"] - 1)
        id0 <- which(x == 0)
        if (length(id0) > 0) {
          mu[id0] <- exp(theta["a0"]) + theta["c"]
        }
        return(mu)
      }
    }
  } else if (algObj$model == "LO") {
    if (algObj$shape == "simple") {
      mortfun <- function(theta, x) {
        exp(theta["b0"] + theta["b1"] * x) / 
          (1 + theta["b2"] * exp(theta["b0"]) / 
             theta["b1"] * (exp(theta["b1"] * x) - 1))
      }
    } else if (algObj$shape == "Makeham") {
      mortfun <- function(theta, x) {
        theta["c"] + exp(theta["b0"] + theta["b1"] * x) / 
          (1 + theta["b2"] * exp(theta["b0"]) / 
             theta["b1"] * (exp(theta["b1"] * x) - 1))
      }
    } else {
      mortfun <- function(theta, x) {
        exp(theta["a0"] - theta["a1"] * x) + theta["c"] + 
          exp(theta["b0"] + theta["b1"] * x) / 
          (1 + theta["b2"] * exp(theta["b0"]) / 
             theta["b1"] * (exp(theta["b1"] * x) - 1))
      }
    }
  }
  return(mortfun)
}

# b) General cumulative hazards functions:
.DefineCumHazMatrix <- function(algObj) {
  if (algObj$model == "EX") {
    cumhazfun <- function(theta, x) c(theta) * x
  } else if (algObj$model == "GO") {
    if (algObj$shape == "simple") {
      cumhazfun <- function(theta, x) {
        exp(theta[, "b0"]) / theta[, "b1"] * 
          (exp(theta[, "b1"] * x) - 1)
      }
    } else if (algObj$shape == "Makeham") {
      cumhazfun <- function(theta, x) {
        theta[, "c"] * x + exp(theta[, "b0"]) / theta[, "b1"] * 
          (exp(theta[, "b1"] * x) - 1)
      }
    } else {
      cumhazfun <- function(theta, x) {
        exp(theta[, "a0"]) / theta[, "a1"] * (1 - exp(-theta[, "a1"] * x)) + 
          theta[, "c"] * x + exp(theta[, "b0"]) / theta[, "b1"] * 
          (exp(theta[, "b1"] * x) - 1)
      }
    }
  } else if (algObj$model == "WE") {
    if (algObj$shape == "simple") {
      cumhazfun <- function(theta, x) {
        (theta[, "b1"] * x)^theta[, "b0"]
      }      
    } else if (algObj$shape == "Makeham") {
      cumhazfun <- function(theta, x) {
        theta[, "c"] * x + (theta[, "b1"] * x)^theta[, "b0"]
      }
    } else {
      cumhazfun <- function(theta, x) {
        exp(theta[, "a0"]) / theta[, "a1"] * (1 - exp(-theta[, "a1"] * x)) +
          theta[, "c"] * x + (theta[, "b1"] * x)^theta[, "b0"]
      }
    }
  } else if (algObj$model == "LO") {
    if (algObj$shape == "simple") {
      cumhazfun <- function(theta, x) {
        log(1 + theta[, "b2"] * exp(theta[, "b0"]) / theta[, "b1"] * 
              (exp(theta[, "b1"] * x) - 1)) * (1 / theta[, "b2"])
      }
    } else if (algObj$shape == "Makeham") {
      cumhazfun <- function(theta, x) {
        theta[, "c"] * x + log(1 + theta[, "b2"] * exp(theta[, "b0"]) / 
                                 theta[, "b1"] * 
                                 (exp(theta[, "b1"] * x) - 1)) * 
          (1 / theta[, "b2"])
      }
    } else {
      cumhazfun <- function(theta, x) {
        exp(theta[, "a0"]) / theta[, "a1"] * (1 - exp(-theta[, "a1"] * x)) +
          theta[, "c"] * x + log(1 + theta[, "b2"] * 
                                   exp(theta[, "b0"]) / theta[, "b1"] * 
                                   (exp(theta[, "b1"] * x) - 1)) *
          (1 / theta[, "b2"])
      }
    }
  }
  return(cumhazfun)
}

.DefineCumHazNumeric <- function(algObj) {
  if (algObj$model == "EX") {
    cumhazfun <- function(theta, x) c(theta) * x
  } else if (algObj$model == "GO") {
    if (algObj$shape == "simple") {
      cumhazfun <- function(theta, x) {
        exp(theta["b0"]) / theta["b1"] * 
          (exp(theta["b1"] * x) - 1)
      }
    } else if (algObj$shape == "Makeham") {
      cumhazfun <- function(theta, x) {
        theta["c"] * x + exp(theta["b0"]) / theta["b1"] * 
          (exp(theta["b1"] * x) - 1)
      }
    } else {
      cumhazfun <- function(theta, x) {
        exp(theta["a0"]) / theta["a1"] * (1 - exp(-theta["a1"] * x)) + 
          theta["c"] * x + exp(theta["b0"]) / theta["b1"] * 
          (exp(theta["b1"] * x) - 1)
      }
    }
  } else if (algObj$model == "WE") {
    if (algObj$shape == "simple") {
      cumhazfun <- function(theta, x) {
        (theta["b1"] * x)^theta["b0"]
      }      
    } else if (algObj$shape == "Makeham") {
      cumhazfun <- function(theta, x) {
        theta["c"] * x + (theta["b1"] * x)^theta["b0"]
      }
    } else {
      cumhazfun <- function(theta, x) {
        exp(theta["a0"]) / theta["a1"] * (1 - exp(-theta["a1"] * x)) +
          theta["c"] * x + (theta["b1"] * x)^theta["b0"]
      }
    }
  } else if (algObj$model == "LO") {
    if (algObj$shape == "simple") {
      cumhazfun <- function(theta, x) {
        log(1 + theta["b2"] * exp(theta["b0"]) / theta["b1"] * 
              (exp(theta["b1"] * x) - 1)) * (1 / theta["b2"])
      }
    } else if (algObj$shape == "Makeham") {
      cumhazfun <- function(theta, x) {
        theta["c"] * x + log(1 + theta["b2"] * exp(theta["b0"]) / 
                               theta["b1"] * 
                               (exp(theta["b1"] * x) - 1)) * 
          (1 / theta["b2"])
      }
    } else {
      cumhazfun <- function(theta, x) {
        exp(theta["a0"]) / theta["a1"] * (1 - exp(-theta["a1"] * x)) +
          theta["c"] * x + log(1 + theta["b2"] * 
                                 exp(theta["b0"]) / theta["b1"] * 
                                 (exp(theta["b1"] * x) - 1)) *
          (1 / theta["b2"])
      }
    }
  }
  return(cumhazfun)
}

# c) General survival functions:
.DefineSurvMatrix <- function(algObj) {
  if (algObj$model == "EX") {
    survfun <- function(theta, x) exp(- c(theta) * x)
  } else if (algObj$model == "GO") {
    if (algObj$shape == "simple") {
      survfun <- function(theta, x) {
        exp(exp(theta[, "b0"]) / theta[, "b1"] * 
              (1 - exp(theta[, "b1"] * x)))
      }
    } else if (algObj$shape == "Makeham") {
      survfun <- function(theta, x) {
        exp(-theta[, "c"] * x + exp(theta[, "b0"]) / theta[, "b1"] * 
              (1 - exp(theta[, "b1"] * x)))
      }
    } else {
      survfun <- function(theta, x) {
        exp(exp(theta[, "a0"]) / theta[, "a1"] * (exp(-theta[, "a1"] * x) - 1) - 
              theta[, "c"] * x + exp(theta[, "b0"]) / theta[, "b1"] * 
              (1 - exp(theta[, "b1"] * x)))
      }
    }
  } else if (algObj$model == "WE") {
    if (algObj$shape == "simple") {
      survfun <- function(theta, x) {
        exp(-(theta[, "b1"] * x)^theta[, "b0"])
      }      
    } else if (algObj$shape == "Makeham") {
      survfun <- function(theta, x) {
        exp(-theta[, "c"] * x - (theta[, "b1"] * x)^theta[, "b0"])
      }
    } else {
      survfun <- function(theta, x) {
        exp(exp(theta[, "a0"]) / theta[, "a1"] * (exp(-theta[, "a1"] * x) - 1) -
              theta[, "c"] * x - (theta[, "b1"] * x)^theta[, "b0"])
      }
    }
  } else if (algObj$model == "LO") {
    if (algObj$shape == "simple") {
      survfun <- function(theta, x) {
        (1 + theta[, "b2"] * exp(theta[, "b0"]) / theta[, "b1"] * 
           (exp(theta[, "b1"] * x) - 1))^(-1 / theta[, "b2"])
      }
    } else if (algObj$shape == "Makeham") {
      survfun <- function(theta, x) {
        exp(-theta[, "c"] * x) * (1 + theta[, "b2"] * exp(theta[, "b0"]) / 
                                    theta[, "b1"] * (exp(theta[, "b1"] * x) - 1))^(-1 / theta[, "b2"])
      }
    } else {
      survfun <- function(theta, x) {
        exp(exp(theta[, "a0"]) / theta[, "a1"] * (exp(-theta[, "a1"] * x) - 1) -
              theta[, "c"] * x) * (1 + theta[, "b2"] * 
                                     exp(theta[, "b0"]) / theta[, "b1"] * 
                                     (exp(theta[, "b1"] * x) - 1))^(-1 / theta[, "b2"])
      }
    }
  }
  return(survfun)
}

.DefineSurvNumeric <- function(algObj) {
  if (algObj$model == "EX") {
    survfun <- function(theta, x) exp(- c(theta) * x)
  } else if (algObj$model == "GO") {
    if (algObj$shape == "simple") {
      survfun <- function(theta, x) {
        exp(exp(theta["b0"]) / theta["b1"] * 
              (1 - exp(theta["b1"] * x)))
      }
    } else if (algObj$shape == "Makeham") {
      survfun <- function(theta, x) {
        exp(-theta["c"] * x + exp(theta["b0"]) / theta["b1"] * 
              (1 - exp(theta["b1"] * x)))
      }
    } else {
      survfun <- function(theta, x) {
        exp(exp(theta["a0"]) / theta["a1"] * (exp(-theta["a1"] * x) - 1) - 
              theta["c"] * x + exp(theta["b0"]) / theta["b1"] * 
              (1 - exp(theta["b1"] * x)))
      }
    }
  } else if (algObj$model == "WE") {
    if (algObj$shape == "simple") {
      survfun <- function(theta, x) {
        exp(-(theta["b1"] * x)^theta["b0"])
      }      
    } else if (algObj$shape == "Makeham") {
      survfun <- function(theta, x) {
        exp(-theta["c"] * x - (theta["b1"] * x)^theta["b0"])
      }
    } else {
      survfun <- function(theta, x) {
        exp(exp(theta["a0"]) / theta["a1"] * (exp(-theta["a1"] * x) - 1) -
              theta["c"] * x - (theta["b1"] * x)^theta["b0"])
      }
    }
  } else if (algObj$model == "LO") {
    if (algObj$shape == "simple") {
      survfun <- function(theta, x) {
        (1 + theta["b2"] * exp(theta["b0"]) / theta["b1"] * 
           (exp(theta["b1"] * x) - 1))^(-1 / theta["b2"])
      }
    } else if (algObj$shape == "Makeham") {
      survfun <- function(theta, x) {
        exp(-theta["c"] * x) * (1 + theta["b2"] * exp(theta["b0"]) / 
                                  theta["b1"] * (exp(theta["b1"] * x) - 1))^(-1 / theta["b2"])
      }
    } else {
      survfun <- function(theta, x) {
        exp(exp(theta["a0"]) / theta["a1"] * (exp(-theta["a1"] * x) - 1) -
              theta["c"] * x) * (1 + theta["b2"] * 
                                   exp(theta["b0"]) / theta["b1"] * 
                                   (exp(theta["b1"] * x) - 1))^(-1 / theta["b2"])
      }
    }
  }
  return(survfun)
}

# ------------------------------------------ #
# Likelihood, posterior and prior functions:
# ------------------------------------------ #
# Build initial likelihood object:
.BuildLikeObj <- function(parCovObj, parObj, fullParObj, ageObj,  
                          dataObj) {
  df <- data.frame(mort = rep(0, dataObj$n), ages = rep(0, dataObj$n))
  likeObj <- list(df = df, pars = c())
  likeObj <- .CalcLike(dataObj, parCovObj, parObj, fullParObj, ageObj, likeObj, 
                       ind = 1:dataObj$n)
}

# Likelihood function:
.CalcLike <- function(dataObj, ...) UseMethod(".CalcLike")

.CalcLike.bastacmr <- function(dataObj, parCovObj, parObj, fullParObj, 
                               ageObj, likeObj, ind) {
  # First time created:
  if (inherits(parObj, "lambda")) {
    
    # log-interval censored or survival before min age:
    likeMbef <- log(exp(- parObj$lambda * ageObj$ages$ageBef[ind]) -
                      exp(- parObj$lambda * ageObj$ages$ageBef[ind] + 
                            dataObj$Dx) *
                      (1 - ageObj$inds$ageBef[ind]))
    
    # log-survival truncation before min age:
    likeHbefTr <- parObj$lambda * ageObj$ages$truBef[ind]
  } else {
    likeMbef <- 0
    likeHbefTr <- 0
  }
  # log-mortality after min age:
  likeMaft <- log(.CalcSurv(parCovObj$theta[ind, ], 
                            ageObj$ages$ageAft[ind])^exp(parCovObj$gamma[ind]) -
                    .CalcSurv(parCovObj$theta[ind, ], 
                              ageObj$ages$ageAft[ind] + 
                                dataObj$Dx)^exp(parCovObj$gamma[ind])) * 
    ageObj$inds$uncens[ind] * (1 - ageObj$inds$ageBef[ind])
  
  # log-survival at truncation after min age: 
  # (CHANGE 2021-01-15, Missed prop. hazards)
  likeHaftTr <- exp(parCovObj$gamma[ind]) * 
    .CalcCumHaz(parCovObj$theta[ind, ], ageObj$ages$truAft[ind])
  
  # log-likelihood for ages (without recapture):
  likeAges <- likeMbef + likeMaft
  
  # log-likelihood for mortality parameters:
  likeObj$df$mort[ind] <- likeAges + likeHbefTr + likeHaftTr
  
  # Including recapture probability before first obs. and after last obs.:
  likeRecap <- c((ageObj$alive - dataObj$obsMat)[ind, ] %*% 
                   log(1 - parObj$pi[fullParObj$pi$idpi]))
  
  # log-likelihood for ages:
  likeObj$df$ages[ind] <- likeAges + likeRecap
  
  # Sum of log-like for parameters:
  likeObj$pars <- sum(likeObj$df$mort)
  
  # Assign class to likelihood object:
  class(likeObj) <- "likecmr"
  return(likeObj)
}

.CalcLike.bastacensus <- function(dataObj, parCovObj, parObj, fullParObj, 
                                  ageObj, likeObj, 
                                  ind) {
  # First time created:
  if (inherits(parObj, "lambda")) {
    
    # log-mortality before min age:
    likeMbef <- log(parObj$lambda) * ageObj$inds$ageBef[ind] * 
      ageObj$inds$uncens[ind]
    
    # log-survival before min age:
    likeHbef <- - parObj$lambda * ageObj$ages$ageBef[ind]
    
    # log-survival truncation before min age:
    likeHbefTr <- parObj$lambda * ageObj$ages$truBef[ind]
  } else {
    likeMbef <- 0
    likeHbef <- 0
    likeHbefTr <- 0
  }
  # log-mortality after min age:
  likeMaft <- (parCovObj$gamma[ind] + 
                 log(.CalcMort(parCovObj$theta[ind, ], 
                               ageObj$ages$ageAft[ind]))) * 
    ageObj$inds$uncens[ind] * (1 - ageObj$inds$ageBef[ind])
  
  # log-survival after min age:
  likeHaft <- - exp(parCovObj$gamma[ind]) * 
    .CalcCumHaz(parCovObj$theta[ind, ], ageObj$ages$ageAft[ind])
  
  # log-survival at truncation after min age: 
  # (CHANGE 2021-01-15, Missed prop. hazards)
  likeHaftTr <- exp(parCovObj$gamma[ind]) * 
    .CalcCumHaz(parCovObj$theta[ind, ], ageObj$ages$truAft[ind])
  
  # log-likelihood for ages (without recapture):
  likeAges <- likeMbef + likeHbef + likeMaft + likeHaft
  
  # log-likelihood for mortality parameters:
  likeObj$df$mort[ind] <- likeAges + likeHbefTr + likeHaftTr
  
  # log-likelihood for ages:
  likeObj$df$ages[ind] <- likeAges
  
  # Sum of log-like for parameters:
  likeObj$pars <- sum(likeObj$df$mort)
  
  # Assign class to likelihood object:
  class(likeObj) <- "likecensus"
  return(likeObj)
}


# Create initial posterior object:
.CreatePost <- function(likeObj, fullParObj, parObj, dataObj,
                        ageObj, covObj, priorAgeObj = NULL) {
  # Create posterior list object:
  postObj <- list(priorThe = 0, priorGam = 0, priorLam = 0, postPars = 0, 
                  postAge = rep(0, dataObj$n), priorAge = rep(0, dataObj$n))
  
  # Prior for theta parameters:
  # ---------------------------
  postObj$priorThe <- sum(.dtnorm(c(parObj$theta), 
                                  c(fullParObj$theta$priorMean), 
                                  c(fullParObj$theta$priorSd), 
                                  low = c(fullParObj$theta$low), log = TRUE))
  
  # Prior for gamma parameters:
  # ----------------------------
  if (inherits(parObj, "theGam")) {
    postObj$priorGam <- sum(.dtnorm(parObj$gamma, fullParObj$gamma$priorMean,
                                    fullParObj$gamma$priorSd, log = TRUE))
  } 
  
  # Prior for lambda parameter:
  # ---------------------------
  if (inherits(parObj, "lambda")) {
    postObj$priorLam <- .dtnorm(parObj$lambda, fullParObj$lambda$priorMean,
                                fullParObj$lambda$priorSd, low = 0, 
                                log = TRUE)
  } 
  
  # Posterior for all parameters:
  # ------------------------------
  postObj$postPars <- likeObj$pars + postObj$priorThe + postObj$priorGam + 
    postObj$priorLam
  
  # Posterior for ages:
  # -------------------
  if (inherits(dataObj, "bastacmr")) {
    postObj$priorAge <- .CalcPriorAgeDist(ageObj, priorAgeObj, 
                                          ind = 1:dataObj$n)
  }
  postObj$postAge <- likeObj$df$ages + postObj$priorAge
  
  return(postObj)
}

# Update posterior:
.CalcPost <- function(likeObj, fullParObj, parObj, postObj, dataObj, ageObj,
                      covObj, priorAgeObj = NULL, type = "theta", 
                      ind = NA) {
  if (type == "theta") {
    # Theta parameters:
    # ------------------
    postObj$priorThe <- sum(.dtnorm(parObj$theta, fullParObj$theta$priorMean,
                                    fullParObj$theta$priorSd, 
                                    low = fullParObj$theta$low, log = TRUE))
  } else if (type == "gamma") {
    # Gamma parameters:
    # -----------------
    postObj$priorGam <- sum(.dtnorm(parObj$gamma, fullParObj$gamma$priorMean,
                                    fullParObj$gamma$priorSd, log = TRUE))
  } else if (type == "lambda") {
    # Lambda parameter:
    # -----------------
    postObj$priorLam <- .dtnorm(parObj$lambda, fullParObj$lambda$priorMean,
                                fullParObj$lambda$priorSd, low = 0, 
                                log = TRUE)
  } else {
    # posterior for ages:
    # -------------------
    if (inherits(dataObj, "bastacmr")) {
      postObj$priorAge[ind] <- .CalcPriorAgeDist(ageObj, priorAgeObj, 
                                            ind = ind)
    }
  }
  # Posterior for all parameters:
  # ------------------------------
  postObj$postPars <- likeObj$pars + postObj$priorThe + postObj$priorGam + 
    postObj$priorLam
  postObj$postAge <- likeObj$df$ages + postObj$priorAge
  return(postObj)
}

# Create object to calculate the prior age distribution:
.SetPriorAgeDist <- function(fullParObj, dataObj, covObj) {
  if (inherits(dataObj, "bastacmr")) {
    dxx <- 0.001
    xx <- seq(0,100,dxx)
    parsPrior <- list()
    parsPrior$theta <- fullParObj$theta$priorMean
    if (inherits(fullParObj, "theGam")) {
      parsPrior$gamma <- fullParObj$gamma$priorMean
    }
    class(parsPrior) <- class(fullParObj)
    priorParsCov <- .CalcCovPars(parsPrior, parCov = NULL, covObj, dataObj, 
                                 type = "both")
    ExVec <- sapply(1:dataObj$n, function(pi) {
      the <- priorParsCov$theta[pi, ]
      gam <- priorParsCov$gamma[pi]
      sum(.CalcSurv(the, xx)^exp(gam) * dxx)
    })
    priorAgeObj <- list(lifeExp = ExVec, pars = parsPrior, 
                        parsCov = priorParsCov)
  } else {
    priorAgeObj <- NA
  }
  return(priorAgeObj)
}

# function to calculate prior age distribution for basta cmr:
.CalcPriorAgeDist <- function(ageObj, priorAgeObj, ind) {
  the <- priorAgeObj$parsCov$theta[ind, ]
  gam <- priorAgeObj$parsCov$gamma[ind]
  lifeExp <- priorAgeObj$lifeExp[ind]
  ages <- ageObj$ages$ageAft[ind]
  priorAge <- (.CalcSurv(the, ages)^exp(gam) / lifeExp) *
    (1 - ageObj$inds$ageBef[ind])
  return(priorAge)
}

# Calculate Hastings ratio for Metropolis-Hastings sampling:
.CalcHastRatio <- function(pNow, pNew, parJumps, type = "theta", pp) {
  if (type == 'theta') {
    hRatio <- .dtnorm(pNow$theta[pp], pNew$theta[pp], parJumps$theta[pp], 
                      low = fullParObj$theta$low[pp], log = TRUE) -
      .dtnorm(pNew$theta[pp], pNow$theta[pp], parJumps$theta[pp], 
              low = fullParObj$theta$low[pp], log = TRUE)
  } else if (type == 'gamma') {
    if (inherits(pNow, "theGam")) {
      hRatio <- .dtnorm(pNow$gamma[pp], pNew$gamma[pp], 
                        parJumps$gamma[pp],
                        low = fullParObj$gamma$low[pp], log = TRUE) -
        .dtnorm(pNew$gamma[pp], pNow$gamma[pp], 
                parJumps$gamma[pp],
                low = fullParObj$gamma$low[pp], log = TRUE)
    }
  } else {
    if (type == "lambda") {
      hRatio <- .dtnorm(pNow$lambda, pNew$lambda, parJumps$lambda, 
                        low = fullParObj$lambda$low, log = TRUE) -
        .dtnorm(pNew$lambda, pNow$lambda, parJumps$lambda, 
                low = fullParObj$lambda$low, log = TRUE)
    }
  }
  return(hRatio)
}

# ============================================= #
# ==== SAMPLING, JUMP, AND MCMC FUNCTIONS: ====
# ============================================= #
# -------------------- #
# Sampling parameters:
# -------------------- #
# Jitter parameters for different starting values in the 
# MCMC chains:
.JitterPars <- function(parObj, fullParObj) {
  parObjn <- parObj
  nthe <- fullParObj$theta$len
  parObjn$theta[1:nthe] <- .rtnorm(nthe, c(parObj$theta), 
                                   c(fullParObj$theta$jitter), 
                                   low = c(fullParObj$theta$low))
  if (inherits(parObj,  "theGam")) {
    parObjn$gamma <- .rtnorm(fullParObj$gamma$len,
                             fullParObj$gamma$start, 
                             fullParObj$gamma$jitter)
  }
  if (inherits(parObj,  "lambda")) {
    parObjn$lambda <- .rtnorm(fullParObj$lambda$len, 
                              fullParObj$lambda$start, 
                              fullParObj$lambda$jitter, 
                              lower = fullParObj$lambda$low)
  }
  if (inherits(parObj, "pi")) {
    rho2 <- fullParObj$pi$prior2 + 
      t(t(ageObj$alive - dataObj$Y) %*% rep(1, dataObj$n))
    Rho2 <- tapply(rho2, fullParObj$pi$idpi, sum)
    parObjn$pi <- rbeta(fullParObj$pi$n, fullParObj$pi$Prior1, Rho2)
  }
  return(parObjn)
}

# Sample parameters:
.SamplePars <- function(parObj, parJumps, fullParObj, dataObj,
                        ageObj, type = "theta", pp) {
  parObjn <- parObj
  if (type == 'theta') {
    nthe <- fullParObj$theta$len
    parObjn$theta[pp] <- .rtnorm(1, parObj$theta[pp], parJumps$theta[pp], 
                                 low = fullParObj$theta$low[pp])
  } else if (type == 'gamma') {
    if (inherits(parObj, "theGam")) {
      parObjn$gamma[pp] <- .rtnorm(1, parObj$gamma[pp], parJumps$gamma[pp],
                                   low = fullParObj$gamma$low[pp])
    }
  } else if (type == "lambda") {
    if (inherits(parObj, "lambda")) {
      parObjn$lambda <- .rtnorm(1, parObj$lambda, parJumps$lambda, 
                                low = fullParObj$lambda$low)
    }
  } else {
    if (inherits(parObj, "pi")) {
      rho2 <- fullParObj$pi$prior2 + 
        t(t(ageObj$alive - dataObj$Y) %*% rep(1, dataObj$n))
      Rho2 <- tapply(rho2, fullParObj$pi$idpi, sum)
      parObjn$pi <- rbeta(fullParObj$pi$n, fullParObj$pi$Prior1, Rho2)
    }
  }
  return(parObjn)
}

# -------------- #
# Sampling ages:
# -------------- #
# Sample age object:
.SampleAges <- function(ageObj, ...) UseMethod(".SampleAges")

.SampleAges.agecmr <- function(ageObj, dataObj, algObj) {
  if (dataObj$updA) {
    
    # SAMPLE BIRTH DATES:
    if (dataObj$updB) {
      bi <- ageObj$ages$birth
      
      # Sample births:
      bi[dataObj$idNoB] <- bi[dataObj$idNoB] + 
        sample(-1:1, dataObj$nUpdB, replace = TRUE)
      
      # Find which first obs is not 0:
      idFirst <- dataObj$idNoB[which(dataObj$firstObs[dataObj$idNoB] > 0)]
      # bi[idFirst] <- apply(cbind(bi[idFirst], 
      #                            dataObj$firstObs[idFirst]), 1, min)
      
      if (length(idFirst) > 0) {
        # Find the min between bi and first obs:
        idBiFi <- which(bi[idFirst] > dataObj$firstObs[idFirst] - 1)
        if (length(idBiFi) > 0) {
          bi[idFirst[idBiFi]] <- dataObj$firstObs[idFirst[idBiFi]] - 1
        }
      }
      
      # Find which do not have a first obs:
      idNoFirst <- dataObj$idNoB[which(dataObj$firstObs[dataObj$idNoB] == 0)]
      if (length(idNoFirst) > 0) {
        idBiFi <- which(bi[idNoFirst] > dataObj$ages$death[idNoFirst])
        if (length(idBiFi) > 0) {
          bi[idNoFirst[idBiFi]] <- ageObj$ages$death[idNoFirst[idBiFi]]
        }
        # bi[idNoFirst] <- apply(cbind(bi[idNoFirst], 
        #                              ageObj$ages$death[idNoFirst]), 1, min)
      }
      ageObj$ages$birth <- bi
    }
    
    # SAMPLE DEATH DATES:
    if (dataObj$updD) {
      di <- ageObj$ages$death
      # Sample death year:
      di[dataObj$idNoD] <- di[dataObj$idNoD] +
        sample(-1:1, dataObj$nUpdD, replace = TRUE)
      
      # Find which have a last obs:
      idLast <- dataObj$idNoD[which(dataObj$lastObs[dataObj$idNoD] > 0)]
      if (length(idLast) > 0) {
        # Find which have earlier death than last obs:
        idDiLi <- which(di[idLast] < dataObj$lastObs[idLast])
        if (length(idDiLi) > 0) {
          di[idLast[idDiLi]] <- dataObj$lastObs[idLast[idDiLi]]
        }
      }
      # di[idLast] <- apply(cbind(di[idLast], dataObj$lastObs[idLast]),
      #                     1, max)
      # Find which do no have a last obs:
      idNoLast <- dataObj$idNoD[which(dataObj$lastObs[dataObj$idNoD] == 0)]
      if (length(idNoLast) > 0) {
        # Find which have earlier death than birth:
        idDiLi <- which(di[idNoLast] < ageObj$ages$birth[idNoLast])
        if (length(idDiLi) > 0) {
          di[idNoLast[idDiLi]] <- ageObj$ages$birth[idNoLast[idDiLi]]
        }
      }
      # di[idNoLast] <- apply(cbind(di[idNoLast], ageObj$ages$birth[idNoLast]),
      #                       1, max)
      ageObj$ages$death <- di
    }
    ageObj$ages$age <- ageObj$ages$death - ageObj$ages$birth
    ageObj$ages$ageTr <- algObj$start - ageObj$ages$birth
    ageObj$ages$ageTr[which(ageObj$ages$ageTr < 0)] <- 0
    
    # Create alive matrix:
    firstObs <- ageObj$ages$birth + 1
    firstObs[which(firstObs < algObj$start)] <- algObj$start
    # firstObs <- c(apply(cbind(algObj$start, ageObj$ages$birth + 1), 1, max))
    lastObs <- ageObj$ages$death
    # lastObs[which(lastObs > algObj$end)] <- algObj$end
    idlc <- which(lastObs > dataObj$censTime)
    lastObs[idlc] <- dataObj$censTime[idlc]
    # lastObs <- c(apply(cbind(algObj$end, dataObj$censTime, 
    #                          ageObj$ages$death), 1, min))
    ageObj$alive <- .BuildAliveMatrix(firstObs, lastObs, dataObj)
    
    # Assign ages if minAge > 0:
    ageObj <- .AssignMinAge(ageObj, dataObj, algObj)
  }
  return(ageObj)
}

.SampleAges.agecensus <- function(ageObj, dataObj, algObj) {
  if (dataObj$updB) {
    ageObj$ages$birth[dataObj$idNoB] <-     
      round(.rtnorm(dataObj$nUpdB, ageObj$ages$birth[dataObj$idNoB], 
                    0.2, lower = dataObj$bil[dataObj$idNoB], 
                    upper = dataObj$biu[dataObj$idNoB]), 2)
    ageObj$ages$age <- dataObj$lastObs - ageObj$ages$birth
    ageObj$ages$ageTr <- dataObj$firstObs - ageObj$ages$birth
    ageObj$ages$ageTr[ageObj$ages$ageTr < 0] <- 0
    
    # Assign ages if minAge > 0:
    ageObj <- .AssignMinAge(ageObj, dataObj, algObj)
  }
  
  return(ageObj)
}

# Recalculate ages based on minAge:
.AssignMinAge <- function(ageObj, ...) UseMethod(".AssignMinAge")

.AssignMinAge.minAge <- function(ageObj, dataObj, algObj) {
  # Ages after min age:
  ageObj$ages$ageAft <- ageObj$ages$age - algObj$minAge
  ageObj$ages$ageAft[ageObj$ages$ageAft < 0] <- 0
  
  # Ages and indicator before min age:
  ageObj$ages$ageBef <- ageObj$ages$age
  ageObj$inds$ageBef <- rep(0, dataObj$n)
  ageObj$inds$ageBef[which(ageObj$ages$ageBef < algObj$minAge)] <- 1
  ageObj$ages$ageBef[which(ageObj$ages$age >= algObj$minAge)] <- algObj$minAge
  
  # Ages at truncation after min age:
  ageObj$ages$truAft <- ageObj$ages$ageTr - algObj$minAge
  ageObj$ages$truAft[which(ageObj$ages$truAft < 0)] <- 0
  
  # Ages at truncation before min age:
  ageObj$ages$truBef <- ageObj$ages$ageTr
  ageObj$ages$truBef[which(ageObj$ages$ageTr >= algObj$minAge)] <- algObj$minAge
  
  return(ageObj)
}

.AssignMinAge.noMinAge <- function(ageObj, dataObj, algObj) {
  ageObj$ages$ageAft <- ageObj$ages$age
  ageObj$ages$truAft <- ageObj$ages$ageTr
  return(ageObj)
}

# ------------------------------ #
# Dynamic jump update functions:
# ------------------------------ #
# Update jumps:
.UpdateJumps <- function(parJumps, jumpsMat, iter, iterUpd, updTarg) {
  parnames <- names(jumpsMat)
  for (pn in parnames) {
    if (ncol(jumpsMat[[pn]]) > 1) {
      updRate <- apply(jumpsMat[[pn]][iter - ((iterUpd - 1):0), ], 2, sum) / 
        iterUpd  
    } else {
      updRate <- sum(jumpsMat[[pn]][iter - ((iterUpd - 1):0), ]) / 
        iterUpd 
    }
    updRate[updRate == 0] <- 1e-2
    if (is.matrix(parJumps[[pn]])) {
      parJumps[[pn]] <- parJumps[[pn]] * 
        matrix(updRate, nrow(parJumps[[pn]]), ncol(parJumps[[pn]])) / updTarg
    } else {
      parJumps[[pn]] <- parJumps[[pn]] * updRate / updTarg
    }
  }
  return(parJumps)
}

# Function to create parameter jump sd object:
.SetParJumps <- function(parObj, fullParObj) {
  parJump <- parObj
  nthe <- fullParObj$theta$len
  parJump$theta[1:nthe] <- rep(0.1, nthe)
  if (inherits(parObj, "theGam")) {
    parJump$gamma <- rep(0.1, fullParObj$gamma$len)
  }
  if (inherits(parObj, "lambda")) {
    parJump$lambda <- 0.001
  }
  return(parJump)
}

# --------------------------- #
# Functions for MCMC outputs:
# --------------------------- #
# Function to create output and jumps matrices:
.CreateMCMCoutObj <- function(fullParObj, dataObj, algObj, parObj, ageObj,
                              likeObj, postObj, type = "mcmc", matrows) {
  McmcOutObj <- list()
  # Create matrices for parameters:
  # matrows <- ifelse(type == "mcmc", algObj$niter, algObj$burnin)
  McmcOutObj$theta <- matrix(0, matrows, fullParObj$theta$len,
                             dimnames = list(NULL, fullParObj$theta$names))
  if (inherits(fullParObj, "theGam")) {
    McmcOutObj$gamma <- matrix(0, matrows, fullParObj$gamma$len, 
                               dimnames = list(NULL, fullParObj$gamma$names))
  }
  if (inherits(fullParObj, "lambda")) {
    McmcOutObj$lambda <- matrix(0, matrows, 1)
  }
  if (inherits(fullParObj, "pi") & type == "mcmc") {
    McmcOutObj$pi <- matrix(0, matrows, fullParObj$pi$len)
  }
  McmcOutObj <- .FillMCMCoutObj(McmcOutObj, parObj, dataObj, ageObj, likeObj,
                                postObj, iter = 1, variables = "params", 
                                type = type)
  # Fill up initial values:
  if (type == "mcmc") {
    nkeep <- ceiling((algObj$niter - algObj$burnin) / algObj$thinning)
    # Create matrices for unknown births and deaths:
    if (dataObj$updB) {
      McmcOutObj$B <- matrix(0, nkeep, dataObj$nUpdB)
      # McmcOutObj$B[1, ] <- ageObj$ages[dataObj$idNoB, "birth"]
    } else {
      McmcOutObj$B <- NA
    }
    if (inherits(dataObj, "bastacmr")) {
      if (dataObj$updD) {
        McmcOutObj$D <- matrix(0, nkeep, dataObj$nUpdD)
        # McmcOutObj$D <- ageObj$ages[dataObj$idNoD, "death"]
      } else {
        McmcOutObj$D <- NA
      }
    } else {
      McmcOutObj$D <- NA
    }
    # Create matrix for likelihood and posterior:
    McmcOutObj$likePost <- matrix(NA, nrow = matrows, ncol = 2,
                                  dimnames = list(NULL, c("Like", "Post")))
    McmcOutObj$likePost[1, ] <- c(likeObj$pars, postObj$postPars)
  }
  return(McmcOutObj)
}

# Function to fill in MCMC output list:
.FillMCMCoutObj <- function(McmcOutObj, parObj, dataObj, ageObj, likeObj,
                            postObj, iter, variables = "params", 
                            type = "mcmc") {
  if (variables == "params") {
    McmcOutObj$theta[iter, ] <- c(parObj$theta)
    if (inherits(parObj, "theGam")) {
      McmcOutObj$gamma[iter, ] <- parObj$gamma    
    }
    if (inherits(parObj, "lambda")) {
      McmcOutObj$lambda[iter, ] <- parObj$lambda
    }
    if (inherits(parObj, "pi") & type == "mcmc") {
      McmcOutObj$pi[iter, ] <- parObj$pi
    }
  } else if (variables == "ages") {
    if (dataObj$updB) {
      McmcOutObj$B[iter, ] <- ageObj$ages[dataObj$idNoB, "birth"]
    }
    if (inherits(dataObj, "bastacmr")) {
      if (dataObj$updD) {
        McmcOutObj$D[iter, ] <- ageObj$ages[dataObj$idNoD, "death"]
      }
    }
  } else {
    McmcOutObj$likePost[iter, ] <- c(likeObj$pars, postObj$postPars)
  }
  return(McmcOutObj)
}

# -------------- #
# MCMC function:
# -------------- #
.RunMCMC <- function(sim, UpdJumps = TRUE, parJumps = NA) {
  # Create initial age object:
  ageNow <- ageObj

  # Create intial parameters:
  rm(".Random.seed", envir = .GlobalEnv); runif(1)
  parNow <- .JitterPars(parObj, fullParObj)
  parCovNow <- .CalcCovPars(parObj = parNow, parCov = NULL, 
                            covObj = covObj, dataObj = dataObj, 
                            type = "both")
  likeNow <- .BuildLikeObj(parCovObj = parCovNow, parObj = parNow, 
                           fullParObj = fullParObj, ageObj = ageNow, 
                           dataObj = dataObj)
  while (is.na(likeNow$pars) | likeNow$pars == -Inf) {
    parNow <- .JitterPars(parObj, fullParObj)
    parCovNow <- .CalcCovPars(parObj = parNow, parCov = NULL, 
                              covObj = covObj, dataObj = dataObj, 
                              type = "both")
    likeNow <- .CalcLike(dataObj = dataObj, parCovObj = parCovNow, 
                         parObj = parNow, fullParObj = fullParObj, 
                         ageObj = ageNow, likeObj = NULL, ind = 1:dataObj$n)
  }
  postNow <- .CreatePost(likeObj = likeNow, fullParObj = fullParObj, 
                         parObj = parNow, dataObj = dataObj, ageObj = ageNow, 
                         covObj = covObj, priorAgeObj = priorAgeObj)
  # cat("2323, ")
  if (UpdJumps) {
    # Start jumps for Metropolis algorithm:
    niter <- 7500
    burnin <- niter
    iterUpd <- 50
    updTarg <- 0.25
    updSeq <- seq(1, niter, iterUpd)
    nUpdSeq <- length(updSeq)
    parJumps <- .SetParJumps(parObj, fullParObj)
    jumpsMat <- .CreateMCMCoutObj(fullParObj, dataObj, algObj, 
                                  parObj = parJumps,
                                  ageNow, type = "jumps", matrows = nUpdSeq)
    jumpsUpdMat <- .CreateMCMCoutObj(fullParObj, dataObj, algObj, 
                                     parObj = parJumps,
                                     ageNow, type = "jumps", matrows = niter)
    for (ju in 1:length(jumpsMat)) {
      jumpsUpdMat[[ju]] <- jumpsUpdMat[[ju]] * 0
    }
    
    # Create MCMC output object:
    McmcOutObj <- list()
  } else {
    niter <- algObj$niter
    thinSeq <- seq(burnin, niter, thinning)
    
    # Create MCMC output object:
    McmcOutObj <- .CreateMCMCoutObj(fullParObj, dataObj, algObj, 
                                    parObj = parNow, ageNow, likeNow, 
                                    postNow, type = "mcmc", 
                                    matrows = niter)
    
    # Counter for updated ages:
    indAge <- 0
  }
  # cat("2358, ")
  # Start MCMC:
  for (iter in 2:niter) {
    for (tg in c("theta", 'gamma', "lambda")) {
      if (tg %in% names(fullParObj)) {
        ntg <- fullParObj[[tg]]$len
        for (pp in 1:ntg) {
          parNew <- .SamplePars(parNow, parJumps, fullParObj, dataObj, 
                                ageNow, type = tg, pp = pp)
          parCovNew <- .CalcCovPars(parObj = parNew, parCov = NULL, 
                                    covObj = covObj, dataObj = dataObj, 
                                    type = "both")
          likeNew <- .CalcLike(dataObj, parCovNew, parNew, fullParObj, ageNow, 
                               likeObj = likeNow, ind = 1:dataObj$n)
          postNew <- .CalcPost(likeNew, fullParObj, parNew, postNow, 
                               dataObj, ageNow, covObj, priorAgeObj, 
                               type = tg, ind = NA)
          hRatio <- .CalcHastRatio(parNow, parNew, parJumps, type = tg,
                                   pp)
          postRatio <- exp(postNew$postPars - postNow$postPars + hRatio)
          if (!is.na(postRatio) & postRatio > runif(1)) {
            parNow <- parNew
            parCovNow <- parCovNew
            likeNow <- likeNew
            postNow <- postNew
            if (UpdJumps & iter <= burnin) {
              jumpsUpdMat[[tg]][iter, pp] <- 1
            }
          }
        }
      }
    }
    # if (iter == 2) cat("2390, ")
    # Sample recapture probabilities:
    if (inherits(dataObj, "bastacmr")) {
      parNow <- .SamplePars(parNow, parJumps, fullParObj, dataObj, 
                            ageNow, type = "pi", pp = pp)
      likeNow <- .CalcLike(dataObj, parCovNow, parNow, fullParObj, ageNow, 
                           likeObj = likeNow, ind = 1:dataObj$n)
      postNow <- .CalcPost(likeNow, fullParObj, parNow, postNow, dataObj, 
                           ageNow, covObj, priorAgeObj, type = "ages", 
                           ind = dataObj$idNoA)
    }
    # if (iter == 2) cat("2401, ")
    # Fill up paramters in output object:
    if (!UpdJumps) {
      McmcOutObj <- .FillMCMCoutObj(McmcOutObj, parNow, dataObj, ageNow,
                                    likeNow, postNow, iter = iter,
                                    variables = "params",
                                    type = "mcmc")
    }
    
    # Sample ages:
    if (dataObj$updA) {
      ageNew <- .SampleAges(ageNow, dataObj, algObj)
      likeNew <- .CalcLike(dataObj, parCovNow, parNow, fullParObj, ageNew, 
                           likeNow, ind = dataObj$idNoA)
      # CHANGE 2021-08-10: Find subset of times of birth and death that 
      #                    did change.
      if (inherits(ageNew, "agecmr")) {
        idNew <- dataObj$idNoA[which(ageNew$ages[dataObj$idNoA, "birth"] != 
                                       ageNow$ages[dataObj$idNoA, "birth"] | 
                                       ageNew$ages[dataObj$idNoA, "death"] != 
                                       ageNow$ages[dataObj$idNoA, "death"])]
      } else {
        idNew <- dataObj$idNoA[which(ageNew$ages[dataObj$idNoA, "birth"] != 
                                       ageNow$ages[dataObj$idNoA, "birth"])]
      }
      
      likeNew <- .CalcLike(dataObj, parCovNow, parNow, fullParObj, ageNew, 
                           likeNow, ind = idNew)
      # CHANGE 2021-08-10: Replace ageObj by ageNew:
      postNew <- .CalcPost(likeNew, fullParObj, parNow, postNow, dataObj, 
                           ageNew, covObj, priorAgeObj, type = "age", 
                           ind = idNew)
      postRatio <- exp(postNew$postAge[idNew] - 
                         postNow$postAge[idNew])
      idUpdAges <- idNew[!is.na(postRatio) & 
                           postRatio > runif(length(idNew))]

      if (length(idUpdAges) > 0) {
        ageNow$ages[idUpdAges, ] <- ageNew$ages[idUpdAges, ]
        ageNow$inds[idUpdAges, ] <- ageNew$inds[idUpdAges, ]
        if (inherits(dataObj, "bastacmr")) {
          ageNow$alive[idUpdAges, ] <- ageNew$alive[idUpdAges, ]
        }
        # CHANGE 2021-08-10: Missed updating likelihood and posterior:
        likeNow <- .CalcLike(dataObj, parCovNow, parNow, fullParObj, ageNow, 
                             likeNow, ind = idUpdAges)
        postNow <- .CalcPost(likeNow, fullParObj, parNow, postNow, dataObj, 
                             ageNow, covObj, priorAgeObj, type = "age", 
                             ind = idUpdAges)
      }
      
      # Fill up paramters in output object:
      if (!UpdJumps) {
        if (iter %in% thinSeq) {
          indAge <- indAge + 1
          McmcOutObj <- .FillMCMCoutObj(McmcOutObj, parNow, dataObj, ageNow,
                                        likeNow, postNow, iter = indAge,
                                        variables = "ages", type = "mcmc")
        }
      }
    }
    # if (iter == 2) cat("2439.\n")
    if (!UpdJumps) {
      McmcOutObj <- .FillMCMCoutObj(McmcOutObj, parNow, dataObj, ageNow,
                                    likeNow, postNow, iter = iter,
                                    variables = "likepost", type = "mcmc")
    }
    
    # Dynamic Metropolis to update jumps sd:
    if (UpdJumps) {
      if (iter %in% updSeq) {
        idpar <- which(updSeq == iter)
        parJumps <- .UpdateJumps(parJumps, jumpsUpdMat, iter, iterUpd, updTarg)
        # jumpsMat <- .FillMCMCoutObj(jumpsMat, parJumps, dataObj, ageNow, 
        #                             idpar, variables = "params", 
        #                             type = "jumps")
        for (pj in 1:length(jumpsMat)) {
          parName <- names(jumpsMat)[pj]
          jumpsMat[[parName]][idpar, ] <- c(parJumps[[parName]])
        }
        if (iter == max(updSeq)) {
          for (pj in 1:length(jumpsMat)) {
            parName <- names(jumpsMat)[pj]
            len <- fullParObj[[parName]]$len
            njumps <- floor(nUpdSeq / 2):nUpdSeq
            mat <- jumpsMat[[parName]][njumps, ]
            if (is.matrix(mat)) {
              parJumps[[parName]][1:len] <- apply(mat, 2, mean)
            } else {
              parJumps[[parName]][1:len] <- mean(mat)
            }
          }
        }
      }
    }
  }
  McmcOutObj$jumps <- parJumps
  return(McmcOutObj)
}

# ===================================== #
# ==== FUNCTIONS FOR MCMC OUTPUTS: ====
# ===================================== #
# Extract thinned sequences from multiple runs, calculate coefficients,
# DIC and quantiles for mortality and survival:
.ExtractParalOut <- function(bastaOut, keep, fullParObj, covsNames, nsim, 
                             dataObj, algObj, defTheta, .CalcMort, 
                             .CalcMort.numeric, .CalcMort.matrix, 
                             .CalcSurv, .CalcSurv.matrix, 
                             .CalcSurv.numeric, covObj) {
  parMat <- bastaOut[[1]]$theta[keep, ]
  fullParMat <- bastaOut[[1]]$theta
  parnames <- fullParObj$theta$names
  theMat <- parMat
  likePost <- bastaOut[[1]]$likePost[keep, ]
  if (dataObj$updB) {
    birthMat <- bastaOut[[1]]$B
  }
  if (covsNames$class %in% c("propHaz", "fused")) {
    parMat <- cbind(parMat, bastaOut[[1]]$gamma[keep, ])
    parnames <- c(parnames, fullParObj$gamma$names)
    fullParMat <- cbind(fullParMat, bastaOut[[1]]$gamma)
  }
  if (inherits(fullParObj, "lambda")) {
    parMat <- cbind(parMat, bastaOut[[1]]$lambda[keep])
    parnames <- c(parnames, "lambda")
    fullParMat <- cbind(fullParMat, bastaOut[[1]]$lambda)
  }
  if (inherits(fullParObj, "pi")) {
    parMat <- cbind(parMat, bastaOut[[1]]$pi[keep, ])
    parnames <- c(parnames, fullParObj$pi$names)
    fullParMat <- cbind(fullParMat, bastaOut[[1]]$pi)
  }
  for (i in 2:nsim) {
    if (is.matrix(theMat)) {
      theMat <- rbind(theMat, bastaOut[[i]]$theta[keep, ])
    } else {
      theMat <- c(theMat, bastaOut[[i]]$theta[keep, ])
    }
    if (covsNames$class == "inMort") {
      pmat <- bastaOut[[i]]$theta[keep, ]
      fullPmat <- bastaOut[[i]]$theta
    } else {
      pmat <- cbind(bastaOut[[i]]$theta[keep, ],
                    bastaOut[[i]]$gamma[keep, ])
      fullPmat <- cbind(bastaOut[[i]]$theta, bastaOut[[i]]$gamma)
    }
    if (inherits(fullParObj, "lambda")) {
      pmat <- cbind(pmat, bastaOut[[i]]$lambda[keep])
      fullPmat <- cbind(fullPmat, bastaOut[[i]]$lambda)
    }
    if (inherits(fullParObj, "pi")) {
      pmat <- cbind(pmat, bastaOut[[i]]$pi[keep, ])
      fullPmat <- cbind(fullPmat, bastaOut[[i]]$pi)
    }
    
    if (is.matrix(parMat)) {
      parMat <- rbind(parMat, pmat)
      fullParMat <- rbind(fullParMat, fullPmat)
    } else {
      parMat <- c(parMat, pmat)
      parMat <- matrix(parMat, ncol = 1)
      fullParMat <- c(fullParMat, fullPmat)
    }
    likePost <- rbind(likePost, bastaOut[[i]]$likePost[keep, ])
    if (dataObj$updB) {
      birthMat <- rbind(birthMat, bastaOut[[i]]$B)
    }
  }
  colnames(parMat) <- parnames
  coeffs <- cbind(apply(parMat, 2, mean), apply(parMat, 2, sd), 
                  t(apply(parMat, 2, quantile, c(0.025, 0.975))), 
                  apply(parMat, 2, 
                        function(x) cor(x[-1], x[-length(x)], 
                                        use = "complete.obs")),
                  apply(fullParMat, 2, function(x) 
                    length(which(diff(x) != 0)) / (length(x) -1)))
  colnames(coeffs) <- c("Mean", "StdErr", "Lower95%CI", "Upper95%CI", 
                        "SerAutocorr", "UpdateRate")
  if (nsim > 1) {
    nthin <- length(keep)
    idSims <- rep(1:algObj$nsim, each = nthin)
    Means <- apply(parMat, 2, function(x) 
      tapply(x, idSims, mean))
    Vars <- apply(parMat, 2, function(x) 
      tapply(x, idSims, var))
    meanall <- apply(Means, 2, mean)
    B <- nthin / (algObj$nsim - 1) * apply(t((t(Means) - meanall)^2), 2, sum)
    W <- 1 / algObj$nsim * apply(Vars, 2, sum)
    Varpl <- (nthin - 1) / nthin * W + 1 / nthin * B
    Rhat <- sqrt(Varpl / W)
    Rhat[Varpl==0] <- 1
    conv <- cbind(B, W, Varpl, Rhat)
    rownames(conv) <- colnames(parMat)
    coeffs <- cbind(coeffs, conv[, 'Rhat'])
    colnames(coeffs) <- c(colnames(coeffs)[-ncol(coeffs)], "PotScaleReduc")
    idnconv <- which(conv[, 'Rhat'] > 1.1)
    if (length(idnconv) == 0) {
      # DIC:
      Dave <- mean(- 2 * likePost[, 'Like'])
      pD <- 1/2 * var(-2 * likePost[, 'Like'])
      DIC <- pD + Dave
      Dmode <- Dave - 2 * pD
      k <- fullParObj$theta$len + 
        ifelse("gamma" %in% names(fullParObj), fullParObj$gamma$len, 0) +
        ifelse("lambda" %in% names(fullParObj), fullParObj$lambda$len, 0)
      modSel <- c(Dave, Dmode, pD, k, DIC)
      names(modSel) <- c("D.ave", "D.mode", "pD", "k", "DIC")
      convmessage <- "All parameters converged properly.\n"
    } else {
      modSel <- NA
      convmessage <- "Convergence not reached for some parameters.\n"
    }
  } else {
    modSel <- NA
    convmessage <- "Convergence not calculated due to\ninsuficcient number of simulations.\n"
  }
  cat(convmessage)
  # Kullback-Leibler:
  kulLeib <- .CalcKulbackLeibler(coeffs, covObj, defTheta, fullParObj, 
                                 algObj)
  
  # Mortality, survival and density quantiles:
  bQuan <- dataObj$bi
  if (dataObj$updB) {
    bQuan[dataObj$idNoB] <- c(apply(birthMat, 2, quantile, 0.975))
  }
  maxAge <- max(dataObj$lastObs - bQuan)
  xv <- seq(0, maxAge * 4, length = 1000)
  tempOut <- list(params = parMat, theta = theMat)
  mortQuan <- .CalcDemoFunQuan(tempOut, xv, covsNames, defTheta, 
                               funtype = "mort", .CalcMort, 
                               .CalcMort.numeric, .CalcMort.matrix, 
                               .CalcSurv, .CalcSurv.matrix, 
                               .CalcSurv.numeric)
  survQuan <- .CalcDemoFunQuan(tempOut, xv, covsNames, defTheta, 
                               funtype = "surv", .CalcMort, 
                               .CalcMort.numeric, .CalcMort.matrix, 
                               .CalcSurv, .CalcSurv.matrix, 
                               .CalcSurv.numeric)
  densQuan <- .CalcDemoFunQuan(tempOut, xv, covsNames, defTheta, 
                               funtype = "dens", .CalcMort, 
                               .CalcMort.numeric, .CalcMort.matrix, 
                               .CalcSurv, .CalcSurv.matrix, 
                               .CalcSurv.numeric)
  PSQuan <- .CalcDemoFunQuan(tempOut, xv, covsNames, defTheta, 
                             funtype = "PS", .CalcMort, 
                             .CalcMort.numeric, .CalcMort.matrix, 
                             .CalcSurv, .CalcSurv.matrix, 
                             .CalcSurv.numeric)
  cuts <- list()
  for (nta in names(survQuan)) {
    cuts[[nta]] <- which(survQuan[[nta]][1, ] > 0.05)
  }
  out <- list(params = parMat, theta = theMat, coefficients = coeffs, 
              names = parnames, DIC = modSel, KullbackLeibler = kulLeib, 
              PS = PSQuan, mort = mortQuan, surv = survQuan, dens = densQuan, 
              x = xv, cuts = cuts, convergence = conv, 
              convmessage = convmessage)
  
  return(out)
}

# Calculate Kulback-Leibler discrepancies between parameters:
.CalcKulbackLeibler <- function(coef, covObj, defTheta, fullParObj, algObj,
                                dataObj) {
  if (!is.null(covObj$cat) & 
      !(length(covObj$cat) == 2 & inherits(covObj, "propHaz"))) {
    if (length(covObj$cat) > 1) {
      if (inherits(covObj, c("fused", "inMort"))) {
        parNames <- defTheta$name
        nPar <- defTheta$length
        low <- defTheta$low
        nCat <- length(covObj$cat)
        namesCat <- names(covObj$cat)
      } else {
        parNames <- "gamma"
        nCat <- length(covObj$cat) - 1
        namesCat <- names(covObj$cat)[-1]
        nPar <- 1
        low <- -Inf
      }
      nComb <- (nCat - 1)^2 - ((nCat - 1)^2 - (nCat - 1)) / 2
      covComb1 <- c()
      covComb2 <- c()
      klMat1 <- matrix(0, nPar, nComb, dimnames = list(parNames, NULL))
      klMat2 <- klMat1
      comb <- 0
      for (i in 1:nCat) {
        for (j in 1:nCat) {
          if (i > j) {
            comb <- comb + 1
            covComb1 <- c(covComb1, 
                          sprintf("%s - %s", namesCat[i], namesCat[j]))
            covComb2 <- c(covComb2, 
                          sprintf("%s - %s", namesCat[j], namesCat[i]))
            for (p in 1:nPar) {
              if (inherits(covObj, c("fused", "inMort"))) {
                idP <- sapply(c(i, j), function(ij) 
                  which(rownames(coef) == 
                          sprintf("%s.%s", parNames[p], namesCat[ij])))
              } else {
                idP <- sapply(c(i, j), function(ij) 
                  which(rownames(coef) == namesCat[ij]))
              }
              parRan <- range(sapply(1:2, 
                                     function(pp) .qtnorm(c(0.001, 0.999), 
                                                          coef[idP[pp], 1],
                                                          coef[idP[pp], 2], 
                                                          lower = low[p])))
              parVec <- seq(parRan[1], parRan[2], length = 100)
              dp <- parVec[2] - parVec[1]
              parDens <- sapply(1:2, function(pp) 
                .dtnorm(seq(parRan[1], parRan[2], length = 100), 
                        coef[idP[pp], 1], coef[idP[pp], 2], lower = low[p]))
              klMat1[p, comb] <- sum(parDens[, 1] * 
                                       log(parDens[, 1] / parDens[, 2]) * dp)
              klMat2[p, comb] <- sum(parDens[, 2] * 
                                       log(parDens[, 2] / parDens[, 1]) * dp)
            }
          }
        }
      }
      colnames(klMat1) <- covComb1
      colnames(klMat2) <- covComb2
      qKlMat1 <- (1 + (1 - exp(-2 * klMat1)^(1 / 2))) / 2
      qKlMat2 <- (1 + (1 - exp(-2 * klMat2)^(1 / 2))) / 2
      mqKl <- (qKlMat1 + qKlMat2) / 2
      outList <- list(kl1 = klMat1, kl2 = klMat2, qkl1 = qKlMat1, 
                      qkl2 = qKlMat2, mqKl = mqKl)
    } else {
      outList <- "Not calculated"
    }
  } else {
    outList <- "Not calculated"
  }
  return(outList)
}

# Calculate survival and mortality quantiles:
.CalcDemoFunQuan <- function(out, x, covsNames, defTheta, 
                             funtype = "mort", .CalcMort, 
                             .CalcMort.numeric, .CalcMort.matrix, 
                             .CalcSurv, .CalcSurv.matrix, 
                             .CalcSurv.numeric) {
  covinf <- list(th = list(), ga = list())
  fullm <- list(th = list(), ga = list())
  if (covsNames$class == "inMort") {
    fullm$th <- out$theta
    covinf$th$num <- ifelse(length(covsNames$cat) == 0, 0, 
                            length(covsNames$cat))
    if (is.na(covsNames$cat)[1]) {
      covinf$th$name <- ""
    } else {
      covinf$th$name <- names(covsNames$cat)
    }
    fullm$ga$cat <- 0
    fullm$ga$con <- 0
    covinf$ga$caname <- ""
    covinf$ga$coname <- ""
    covinf$ga$canum <- 0
    covinf$ga$conum <- 0
  } else if (covsNames$class == "fused") {
    fullm$th <- out$theta
    covinf$th$num <- ifelse(length(covsNames$cat) == 0, 0, 
                            length(covsNames$cat))
    if (is.na(covsNames$cat)[1]) {
      covinf$th$name <- ""
    } else {
      covinf$th$name <- names(covsNames$cat)
    }
    covinf$ga$caname <- ""
    covinf$ga$canum <- 0
    covinf$ga$coname <- names(covsNames$con)
    covinf$ga$conum <- length(covsNames$con)
    fullm$ga$cat <- 0
    concol <- colnames(out$params)[(ncol(out$theta)+1):ncol(out$params)]
    idlambda <- grep("lambda", concol)
    if (length(idlambda) > 0) {
      concol <- concol[-idlambda]
    }
    fullm$ga$con <- out$params[, concol]
  } else {
    fullm$th <- out$theta
    covinf$th$num <- 0
    covinf$th$name <- ""
    if (is.na(covsNames$cat)[1]) {
      covinf$ga$caname <- ""
      covinf$ga$canum <- 0
      fullm$ga$cat <- 0
    } else {
      covinf$ga$caname <- names(covsNames$cat)
      covinf$ga$canum <- length(covsNames$cat)
      concol <- colnames(out$params)[(ncol(out$theta)+1):ncol(out$params)]
      idlambda <- grep("lambda", concol)
      if (length(idlambda) > 0) {
        concol <- concol[-idlambda]
      }
      fullm$ga$cat <- cbind(0, out$params[, covinf$ga$caname[-1]])
      colnames(fullm$ga$cat) <- covinf$ga$caname
    }
    if (is.na(covsNames$con)[1]) {
      covinf$ga$coname <- ""
      covinf$ga$conum <- 0
      fullm$ga$con <- 0
    } else {
      covinf$ga$coname <- names(covsNames$con)
      covinf$ga$conum <- length(covsNames$con)
      fullm$ga$con <- out$params[, covinf$ga$coname]
    }
  }
  idlambda <- grep("lambda", colnames(out$params))
  # if (length(idlambda) == 1) {
  #   etav <- out$params[, idlambda]
  # } else {
  #   etav <- rep(0, nrow(out$params))
  # }
  demoQuan <- list()
  for (nta in covinf$th$name) {
    if (nta == "") {
      ntal <- "nocov"
      catname <- "nocov"
    } else {
      ntal <- nta
      catname <- nta
    }
    #demoQuan[[ntal]] <- list()
    for (nga in covinf$ga$caname) {
      if (nga == "") {
        ngal <- "nocov"
      } else {
        ngal <- nga
        catname <- nga
      }
      if (is.matrix(fullm$th)) {
        thm <- fullm$th[, grep(nta, colnames(fullm$th))]
      } else {
        thm <- matrix(fullm$th, ncol = 1)
      }
      if (!is.matrix(thm)) thm <- matrix(thm, ncol = 1)
      colnames(thm) <- defTheta$name
      if (covinf$ga$canum <= 1) {
        gaca <- fullm$ga$cat
      } else {
        gaca <- fullm$ga$cat[, grep(nga, colnames(fullm$ga$cat))]
      }
      if (covinf$ga$conum == 0) {
        gaco <- fullm$ga$con
      } else {
        if (length(covObj$con) > 1) {
          gaco <- sum(fullm$ga$con * apply(covObj$propHaz[, names(covObj$cont)],
                                           2, mean))
        } else {
          gaco <- fullm$ga$con * mean(covObj$propHaz[, names(covObj$cont)])
        }
      }
      gam <- gaca + gaco
      if (length(gam) == 1) {
        gam <- rep(gam, nrow(thm))
      }
      if (funtype %in% c("mort", "surv", "dens")) {
        demof <- sapply(1:nrow(thm), function(ii) {
          if (funtype == "mort") {
            demf <- .CalcMort(thm[ii, ], x) * exp(gam[ii])
          } else if (funtype == "surv") {
            demf <- .CalcSurv(thm[ii, ], x)^{exp(gam[ii])}
          } else if (funtype == "dens") {
            demf <- .CalcMort(thm[ii, ], x) * exp(gam[ii]) * 
              .CalcSurv(thm[ii, ], x)^{exp(gam[ii])}
          }
          return(demf)
        })
        demofave <- apply(demof, 1, mean)
        demofci <- apply(demof, 1, quantile, c(0.025, 0.975), na.rm = TRUE)
        demoffin <- rbind(demofave, demofci)
        rownames(demoffin) <- c("Mean", "2.5%", "97.5%")
        demoQuan[[catname]] <- demoffin
      } else {
        Deltax <- x[2] - x[1]
        surv <- sapply(1:nrow(thm), function(ii) {
          sx <- .CalcSurv(thm[ii, ], x)^{exp(gam[ii])}
          return(sx)
        })
        Ex <- apply(surv, 2, .CalcEx, dx = Deltax)
        Hx <- apply(surv, 2, .CalcHx, dx = Deltax)
        Epx <- - log(Hx)
        Gx <- apply(surv, 2, .CalcGx, dx = Deltax)
        PSq <- rbind(c(mean(Ex), quantile(Ex, c(0.025, 0.975))),
                     c(mean(Hx), quantile(Hx, c(0.025, 0.975))),
                     c(mean(Epx), quantile(Epx, c(0.025, 0.975), na.rm = T)),
                     c(mean(Gx), quantile(Gx, c(0.025, 0.975))))
        dimnames(PSq) <- list(c("LifeExp", "LifeTableEntropy", "LifespanEqual",
                                "Gini"), 
                              c("Mean", "2.5%", "97.5%"))
        demoQuan[[catname]] <- list(PS = PSq, Ex = Ex, Hx = Hx, Epx = Epx,
                                    Gx = Gx)
      }
    }
  }
  return(demoQuan)
}

# ================================ #
# ==== DEMOGRAPHIC FUNCTIONS: ====
# ================================ #
# ------------------------------------------- #
# Functions to calculate pace-shape measures:
# ------------------------------------------- #
# life expectancy:
.CalcEx <- function(Sx, dx) sum(Sx * dx) / Sx[1]

# Keyfitz's entropy:
.CalcHx <- function(Sx, dx) {
  Sx1 <- Sx[Sx > 0]; Sx1 <- Sx1 / Sx1[1]
  -sum(Sx1 * log(Sx1) * dx) / sum(Sx1 * dx)
}

# Gini coefficient:
.CalcGx <- function(Sx, dx) {
  Sx <- Sx / Sx[1]
  Sx <- Sx[Sx > 0]
  return(1 - 1 / sum(Sx * dx) * sum(Sx^2 * dx))
}

# Coefficient of variation:
.CalcCVx <- function(x, Sx, dx) {
  Sx <- Sx / Sx[1]
  idd <- which(Sx > 0)
  Sx <- Sx[idd]
  x <- (x - x[1])[idd]
  dS <- -diff(Sx)
  dS <- dS / sum(dS)
  ex <- sum(Sx * dx)
  return(sqrt(sum((x[-length(x)] + dx/2 - ex)^2 * dS)) / ex)
}

# Function to calculate final pace-shape measures 
# (REQUIRES CHECKING THE STRUCTURE OF THE COVARIATES OUTPUTS):
.CalcPS <- function(object) {
  xv <- object$x
  dx <- xv[2] - xv[1]
  # NOTE: CHECK HOW TO CONSTRUCT THE COVARIATE OUTPUTS....
  PSqMat <- lapply(names(object$covs$cat), function(ss) {
    idcov <- grep(ss, colnames(object$theta))
    thetaMat <- object$theta[, idcov]
    colnames(thetaMat) <- defTheta$name
    surv <- apply(thetaMat, 1, function(th) .CalcSurv(th, xv))
    Ex <- apply(surv, 2, .CalcEx, dx = dx)
    Hx <- apply(surv, 2, .CalcHx, dx = dx)
    Epx <- - log(Hx)
    Gx <- apply(surv, 2, .CalcGx, dx = dx)
    PSq <- rbind(c(mean(Ex), quantile(Ex, c(0.025, 0.975))),
                 c(mean(Hx), quantile(Hx, c(0.025, 0.975))),
                 c(mean(Epx), quantile(Epx, c(0.025, 0.975), na.rm = T)),
                 c(mean(Gx), quantile(Gx, c(0.025, 0.975))))
    dimnames(PSq) <- list(c("LifeExp", "LifeTableEntropy", "LifespanEqual",
                            "Gini"), 
                          c("Mean", "2.5%", "97.5%"))
    return(PSq)
  })
  names(PSqMat) <- names(object$cov$cat)
  return(PSqMat)
}

# --------------------- #
# Construct life table:
# --------------------- #
# Function to calculate lifetable output:
.CalcLifeTable <- function(bastaOut, lifeTable, covObj, algObj, dataObj) {
  cat("Constructing life table... ")
  if (dataObj$updA) {
    nthin <- ceiling((algObj$niter - algObj$burnin + 1) / algObj$thinning)
    if (dataObj$updB) {
      bMat <- matrix(dataObj$bi, dataObj$n, nthin * algObj$nsim)
    }
    if (dataObj$updD) {
      dMat <- matrix(dataObj$di, dataObj$n, nthin * algObj$nsim)
    } else {
      if (algObj$dataType == "census") {
        dMat <- dataObj$lastObs
      } else {
        dMat <- dataObj$di
      }
    }
    for (sim in 1:algObj$nsim) {
      if (dataObj$updB) {
        bMat[dataObj$idNoB, 1:nthin + (sim - 1) * nthin] <- 
          t(bastaOut[[sim]]$B)
      }
      if (dataObj$updD) {
        dMat[dataObj$idNoD, 1:nthin + (sim - 1) * nthin] <- 
          t(bastaOut[[sim]]$D)
      }
    }
  } else {
    bMat <- dataObj$bi
    if (algObj$dataType == "census") {
      dMat <- dataObj$lastObs
    } else {
      dMat <- dataObj$di
    }
  }
  
  # Age at death or last obs:
  ageLast <- dMat - bMat
  
  # Age first:
  ageFirst <- dataObj$firstObs - bMat
  
  # Quantiles for last ages:
  qAgeLast <- cbind(Mean = apply(ageLast, 1, mean),
                    Median = apply(ageLast, 1, quantile, 0.5),
                    Lower = apply(ageLast, 1, quantile, 0.025),
                    Upper = apply(ageLast, 1, quantile, 0.975))
  
  # Quantiles for first ages:
  qAgeFirst <- cbind(Mean = apply(ageFirst, 1, mean),
                     Median = apply(ageFirst, 1, quantile, 0.5),
                     Lower = apply(ageFirst, 1, quantile, 0.025),
                     Upper = apply(ageFirst, 1, quantile, 0.975))
  qAgeFirst[qAgeFirst < 0] <- 0
  
  # Number of levels:
  nq <- ncol(qAgeFirst)
  
  # level names:
  qnames <- colnames(qAgeFirst)
  
  # Indicator for censored or death:
  departType <- rep("D", dataObj$n)
  if (algObj$dataType == "census") {
    departType[dataObj$idCens] <- "C"
  }
  
  # Find maximum age:
  maxAge <- max(qAgeLast)
  ageFact <- 1
  
  # Verify that maxAge is larger than 3, otherwise change the scale to months:
  # if (maxAge < 3) {
  #   ageFact <- 12
  #   maxAge <- maxAge * ageFact
  #   ageFirst <- ageFirst * ageFact
  #   ageLast <- ageLast * ageFact
  # } else {
  #   ageFact <- 1
  # }
  
  # Age vector:
  agev <- (algObj$minAge * ageFact):maxAge
  nage <- length(agev)
  
  # Create Life tables:
  LT  <- list()
  if (is.null(covObj$cat)) {
    covNames <- c("noCov")
  } else {
    covNames <- names(covObj$cat)
  }
  for (covar in covNames) {
    if (covar == "noCov") {
      idx <- 1:dataObj$n
    } else {
      if (inherits(covObj, c("fused", "inMort"))) {
        covcat <- "inMort"
      } else {
        covcat <- "propHaz"
      }
      if (covcat == "propHaz" & !covar %in% colnames(covObj[[covcat]])) {
        if (length(covNames) > 2) {
          idx <- which(apply(covObj[[covcat]][, covNames[-1]], 1, sum) == 0)
        } else {
          idx <- which(covObj[[covcat]][, covNames[-1]] == 0)
        }
      } else {
        idx <- which(covObj[[covcat]][, covar] == 1)
      }
    }
    LT[[covar]] <- list()
    for (ltt in 1:nq) {
      x <- qAgeLast[idx, ltt]
      xt <- qAgeFirst[idx, ltt]
      depType <- departType[idx]
      
      # Outputs:
      Nx <- Dx <- ax <- rep(0, nage)
      for (xx in 1:nage) {
        # A) EXPOSURES:
        # Find how many entered the interval (including truncated):
        idNx <- which(xt < agev[xx] + 1 & x >= agev[xx])
        
        # Extract ages and departType:
        xf <- xt[idNx]
        xl <- x[idNx]
        dt <- depType[idNx]
        
        # proportion of truncation in interval:
        trp <- xf - agev[xx]
        trp[trp < 0] <- 0
        
        # proportion of censoring:
        cep <- agev[xx] + 1 - xl
        cep[cep < 0] <- 0
        cep[dt == "D"] <- 0
        
        # Calculate exposures:
        nexp <- 1 - trp - cep
        Nx[xx] <- sum(nexp)
        
        # B) DEATHS:
        # Calculate total deaths in the interval:
        idDx <- which(dt == "D" & xl < agev[xx] + 1)
        Dx[xx] <- sum(nexp[idDx])
        
        # C) PROPORTION LIVED BY THOSE THAT DIED IN INTERVAL:
        if (Dx[xx] > 1) {
          ylived <- xl[idDx] - agev[xx]
          ax[xx] <- sum(ylived) / length(idDx)
        } else {
          ax[xx] <- 0
        }
      }
      # Age-specific mortality probability:
      qx <- Dx / Nx
      
      # Age-specific survival probability:
      px <- 1 - qx
      
      # Survivorship (or cumulative survival):
      lx <- c(1, cumprod(px))[1:nage]
      
      # Number of individual years lived within the interval:
      Lx <- lx * (1 + ax * qx)
      Lx[is.na(Lx)] <- 0
      
      # Total number of individual years lived after age x:
      Tx <- rev(cumsum(rev(Lx)))
      
      # Remaining life expectancy after age x:
      ex <- Tx / lx
      ex[is.na(ex)] <- 0
      
      # Life-table:
      LT[[covar]][[qnames[ltt]]] <- 
        data.frame(Ages = agev, Nx = Nx, Dx = Dx, lx = lx, px = px,
                   qx = qx, Lx = Lx, Tx = Tx, ex = ex)      
    }
  }
  cat("done.\n")
  return(LT)
}


