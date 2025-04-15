# consolidated version from prep_brainchart_test.R
library(tidyverse)

# Lifespan is a copy of the braincharts github repo

getCurrentFileLocation <-  function()
{
  this_file <- commandArgs() %>% 
    tibble::enframe(name = NULL) %>%
    tidyr::separate(col=value, into=c("key", "value"), sep="=", fill='right') %>%
    dplyr::filter(key == "--file") %>%
    dplyr::pull(value)
  if (length(this_file)==0)
  {
    this_file <- rstudioapi::getSourceEditorContext()$path
  }
  return(dirname(this_file))
}

fileLocation <- getCurrentFileLocation()

wd <- setwd(file.path(fileLocation, "Lifespan"))
source("100.common-variables.r")
source("101.common-functions.r")
source("300.variables.r")
source("301.functions.r")
setwd(wd)


calibrateBrainCharts <- function(noveldata, phenotype=c("GMV", "CT", "TCV", "SA", "sGMV", "WMV", "Ventricles"), phenofit = NULL, expandedFunc = Calc.Expanded) {
  phenotype <- match.arg(phenotype)
  # the fit will be loaded from file. The argument option
  # is for the bootstrap case
  if (is.null(phenofit)) {
    phenofit <- readRDS(paste0("./LifespanDist/Share/FIT_", phenotype, ".rds"))
  }
  # Mess with tidyverse style variables etc
  #message(phenotype)
  PHENO <- rlang::sym(phenotype)
  # GMV has 5 6 freesurfer levels - the others have 5
  # some phenotypes have special naming. (totalSA, meanCT ...)
  # special case for CT
  # A bit of trial and error to deduce the combination of freesurfer versions
  # used in each model creation. The seleciton below
  # seems to produce pretty good results across the board!
  if (phenotype == "CT") {
    fsfactor <- factor(c("FS6_T1", "FS6_T1T2", "FSInfant", "FS53", "Custom_T1T2")) # suspect no infant CT measures???
    PHENOTRANSFORMED <- rlang::sym("meanCT2Transformed")
    ggg <- mutate(noveldata, meanCT2 = CT)
    
  } else if (phenotype == "SA") {
    fsfactor <- factor(c("FS6_T1", "FS6_T1T2", "FSInfant", "FS53", "Custom_T1T2")) # suspect no infant SA measures???
    PHENOTRANSFORMED <- rlang::sym("totalSA2Transformed")
    ggg <- mutate(noveldata, totalSA2 = SA)
    
  } else if (phenotype %in% c("WMV", "GMV")) {
    #  the full freesurfer list
    fsfactor <- factor(c("FS6_T1", "FS6_T1T2", "Custom", "FS53", "Custom_T1T2", "FSInfant")) # suspect no infant CT measures???
    PHENOTRANSFORMED <- rlang::sym(paste0(phenotype, "Transformed"))
    ggg <- noveldata
  } else if (phenotype == "Ventricles") {
    fsfactor <- factor(c("FS6_T1", "FS6_T1T2", "FSInfant", "FS53", "Custom_T1T2"))
    PHENOTRANSFORMED <- rlang::sym(paste0(phenotype, "Transformed"))
    ggg <- noveldata
  } else if (phenotype %in% c("sGMV", "TCV")) {
    # seems to perform best without FSInfant
    fsfactor <- factor(c("FS6_T1", "FS6_T1T2", "Custom", "FS53", "Custom_T1T2"))
    PHENOTRANSFORMED <- rlang::sym(paste0(phenotype, "Transformed"))
    ggg <- noveldata  
  }
  ggg <- mutate(ggg, !!PHENOTRANSFORMED := !!PHENO/10000, AgeTransformed = log(age_days), Study=study, 
                sex = factor(sex, levels=c("Female", "Male")))
  
  # No tibbles inside braincharts code!
  # difference in default drop behaviour causes problems
  # need to add INDEX.OB, INDEX.ID, INDEX.TYPE
  ggg <- mutate(ggg, 
                INDEX.OB = 1,
                INDEX.TYPE = factor(dx == "CN", levels=c(TRUE,FALSE),labels=c("CN","notCN")),
                INDEX.ID = factor( paste( study, participant, sex, sep="|" ) )
  )
  
  ggg <- mutate(ggg,
                study = factor(study),
                fs_version = factor(fs_version,
                                    levels=levels(fsfactor)))
  
  ggg <- as.data.frame(ggg)
  novel <- list()
  novel$DATA <- as.data.frame(ggg)
  
  ## this is a sample dataset for matching up types etc
  # Needs to be matched to the specific phenotype
  example_calibrate <- read_csv(file.path("Samples/", paste0(phenotype, ".csv")), show_col_types = FALSE)
  example_calibrate <- mutate(example_calibrate, 
                              fs_version="FS53",
                              fs_version = factor(fs_version, 
                              levels=levels(fsfactor)),
                              study = factor(study, levels = names(phenofit$param$mu$ranef)),
                              sex = factor(sex, levels=c("Female", "Male")))
  contrasts(example_calibrate$fs_version) <- contr.sum(length(levels(example_calibrate$fs_version)))
  
  novel$DATA <- ValidateCleanInput(novel$DATA, 
                                   as.data.frame(example_calibrate), 
                                   attr(phenofit$param, "model"), 
                                   phenofit$param)
  
  # This will give a warning about fs_version. I don't think it is a problem. There
  # is a contrasts entry in the fit parameters that causes the contrasts of fs_version 
  # to be set to contr.sum. The warning pops up when the formula without an fs_version
  # term (like sigma) are processed.
  novel$DATA.PRED <- Apply.Param(NEWData=novel$DATA,
                                 FITParam=phenofit$param,
                                 Reference.Holder=NULL, # Done above - passed to ValidateCleanInput
                                 Pred.Set=NULL, Prefix="", Add.Moments=FALSE, Add.Normalise=FALSE, Add.Derivative=FALSE, 
                                 MissingToZero=TRUE,
                                 verbose=FALSE )
  summary(novel$DATA.PRED) ## see that mu.wre, sigma.wre are NA, but nu.wre is not (as there are no missing random-effects)
  ## this is different to the tutorial because sigma does have a random effect, but the tutorial version doesn't
  
  attr(novel$DATA.PRED,"missing.levels")
  
  novel$SUBSET <- novel$DATA.PRED[attr(novel$DATA.PRED,"logical.selectors")$REFIT.VALID,]
  
  # print("novel$SUBSET")
  # print(novel$SUBSET)
  # print("Cur.Param")
  # print(phenofit$param)
  # Calc.Expanded uses optim and Ranef.MLE.Func
  
  EXPANDED <- do.call(expandedFunc, list(NewData=novel$SUBSET,
                            Cur.Param=phenofit$param,
                            Missing=attr(novel$DATA.PRED,"missing.levels"))
                      )
  
  
  novel$DATA.PRED2 <- Apply.Param(NEWData=novel$DATA,
                                  Reference.Holder=NULL,
                                  FITParam=EXPANDED,
                                  Pred.Set=c("l025"=0.025,"l250"=0.250,"m500"=0.5,"u750"=0.750,"u975"=0.975),
                                  Prefix="",
                                  Add.Moments=FALSE, ## does not make sense for observations
                                  Add.Normalise=TRUE,
                                  Add.Derivative=FALSE,  ## does not make sense for observations
                                  MissingToZero=TRUE, NAToZero=TRUE,
                                  verbose=FALSE )
  novel$expanded <- EXPANDED
  return(novel)
}

calibrateBrainChartsIDQuantilePenalty <- function(noveldata, phenotype=c("GMV", "CT", "TCV", "SA", "sGMV", "WMV", "Ventricles"), phenofit=NULL) {
  return(calibrateBrainCharts(noveldata, phenotype = phenotype, phenofit = phenofit, expandedFunc = Calc.Expanded.ID.QuantilePenalty))
}

Calc.Expanded.ID.QuantilePenalty <- function( NewData, Cur.Param, Missing, Prefix="" ) {
  OPT <- optim(par=Missing$Vector, fn=Ranef.MLE.Func.ID.QuantilePenalty,
               Param=Cur.Param, Missing=Missing$Levels, Novel=NewData, Prefix=Prefix,
               method=if(length(Missing$Vector)==1){"Brent"}else{"Nelder-Mead"},
               lower=if(length(Missing$Vector)==1){-1000}else{-Inf},
               upper=if(length(Missing$Vector)==1){ 1000}else{ Inf},
               control=list(maxit=10000000))
  
  if( OPT$convergence > 0 ) {
    return(NULL)
  } else {
    ##
    ## Append new levels to fit object
    ##
    EXPANDED <- Add.New.Ranefs( new.vector=OPT$par, Fit.Extract=Cur.Param, Missing=Missing$Levels )
    EXPANDED$LL <- Ranef.MLE.Func.ID.QuantilePenalty(OPT$par, Param=Cur.Param, Missing=Missing$Levels, Novel=NewData, Prefix=Prefix, Return="LL")
    return(EXPANDED)
  }
}


##
## MLE refitting functions supporting ID as random effect
##
Ranef.MLE.Func.ID.QuantilePenalty <- function( theta, Param, Missing, Novel, Prefix="", Return="optim" ) {
  
  if( missing(theta)|missing(Param)|missing(Missing)|missing(Novel) ) {stop("Mandatory argmument(s) missing")}
  #print("theta")
  #print(theta)
  #print("Missing")
  #print(Missing)
  #print(sapply(Missing,lengths))
  if(length(theta)!=sum(sapply(Missing,lengths))) {
    
    print("sum(sapply(Missing,lengths))")
    print(sum(sapply(Missing,lengths)))
    print(Missing)
    print(theta)
    
  }  
  if(length(theta)!=sum(sapply(Missing,lengths))){stop("Problem: length(theta) != number of Missing levels")}
  
  LL.ranef <- list()
  #print("length(Missing)")
  #print(length(Missing))
  #print(Missing)
  DIST <- get( Param$family )()
  
  for( lIDX in 1:length(Missing) ) {
    #print("lIDX")
    #print(lIDX)
    LAB <- names(Missing)[lIDX]
    #if( length(Missing[[lIDX]])> 1 ) { stop("FATAL ERROR: MLE function currently assumes a single random-effect per gamlss parameter") }
    JDX <- attr(Missing[[lIDX]][[1]],"index") ## NOTE: making assumption of a single random-effect within each gamlss parameter
    
    LL.ranef[[LAB]] <- dnorm(x=theta[JDX],mean=0,sd=Param[[LAB]]$sigma,log=TRUE) ## this can be a vector of length>1
    
    lMATCH <- match( Novel[ , names(Missing[[lIDX]])[1] ], Missing[[lIDX]][[1]] )
    
    if( any(is.na(lMATCH)) ) { stop("This should not happen") }
    Novel[,sprintf("%s%s.wre",Prefix,LAB)]  <-  Novel[,sprintf("%s%s.pop",Prefix,LAB),drop=TRUE] + theta[JDX[lMATCH]]
    
  }
  #print(theta)
  lARGS <- list()
  CHECK <- DIST$parameters
  for( LAB in names(DIST$parameters) ) {
    lARGS[[LAB]] <- DIST[[sprintf("%s.linkinv",LAB)]]( Novel[,sprintf("%s.wre",LAB)] )
    CHECK[[LAB]] <- DIST[[sprintf("%s.valid",LAB)]]( lARGS[[LAB]] )
  }
  if( !all(unlist(CHECK)) ) {stop("Failed distribution parameter checks")}
  #print(c(lARGS,list(x=Novel[,attr(Param,"model")$covariates$Y])))
  # LL.out
  LL.out <- do.call( what=get(paste0("d",Param$family)), args=c(lARGS,list(x=Novel[,attr(Param,"model")$covariates$Y],log=TRUE)))
  #print(LL.out)
  # get the quantiles of the data in the adjusted distributions
  Novel[,"Quantiles"] <- do.call( what=get(paste0("p",Param$family)), args=c(lARGS,list(q=Novel[,attr(Param,"model")$covariates$Y],log=FALSE)))
  #print(c(lARGS,list(q=Novel[,attr(Param,"model")$covariates$Y])))
  UniqueParticipants <- unique(Novel$participant)
  LL.QuantilePenalty <- list()
  for (CurParticipant in UniqueParticipants) {
    Y <- Novel[Novel$participant == CurParticipant, "Quantiles"]
    if(length(Y) > 1) {
      #LL.QuantilePenalty[CurParticipant] <- sum(abs(Y - Y[1]))
      LL.QuantilePenalty[CurParticipant] <- dnorm(x=sum(abs(Y - Y[1])),mean=0,sd=0.4 * (length(Y) - 1),log=TRUE) ## this can be a vector of length>1  
      # browser()
      #print(Y)
      # print("Abs")
      # print(LL.QuantilePenaltyAbs[CurParticipant])
      # print("DNorm")
      # print(LL.QuantilePenalty[CurParticipant])
    }
  }
  StudyTable <- table(Novel$study)
  
  S.LL.QuantilePenalty <- sum(unlist(LL.QuantilePenalty))
  S.LL.QuantilePenalty <- S.LL.QuantilePenalty * max(StudyTable) / min(StudyTable) / 4
  
  # NEWData[out.WHICH[[lTYPE]],sprintf("%s%s.q.%s",Prefix,OUT.NAME,lTYPE)] <- do.call(what=get(paste0("p",FAMILY$family[1])),
  #                                                                                   args=c(out.SHORT[[lTYPE]],
  #                                                                                          list(q=NEWData[out.WHICH[[lTYPE]],OUT.NAME,drop=TRUE])))
  
  #print(c(sum(LL.out), sum(unlist(LL.ranef)), S.LL.QuantilePenalty))
  if( Return=="optim" ) {
    #print(-1 * ( sum(LL.out) + sum(unlist(LL.ranef)) ))
    -1 * ( sum(LL.out) + sum(unlist(LL.ranef)) + S.LL.QuantilePenalty)
  } else if ( Return=="LL" ) {
    sum(LL.out)
  } else if ( Return=="LL+RE" ) {
    ( sum(LL.out) + sum(unlist(LL.ranef)) )
  }
}
