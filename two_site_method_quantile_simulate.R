source("calibrate_braincharts_centiles.R")
source("study_simulations.R")
library('pracma')
library('parallel')
library('furrr')

library('lme4')
library("future")
library("optparse")


## generates data for 22 "studies" across age range
DD <- genSimData()

## select one study with multiple timepoints
simdata <- filter(DD, Study == "V")

# set up data frame with required columns
simdata <- mutate(simdata, 
                  CT = Wand, # set cortical thickness "CT" to Wand, which is the most CT-like data
                  sex=factor(Grp, labels = c("Male", "Female"), levels = c("M", "F")),
                  dx = as.character(Type),
                  age_days = Time * 365 + 280,
                  study = Study, 
                  participant = as.character(ID),
                  fs_version = "Custom_T1T2"
                  )

nshared <- opt$nshared
nbig <- opt$nbig
nsmall <- opt$nsmall

# select two "sites", timepoints 0 and 4
simdata04 <- filter(simdata, obs == 0 | obs == 4)
simdata04 <- mutate(simdata04, study = paste0("V", INDEX.OB))

FORCE <- TRUE


#if(!file.exists(paste0(opt$out, "base04.rds")) || FORCE) {

# perform two-site calibration using the whole sample

   base04 <- calibrateBrainCharts(simdata04, phenotype = "CT")
   saveRDS(base04, paste0(opt$out, "base04.rds"))
#}


wrapCalibrateBrainCharts <- function(simdata0, phenotype, phenofit) {
  # this returns the estimate of the site effect
  if (is.null(phenofit$param)) {
    return(c(sigranef = NA, muranef = NA))
  }
  res <- calibrateBrainCharts(simdata0, phenotype, phenofit)
  ranefs <- c(sigranef = res$expanded$sigma$ranef[["V"]], muranef = res$expanded$mu$ranef[["V"]])
  return(ranefs)
}

wrapCalibrateBrainChartsTwoTimepoint <- function(simdata0, phenotype) {
  # this returns the estimate of the site effect

  res <- calibrateBrainChartsTwoTimepoint(simdata0, phenotype)
  ranefs <- c(sigranef1 = res$expanded$sigma$ranef[["V1"]], sigranef5 = res$expanded$sigma$ranef[["V4"]], muranef1 = res$expanded$mu$ranef[["V1"]], muranef5 = res$expanded$mu$ranef[["V4"]], LL=res$expanded$LL)
  return(ranefs)
}

wrapCalibrateBrainChartsIDQuantilePenalty <- function(simdata0, phenotype) {
  # this returns the estimate of the site effect

  res <- calibrateBrainChartsIDQuantilePenalty(simdata0, phenotype)
  ranefs <- c(sigranef1 = res$expanded$sigma$ranef[["V1"]], sigranef5 = res$expanded$sigma$ranef[["V4"]], muranef1 = res$expanded$mu$ranef[["V1"]], muranef5 = res$expanded$mu$ranef[["V4"]], LL=res$expanded$LL)
  return(ranefs)
}


# perform 100 resamples of the data
wrapperCIReal <- function(popdf, nshared, nsmall, nbig, wrapperFunc) {

    oneIter <- function(IDX) {
      # needed to avoid repeated sampling
      set.seed(IDX)
      popdf0 <- filter(popdf, obs == 0)
      popdf3 <- filter(popdf, obs == 3)

      P <- randperm(seq(1, nrow(popdf0)))
      #print(P)
      popdf03 <- bind_rows(popdf0[c(P[1:nshared], P[(nshared+1):(nshared+nsmall)]),], popdf3[c(P[1:nshared], P[(nshared+1):(nshared+nbig)]),])
      #popdf03 <- mutate(popdf03, study = paste0("V", INDEX.OB))

      real <- do.call(wrapperFunc, list(popdf03, phenotype = "CT"))
      #print(real)
      #browser()
      #real <- tibble(sigranef1 = real$expanded$sigma$ranef[["V1"]], muranef1 = real$expanded$mu$ranef[["V1"]], sigranef5 = real$expanded$sigma$ranef[["V4"]], muranef5 = real$expanded$mu$ranef[["V4"]])
      #real <- mutate(real, n=S)
    }
    #print("here")
    #boots <- future_map_dfr(bfits, ~wrapCalibrateBrainCharts(testdf, phenotype = "CT", phenofit = .x))
    allReals <- future_map(1:100, oneIter)
    allReals <- bind_rows(allReals)

    #print(allReals)
    boot.out.fake <- list(R=nrow(allReals))
    boot.out.fake$sim <- "ordinary"
    class(boot.out.fake) <- "boot"
    mu0V1.ci <- boot::boot.ci(boot.out = boot.out.fake, t=allReals$muranef1,
                              t0=mean(allReals$muranef1),
                              type=c("norm","basic", "perc"))
    sigma0V1.ci <- boot::boot.ci(boot.out = boot.out.fake, t=allReals$sigranef1,
                               t0=mean(allReals$sigranef1),
                               type=c("norm","basic", "perc"))
    mu0V4.ci <- boot::boot.ci(boot.out = boot.out.fake, t=allReals$muranef5,
                              t0=mean(allReals$muranef5),
                              type=c("norm","basic", "perc"))
    sigma0V4.ci <- boot::boot.ci(boot.out = boot.out.fake, t=allReals$sigranef5,
                                 t0=mean(allReals$sigranef5),
                                 type=c("norm","basic", "perc"))

    mtV1 <- tibble(conf = mu0V1.ci$normal[,"conf"],
                 mu=mu0V1.ci$t0,
                 mu.ci.low.normal = mu0V1.ci$normal[,2],
                 mu.ci.high.normal = mu0V1.ci$normal[,3],
                 mu.ci.low.basic = mu0V1.ci$basic[,4],
                 mu.ci.high.basic = mu0V1.ci$basic[,5],
                 mu.ci.low.percent = mu0V1.ci$percent[,4],
                 mu.ci.high.percent = mu0V1.ci$percent[,5]
    )
    mtV4 <- tibble(conf = mu0V4.ci$normal[,"conf"],
                   mu=mu0V4.ci$t0,
                   mu.ci.low.normal = mu0V4.ci$normal[,2],
                   mu.ci.high.normal = mu0V4.ci$normal[,3],
                   mu.ci.low.basic = mu0V4.ci$basic[,4],
                   mu.ci.high.basic = mu0V4.ci$basic[,5],
                   mu.ci.low.percent = mu0V4.ci$percent[,4],
                   mu.ci.high.percent = mu0V4.ci$percent[,5]
    )
    stV1 <- tibble(conf = sigma0V1.ci$normal[,"conf"],
                 sigma=sigma0V1.ci$t0,
                 sigma.ci.low.normal = sigma0V1.ci$normal[,2],
                 sigma.ci.high.normal = sigma0V1.ci$normal[,3],
                 sigma.ci.low.basic = sigma0V1.ci$basic[,4],
                 sigma.ci.high.basic = sigma0V1.ci$basic[,5],
                 sigma.ci.low.percent = sigma0V1.ci$percent[,4],
                 sigma.ci.high.percent = sigma0V1.ci$percent[,5]
    )
    stV4 <- tibble(conf = sigma0V4.ci$normal[,"conf"],
                   sigma=sigma0V4.ci$t0,
                   sigma.ci.low.normal = sigma0V4.ci$normal[,2],
                   sigma.ci.high.normal = sigma0V4.ci$normal[,3],
                   sigma.ci.low.basic = sigma0V4.ci$basic[,4],
                   sigma.ci.high.basic = sigma0V4.ci$basic[,5],
                   sigma.ci.low.percent = sigma0V4.ci$percent[,4],
                   sigma.ci.high.percent = sigma0V4.ci$percent[,5]

    )

  return(list(sigmaV1=stV1, muV1=mtV1, sigmaV4=stV4, muV4=mtV4))
  #res <- map_df(res, as.vector.data.frame)
  #res <- mutate(res, n=samps)
}

if(!file.exists("two_site_method_quantile_ground.rds")) {
  print("doing 04")
  ranef_estimations_with_CI_04 <- wrapCalibrateBrainChartsTwoTimepoint(simdata03, phenotype = "CT")
  saveRDS(ranef_estimations_with_CI_04, "two_site_method_quantile_ground.rds")
} else {
  ranef_estimations_with_CI_04 <- readRDS("two_site_method_quantile_ground.rds")
}

# use the two timepoint method without the IDs
TFileName <- paste0(opt$out, "T_nshared_", nshared, "_nsmall_", nsmall, "_nbig_", nbig, ".rds")
#if(!file.exists(TFileName)) {
  ranef_estimations_with_CI_04T <- wrapperCIReal(simdata, nshared, nsmall, nbig, wrapCalibrateBrainChartsTwoTimepoint)
  saveRDS(ranef_estimations_with_CI_04T, TFileName)
#}

IDFileName <- paste0(opt$out, "ID_nshared_", nshared, "_nsmall_", nsmall, "_nbig_", nbig, ".rds")
#if(!file.exists(IDFileName)) {
  ranef_estimations_with_CI_04ID <- wrapperCIReal(simdata, nshared, nsmall, nbig, wrapCalibrateBrainChartsIDQuantilePenalty)
  saveRDS(ranef_estimations_with_CI_04ID, IDFileName)
#}
