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

# nshared <- opt$nshared
# nbig <- opt$nbig
# nsmall <- opt$nsmall

# select two "sites", timepoints 0 and 4
#simdata04 <- filter(simdata, obs == 0 | obs == 4)
simdata0 <- filter(simdata, obs == 0 & dx == "CN")
simdata4 <- filter(simdata, obs == 4 & dx == "CN")

simdata0 <- mutate(simdata0, study = paste0("V", INDEX.OB))
simdata4 <- mutate(simdata4, study = paste0("V", INDEX.OB))

# fit for the main site
base0 <- calibrateBrainCharts(simdata0, phenotype = "CT")

# mu and sigma will be in
# $base0$expanded$mu$ranef["V1"], base0$expanded$sigma$ranef["V1"]

#baseparams <- list(muranef = base0$expanded$mu$ranef["V1"], sigmaranef = base0$expanded$sigma$ranef["V1"])

# get the shared IDs
allIDs <- c(simdata0$ID, simdata4$ID)
sharedIDs <- allIDs[duplicated(allIDs)]

# pick random selection of subjects from time point 4
nsmall <- 10
P <- randperm(1:length(sharedIDs), nsmall)

sample4 <- simdata4[simdata4$ID %in% sharedIDs[P],]

base4 <- calibrateBrainChartsIDQuantilePenalty(sample4, phenotype = "CT", largeSiteOutput = base0)

# mu and sigma will be in
# base4$expanded$mu$ranef["V1"], base4$expanded$sigma$ranef["V1"]

# for cortical thickness, the estimated quantiles for the altered distributions after site-effect corrections are
# base4$DATA.PRED2$meanCT2Transformed.q.wre
