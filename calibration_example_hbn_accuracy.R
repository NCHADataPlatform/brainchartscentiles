# ---- Packages ----

library(tidyverse)
library(pracma)
source('calibration_example_hbn.R')
brainChartsDF <- data.frame(brainChartsDF)
brainChartsDF$mu.wre <- NA
brainChartsDF$sigma.wre <- NA
brainChartsDF$nu.wre <- NA
row.names(brainChartsDF) <- brainChartsDF$ID

brainChartsDF[row.names(fullSampleCBICCalibration$DATA.PRED2), 'mu.wre'] <- fullSampleCBICCalibration$DATA.PRED2[row.names(fullSampleCBICCalibration$DATA.PRED2), 'mu.wre']
brainChartsDF[row.names(fullSampleCBICCalibration$DATA.PRED2), 'sigma.wre'] <- fullSampleCBICCalibration$DATA.PRED2[row.names(fullSampleCBICCalibration$DATA.PRED2), "sigma.wre"]
brainChartsDF[row.names(fullSampleCBICCalibration$DATA.PRED2), 'nu.wre'] <- fullSampleCBICCalibration$DATA.PRED2[row.names(fullSampleCBICCalibration$DATA.PRED2), "nu.wre"]

brainChartsDF[row.names(fullSampleCUNYCalibration$DATA.PRED2), 'mu.wre'] <- fullSampleCUNYCalibration$DATA.PRED2[row.names(fullSampleCUNYCalibration$DATA.PRED2), "mu.wre"]
brainChartsDF[row.names(fullSampleCUNYCalibration$DATA.PRED2), 'sigma.wre'] <- fullSampleCUNYCalibration$DATA.PRED2[row.names(fullSampleCUNYCalibration$DATA.PRED2), "sigma.wre"]
brainChartsDF[row.names(fullSampleCUNYCalibration$DATA.PRED2), 'nu.wre'] <- fullSampleCUNYCalibration$DATA.PRED2[row.names(fullSampleCUNYCalibration$DATA.PRED2), "nu.wre"]

brainChartsDF[row.names(fullSampleRUCalibration$DATA.PRED2), 'mu.wre'] <- fullSampleRUCalibration$DATA.PRED2[row.names(fullSampleRUCalibration$DATA.PRED2), "mu.wre"]
brainChartsDF[row.names(fullSampleRUCalibration$DATA.PRED2), 'sigma.wre'] <- fullSampleRUCalibration$DATA.PRED2[row.names(fullSampleRUCalibration$DATA.PRED2), "sigma.wre"]
brainChartsDF[row.names(fullSampleRUCalibration$DATA.PRED2), 'nu.wre'] <- fullSampleRUCalibration$DATA.PRED2[row.names(fullSampleRUCalibration$DATA.PRED2), "nu.wre"]

brainChartsDF[row.names(fullSampleSICalibration$DATA.PRED2), 'mu.wre'] <- fullSampleSICalibration$DATA.PRED2[row.names(fullSampleSICalibration$DATA.PRED2), "mu.wre"]
brainChartsDF[row.names(fullSampleSICalibration$DATA.PRED2), 'sigma.wre'] <- fullSampleSICalibration$DATA.PRED2[row.names(fullSampleSICalibration$DATA.PRED2), "sigma.wre"]
brainChartsDF[row.names(fullSampleSICalibration$DATA.PRED2), 'nu.wre'] <- fullSampleSICalibration$DATA.PRED2[row.names(fullSampleSICalibration$DATA.PRED2), "nu.wre"]

# generate a new table with altered distribution parameters
brainChartsDFFake <- brainChartsDF

# add timepoint effects
muNudge <- 0.1
sigmaNudge <- 0.1
brainChartsDFFake$mu.wre <- brainChartsDF$mu.wre + muNudge
brainChartsDFFake$sigma.wre <- brainChartsDF$sigma.wre + sigmaNudge

# add a year to the ages and some randomness
brainChartsDFFake$age_days <- brainChartsDFFake$age_days + 365.25 + abs(rnorm(nrow(brainChartsDFFake), mean = 0, sd = 10))
# generate variables from altered distributions
brainChartsDFFake$CT <- rGGalt(nrow(brainChartsDFFake), exp(brainChartsDFFake$mu.wre), exp(brainChartsDFFake$sigma.wre), brainChartsDFFake$nu.wre) * 10000

brainChartsDFFake$study <- paste0(brainChartsDFFake$study, "2")

# calibrate the full fake samples
# in a normal example, these would be your small samples
fakeSampleCBICCalibration <- calibrateBrainCharts(filter(brainChartsDFFake, study == "HBNCBIC2"), phenotype = phenotype)
fakeSampleCUNYCalibration <- calibrateBrainCharts(filter(brainChartsDFFake, study == "HBNCUNY2"), phenotype = phenotype)
fakeSampleSICalibration <- calibrateBrainCharts(filter(brainChartsDFFake, study == "HBNSI2"), phenotype = phenotype)
fakeSampleRUCalibration <- calibrateBrainCharts(filter(brainChartsDFFake, study == "HBNRU2"), phenotype = phenotype)


n <- 50
subFakeCBIC <- filter(brainChartsDFFake, study == "HBNCBIC2")
P <- randperm(1:nrow(subFakeCBIC))
subFakeCBIC <- subFakeCBIC[P[1:n],]
subFakeCUNY <- filter(brainChartsDFFake, study == "HBNCUNY2")
P <- randperm(1:nrow(subFakeCUNY))
subFakeCUNY <- subFakeCUNY[P[1:n],]
subFakeRU <- filter(brainChartsDFFake, study == "HBNRU2")
P <- randperm(1:nrow(subFakeRU))
subFakeRU <- subFakeRU[P[1:n],]
subFakeSI <- filter(brainChartsDFFake, study == "HBNSI2")
P <- randperm(1:nrow(subFakeSI))
subFakeSI <- subFakeSI[P[1:n],]

# accuracy evaluation

# run the centile-informed calibration

subFakeCBICQuantileCalibration <- calibrateBrainChartsIDQuantilePenalty(subFakeCBIC, phenotype = "CT", largeSiteOutput = fullSampleCBICCalibration)
subFakeCUNYQuantileCalibration <- calibrateBrainChartsIDQuantilePenalty(subFakeCUNY, phenotype = "CT", largeSiteOutput = fullSampleCUNYCalibration)
subFakeSIQuantileCalibration <- calibrateBrainChartsIDQuantilePenalty(subFakeSI, phenotype = "CT", largeSiteOutput = fullSampleSICalibration)
subFakeRUQuantileCalibration <- calibrateBrainChartsIDQuantilePenalty(subFakeRU, phenotype = "CT", largeSiteOutput = fullSampleRUCalibration)

# run the normal Braincharts calibration on the small sites 
subFakeCBICCalibration <- calibrateBrainCharts(subFakeCBIC, phenotype = "CT")
subFakeCUNYCalibration <- calibrateBrainCharts(subFakeCUNY, phenotype = "CT")
subFakeRUCalibration <- calibrateBrainCharts(subFakeRU, phenotype = "CT")
subFakeSICalibration <- calibrateBrainCharts(subFakeSI, phenotype = "CT")

# make plots 

fullSampleCBICCalibration$DATA.PRED2$fitting <- 'large sample (cross-sectional method)'
subFakeCBICQuantileCalibration$DATA.PRED2$fitting <- 'small site (centile method)'
subFakeCBICCalibration$DATA.PRED2$fitting <- 'small site (cross-sectional method)'
fakeSampleCBICCalibration$DATA.PRED2$fitting <- 'small site (ground truth)'

fullSampleCBICCalibration$DATA.PRED2$sample <- "large"
fullSampleCBICCalibration$DATA.PRED2
subFakeCBICQuantileCalibration$DATA.PRED2$sample <- "small"
subFakeCBICCalibration$DATA.PRED2$sample <- "small"
fakeSampleCBICCalibration$DATA.PRED2$sample <- "small"

T <- bind_rows(fullSampleCBICCalibration$DATA.PRED2, subFakeCBICQuantileCalibration$DATA.PRED2, subFakeCBICCalibration$DATA.PRED2, fakeSampleCBICCalibration$DATA.PRED2)

#mapping <- c(V1 = "large", V5 = "small")
#T$sample <- mapping[T$study]

# dont need to plot the ground truth CT measures for the small site
T$meanCT2Transformed[T$fitting == "small site (ground truth)"] <- NA

# meanCT2Transformed is the transformed cortical thickness
# PRED.l250.wre, PRED.m500.wre, PRED.u750.wre are the estimated 25th, 50th, and 97.5th Quantiles after site-effect corrections

ggplot(T, aes(x = age_days)) +
  geom_point(aes(y = meanCT2Transformed * 10000, shape = sample), alpha = 0.5) +
  geom_line(aes(y = PRED.l250.wre * 10000, colour = fitting)) + 
  geom_line(aes(y = PRED.m500.wre * 10000, colour = fitting), linewidth = 1) +
  geom_line(aes(y = PRED.u750.wre * 10000, colour = fitting)) +
  labs(x = "Age (years)", y = "CT", colour = "Estimation Method")

ggsave('freesurfer_example_hbn_CBIC.png')