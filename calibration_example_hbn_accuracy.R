# ---- Packages ----

library(tidyverse)
library(pracma)
source('calibration_example_hbn.R')

# accuracy evaluation

# calibrate the full fake samples
fakeSampleCBICCalibration <- calibrateBrainCharts(filter(brainChartsDFFake, study == "HBNCBIC2"), phenotype = phenotype)
# fakeSampleCUNYCalibration <- calibrateBrainCharts(filter(brainChartsDFFake, study == "HBNCUNY2"), phenotype = phenotype)
# fakeSampleSICalibration <- calibrateBrainCharts(filter(brainChartsDFFake, study == "HBNSI2"), phenotype = phenotype)
# fakeSampleRUCalibration <- calibrateBrainCharts(filter(brainChartsDFFake, study == "HBNRU2"), phenotype = phenotype)

# run the centile-informed calibration

# run the normal Braincharts calibration on the small sites 
subFakeCBICCalibration <- calibrateBrainCharts(subFakeCBIC, phenotype = "CT")
# subFakeCUNYCalibration <- calibrateBrainCharts(subFakeCUNY, phenotype = "CT")
# subFakeRUCalibration <- calibrateBrainCharts(subFakeRU, phenotype = "CT")
# subFakeSICalibration <- calibrateBrainCharts(subFakeSI, phenotype = "CT")

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