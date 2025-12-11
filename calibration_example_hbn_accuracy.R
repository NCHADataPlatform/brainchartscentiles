# ---- Packages ----

library(tidyverse)
library(pracma)
source('calibration_example_hbn.R')

# accuracy evaluation


# run the centile-informed calibration

# run the normal Braincharts calibration on the small sites 
# ---- Calib1 ----
subFakeCalibration <- calibrateBrainCharts(subFakeDF, phenotype = "CT")

# ---- Calib3 ----
fakeSampleCalibration <- calibrateBrainCharts(brainChartsDFFake, phenotype = phenotype)


# ---- MakePlots ----
fullSampleCalibration$DATA.PRED2$fitting <- 'large sample (cross-sectional method)'
subFakeQuantileCalibration$DATA.PRED2$fitting <- 'small site (centile method)'
subFakeCalibration$DATA.PRED2$fitting <- 'small site (cross-sectional method)'
fakeSampleCalibration$DATA.PRED2$fitting <- 'small site (ground truth)'

fullSampleCalibration$DATA.PRED2$sample <- "large"
fullSampleCalibration$DATA.PRED2
subFakeQuantileCalibration$DATA.PRED2$sample <- "small"
subFakeCalibration$DATA.PRED2$sample <- "small"
fakeSampleCalibration$DATA.PRED2$sample <- "small"

T <- bind_rows(fullSampleCalibration$DATA.PRED2, subFakeQuantileCalibration$DATA.PRED2, subFakeCalibration$DATA.PRED2, fakeSampleCalibration$DATA.PRED2)

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