# ---- Packages ----

library(tidyverse)
library(pracma)

# ---- Setup ----

BrainChartFolder <- dirname(normalizePath("calibrate_braincharts_centiles.R"))
wd <- setwd(BrainChartFolder)
source("calibrate_braincharts_centiles.R")
setwd(wd)

# ---- LoadingStructuralData ----
readFS <- function(fname) {
  df <- read_csv(fname, name_repair = "universal", show_col_types = FALSE)
  df <- rename(df, ID = 1)
  return(df)
}

asegDF <- readFS(file.path('HBN', 'stats_aseg.csv'))
lhCortThicknessDF <- readFS(file.path('HBN', 'stats_aparc_lh_thickness.csv'))
lhCortSurfaceAreaDF <- readFS(file.path('HBN', 'stats_aparc_lh_area.csv'))
rhCortThicknessDF <- readFS(file.path('HBN', 'stats_aparc_rh_thickness.csv'))
rhCortSurfaceAreaDF <- readFS(file.path('HBN', 'stats_aparc_rh_area.csv'))

# ---- MakeIDPColumns ----
# make the IDP columns
brainChartsDF <- select(asegDF, ID, 
                        CerebralWhiteMatterVol, 
                        SubCortGrayVol, CortexVol, 
                        EstimatedTotalIntraCranialVol, 
                        Left.Lateral.Ventricle, Right.Lateral.Ventricle, 
                        SupraTentorialVolNotVent, BrainSegVol.to.eTIV)
brainChartsDF <- mutate(brainChartsDF, 
                        Ventricles = Left.Lateral.Ventricle + 
                          Right.Lateral.Ventricle, BrainSegVol.to.eTIV)

brainChartsDF <- rename(brainChartsDF, 
                        WMV = CerebralWhiteMatterVol,
                        sGMV = SubCortGrayVol,
                        GMV = CortexVol,
                        etiv = EstimatedTotalIntraCranialVol,
                        TCV = SupraTentorialVolNotVent,
                        )
# join surface area measures
brainChartsDF <- left_join(brainChartsDF, 
                           select(lhCortSurfaceAreaDF, ID, 
                                  lh_WhiteSurfArea_area), by="ID")
brainChartsDF <- left_join(brainChartsDF, 
                           select(rhCortSurfaceAreaDF, ID, 
                                  rh_WhiteSurfArea_area), by="ID")
# join cortical thickness measures
brainChartsDF <- left_join(brainChartsDF, 
                           select(lhCortThicknessDF, ID, 
                                  lh_MeanThickness_thickness), by="ID")
brainChartsDF <- left_join(brainChartsDF, 
                           select(rhCortThicknessDF, ID, 
                                  rh_MeanThickness_thickness), by="ID")

brainChartsDF <- mutate(brainChartsDF, 
                        SA = lh_WhiteSurfArea_area + rh_WhiteSurfArea_area)
brainChartsDF <- mutate(brainChartsDF, 
                        CT = ((lh_MeanThickness_thickness * 
                                 lh_WhiteSurfArea_area) + 
                                (rh_MeanThickness_thickness * 
                                   rh_WhiteSurfArea_area))/SA,
                        fs_version = 'Custom_T1T2',
                        run = 1,
                        participant = ID,
                        country = 'Australia',
                        dx = 'CN'
                        )

brainChartsDF <- select(brainChartsDF, -Right.Lateral.Ventricle, 
                        -Left.Lateral.Ventricle, 
                        -ends_with("area"), -ends_with("thickness"))
# ---- MakeOtherColumns ----


brainChartsDF <- mutate(brainChartsDF, 
                        fs_version = 'Custom_T1T2',
                        run = 1,
                        participant = ID,
                        country = 'Australia',
)

brainChartsDF <- filter(brainChartsDF, CT > 2)

# ---- MakeDXColumn ----
brainChartsDF <- mutate(brainChartsDF, 
                        dx = "CN"
)

# make the last 5 subjects non-controls

brainChartsDF$dx[(nrow(brainChartsDF)-4):nrow(brainChartsDF)] <- "notCN"

# ---- LoadDemo ----

demoDF <- read_csv(file.path('HBN', 'participants-CBIC.csv'), 
                     show_col_types = FALSE)

demoDF <- mutate(demoDF, study = 'HBNCBIC')

demoDF <- mutate(demoDF, age_days = Age * 365.25)
demoDF <- rename(demoDF, sex = Sex)

demoDF <- mutate(demoDF, 
                 sex = factor(sex, labels = c("Male", "Female"), 
                                      levels = c(0, 1)))
# ---- JoinDemographics ----
brainChartsDF <- left_join(brainChartsDF, 
                           select(demoDF, participant_id, sex, age_days, study), 
                           by=join_by(ID==participant_id))

# remove mising values
brainChartsDF <- na.omit(brainChartsDF)
brainChartsDF <- mutate(brainChartsDF, study = factor(study, levels=c("HBNCBIC")))

# ---- StudyOverview ----

count(brainChartsDF, study)
count(brainChartsDF, study, sex, sort = FALSE)

# ---- CalibrationPhase1 ----
# convert to dataframe for compatability with non-tidyverse code
brainChartsDF <- data.frame(brainChartsDF)
row.names(brainChartsDF) <- brainChartsDF$ID

phenotype <- "CT"

# do the full sample calibrations per site
fullSampleCalibration <- calibrateBrainCharts(
  brainChartsDF, phenotype = phenotype)

# ---- SimulateSmallSiteSetup ----


simulateSite <- function(fullSampleCalibration, BC) {
  row.names(fullSampleCalibration$DATA.PRED2) <- fullSampleCalibration$DATA.PRED2$ID
  
  BC$mu.wre <- NA
  BC$sigma.wre <- NA
  BC$nu.wre <- NA
  row.names(BC) <- BC$ID
  
  BC[row.names(fullSampleCalibration$DATA.PRED2), 'mu.wre'] <- fullSampleCalibration$DATA.PRED2[row.names(fullSampleCalibration$DATA.PRED2), 'mu.wre']
  BC[row.names(fullSampleCalibration$DATA.PRED2), 'sigma.wre'] <- fullSampleCalibration$DATA.PRED2[row.names(fullSampleCalibration$DATA.PRED2), "sigma.wre"]
  BC[row.names(fullSampleCalibration$DATA.PRED2), 'nu.wre'] <- fullSampleCalibration$DATA.PRED2[row.names(fullSampleCalibration$DATA.PRED2), "nu.wre"]
  
  # generate a new table with altered distribution parameters
  BCFake <- BC
  # add timepoint effects
  muNudge <- 0.1
  sigmaNudge <- 0.1
  BCFake$mu.wre <- BC$mu.wre + muNudge
  BCFake$sigma.wre <- BC$sigma.wre + sigmaNudge
  
  # add a year to the ages and some randomness
  BCFake$age_days <- BCFake$age_days + 365.25 + abs(rnorm(nrow(BCFake), mean = 0, sd = 10))
  # generate variables from altered distributions
  BCFake$CT <- rGGalt(nrow(BCFake), exp(BCFake$mu.wre), exp(BCFake$sigma.wre), BCFake$nu.wre) * 10000
  
  BCFake$study <- paste0(BCFake$study, "2")
  # add some patients too.
  return(BCFake)
 }

simulateSmallSite <- function(fullSampleCalibration, BC) {
  BCFake <- simulateSite(fullSampleCalibration, BC)
  BCFakeOrig <- BCFake
  n <- 50
  I <- which(BC$dx == "CN")
  P <- randperm(1:length(I))
  BCFake <- BCFake[I[P[1:n]],]
  BCFake <- bind_rows(BCFake, BCFakeOrig[BCFakeOrig$dx == "notCN",])
  return(BCFake)
}
# ---- SimulateCompleteSite ----
brainChartsDFFake <- simulateSite(fullSampleCalibration, brainChartsDF)
# ---- SimulateSmallSite ----

subFakeDF <- simulateSmallSite(fullSampleCalibration, brainChartsDF)
# ---- CalibrationPhase2 ----

subFakeQuantileCalibration <- 
  calibrateBrainChartsIDQuantilePenalty(subFakeDF, phenotype = "CT", 
                                        largeSiteOutput = fullSampleCalibration)


# ---- ResultsPrint -----
fullSampleCalibration$DATA.PRED2$sample <- "large"
subFakeQuantileCalibration$DATA.PRED2$sample <- "small"

I <- which(fullSampleCalibration$DATA.PRED2$dx == "CN")
J <- which(fullSampleCalibration$DATA.PRED2$dx == "notCN")
S <- bind_rows(fullSampleCalibration$DATA.PRED2[I[1],], fullSampleCalibration$DATA.PRED2[J[1],])

S[, c('dx', 'age_days', 
      'PRED.l025.pop', 'PRED.l250.pop', 'PRED.m500.pop', 'PRED.u750.pop', 'PRED.u975.pop',
      'PRED.l025.wre', 'PRED.l250.wre', 'PRED.m500.wre', 'PRED.u750.wre', 'PRED.u975.wre',
      'meanCT2Transformed.q.wre', 'meanCT2Transformed.normalised')]

# ---- MakePlots -----
T <- bind_rows(fullSampleCalibration$DATA.PRED2, subFakeQuantileCalibration$DATA.PRED2)
M <- c(CN = 1, notCN = 5)

T$DXSize <- factor(M[T$dx], labels=c("CN", "notCN"))

ggplot(T, aes(x = age_days / 365.25)) +
  geom_point(aes(y = meanCT2Transformed * 10000, colour=sample, size = DXSize), alpha = 0.5) +
  geom_line(aes(y = PRED.l250.wre * 10000, colour=sample)) + 
  geom_line(aes(y = PRED.m500.wre * 10000, colour=sample), linewidth = 1) +
  geom_line(aes(y = PRED.u750.wre * 10000, colour=sample)) +
  labs(x = "Age (years)", y = "CT", colour = "Sample")