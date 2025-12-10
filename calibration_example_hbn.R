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
  colnames(df)[1] <- "ID"
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

# ---- MakeDXColumn ----
brainChartsDF <- mutate(brainChartsDF, 
                        dx = "CN"
)

# ---- LoadDemo ----

demoCBIC <- read_csv(file.path('HBN', 'participants-CBIC.csv'), 
                     show_col_types = FALSE)
demoCUNY <- read_csv(file.path('HBN', 'participants-CUNY.csv'), 
                     show_col_types = FALSE)
demoSI <- read_csv(file.path('HBN', 'participants-SI.csv'), 
                   show_col_types = FALSE)
demoRU <- read_csv(file.path('HBN', 'participants-RU.csv'), 
                   show_col_types = FALSE)

demoCBIC <- mutate(demoCBIC, study = 'HBNCBIC')
demoCUNY <- mutate(demoCUNY,study = 'HBNCUNY')
demoSI <- mutate(demoSI, study = 'HBNSI')
demoRU <- mutate(demoRU , study = 'HBNRU')

demoDF <- bind_rows(demoCBIC, demoCUNY, demoSI, demoRU)
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
brainChartsDF <- mutate(brainChartsDF, study = factor(study, levels=c("HBNCBIC", "HBNCUNY", "HBNRU", "HBNSI")))

# ---- StudyOverview ----

count(brainChartsDF, study)
count(brainChartsDF, study, sex, sort = FALSE)

# ---- CalibrationPhase1 ----
brainChartsDF <- filter(brainChartsDF, study == "HBNCBIC")
# convert to dataframe for compatability with non-tidyverse code
brainChartsDF <- data.frame(brainChartsDF)
row.names(brainChartsDF) <- brainChartsDF$ID

phenotype <- "CT"

# do the full sample calibrations per site
fullSampleCBICCalibration <- calibrateBrainCharts(
  filter(brainChartsDF, study == "HBNCBIC"), phenotype = phenotype)

# ---- SimulateSmallSiteSetup ----


simulateSite <- function(fullSampleCalibration, BC) {
  row.names(fullSampleCalibration$DATA.PRED2) <- fullSampleCBICCalibration$DATA.PRED2$ID
  
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
  n <- 50
  subFakeCBIC <- filter(BCFake, study == "HBNCBIC2")
  P <- randperm(1:nrow(subFakeCBIC))
  subFakeCBIC <- subFakeCBIC[P[1:n],]
  return(subFakeCBIC)
}
# ---- SimulateCompleteSite ----
brainChartsDFFake <- simulateSite(fullSampleCBICCalibration, brainChartsDF)
# ---- SimulateSmallSite ----

subFakeCBIC <- simulateSmallSite(fullSampleCBICCalibration, brainChartsDF)
# ---- CalibrationPhase2 ----

subFakeCBICQuantileCalibration <- 
  calibrateBrainChartsIDQuantilePenalty(subFakeCBIC, phenotype = "CT", 
                                        largeSiteOutput = fullSampleCBICCalibration)
