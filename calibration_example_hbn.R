# ---- Packages ----

library(tidyverse)

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

# ---- StudyOverview ----

brainChartsDF <- mutate(brainChartsDF, study = factor(study, levels=c("HBNCBIC", "HBNCUNY", "HBNRU", "HBNSI")))
count(brainChartsDF, study)
count(brainChartsDF, study, sex, sort = FALSE)

brainChartsDF <- filter(brainChartsDF, study == "HBNCBIC")
brainChartsDF <- data.frame(brainChartsDF)
row.names(brainChartsDF) <- brainChartsDF$ID

phenotype <- "CT"

# do the full sample calibrations per site
fullSampleCBICCalibration <- calibrateBrainCharts(filter(brainChartsDF, study == "HBNCBIC"), phenotype = phenotype)
# fullSampleCUNYCalibration <- calibrateBrainCharts(filter(brainChartsDF, study == "HBNCUNY"), phenotype = phenotype)
# fullSampleSICalibration <- calibrateBrainCharts(filter(brainChartsDF, study == "HBNSI"), phenotype = phenotype)
# fullSampleRUCalibration <- calibrateBrainCharts(filter(brainChartsDF, study == "HBNRU"), phenotype = phenotype)

row.names(fullSampleCBICCalibration$DATA.PRED2) <- fullSampleCBICCalibration$DATA.PRED2$ID
# row.names(fullSampleCUNYCalibration$DATA.PRED2) <- fullSampleCUNYCalibration$DATA.PRED2$ID
# row.names(fullSampleSICalibration$DATA.PRED2) <- fullSampleSICalibration$DATA.PRED2$ID
# row.names(fullSampleRUCalibration$DATA.PRED2) <- fullSampleRUCalibration$DATA.PRED2$ID

brainChartsDF$mu.wre <- NA
brainChartsDF$sigma.wre <- NA
brainChartsDF$nu.wre <- NA
row.names(brainChartsDF) <- brainChartsDF$ID

brainChartsDF[row.names(fullSampleCBICCalibration$DATA.PRED2), 'mu.wre'] <- fullSampleCBICCalibration$DATA.PRED2[row.names(fullSampleCBICCalibration$DATA.PRED2), 'mu.wre']
brainChartsDF[row.names(fullSampleCBICCalibration$DATA.PRED2), 'sigma.wre'] <- fullSampleCBICCalibration$DATA.PRED2[row.names(fullSampleCBICCalibration$DATA.PRED2), "sigma.wre"]
brainChartsDF[row.names(fullSampleCBICCalibration$DATA.PRED2), 'nu.wre'] <- fullSampleCBICCalibration$DATA.PRED2[row.names(fullSampleCBICCalibration$DATA.PRED2), "nu.wre"]

# brainChartsDF[row.names(fullSampleCUNYCalibration$DATA.PRED2), 'mu.wre'] <- fullSampleCUNYCalibration$DATA.PRED2[row.names(fullSampleCUNYCalibration$DATA.PRED2), "mu.wre"]
# brainChartsDF[row.names(fullSampleCUNYCalibration$DATA.PRED2), 'sigma.wre'] <- fullSampleCUNYCalibration$DATA.PRED2[row.names(fullSampleCUNYCalibration$DATA.PRED2), "sigma.wre"]
# brainChartsDF[row.names(fullSampleCUNYCalibration$DATA.PRED2), 'nu.wre'] <- fullSampleCUNYCalibration$DATA.PRED2[row.names(fullSampleCUNYCalibration$DATA.PRED2), "nu.wre"]

# brainChartsDF[row.names(fullSampleRUCalibration$DATA.PRED2), 'mu.wre'] <- fullSampleRUCalibration$DATA.PRED2[row.names(fullSampleRUCalibration$DATA.PRED2), "mu.wre"]
# brainChartsDF[row.names(fullSampleRUCalibration$DATA.PRED2), 'sigma.wre'] <- fullSampleRUCalibration$DATA.PRED2[row.names(fullSampleRUCalibration$DATA.PRED2), "sigma.wre"]
# brainChartsDF[row.names(fullSampleRUCalibration$DATA.PRED2), 'nu.wre'] <- fullSampleRUCalibration$DATA.PRED2[row.names(fullSampleRUCalibration$DATA.PRED2), "nu.wre"]

# brainChartsDF[row.names(fullSampleSICalibration$DATA.PRED2), 'mu.wre'] <- fullSampleSICalibration$DATA.PRED2[row.names(fullSampleSICalibration$DATA.PRED2), "mu.wre"]
# brainChartsDF[row.names(fullSampleSICalibration$DATA.PRED2), 'sigma.wre'] <- fullSampleSICalibration$DATA.PRED2[row.names(fullSampleSICalibration$DATA.PRED2), "sigma.wre"]
# brainChartsDF[row.names(fullSampleSICalibration$DATA.PRED2), 'nu.wre'] <- fullSampleSICalibration$DATA.PRED2[row.names(fullSampleSICalibration$DATA.PRED2), "nu.wre"]

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

n <- 50
subFakeCBIC <- filter(brainChartsDFFake, study == "HBNCBIC2")
P <- randperm(1:nrow(subFakeCBIC))
subFakeCBIC <- subFakeCBIC[P[1:n],]
# subFakeCUNY <- filter(brainChartsDFFake, study == "HBNCUNY2")
# P <- randperm(1:nrow(subFakeCUNY))
# subFakeCUNY <- subFakeCUNY[P[1:n],]
# subFakeRU <- filter(brainChartsDFFake, study == "HBNRU2")
# P <- randperm(1:nrow(subFakeRU))
# subFakeRU <- subFakeRU[P[1:n],]
# subFakeSI <- filter(brainChartsDFFake, study == "HBNSI2")
# P <- randperm(1:nrow(subFakeSI))
# subFakeSI <- subFakeSI[P[1:n],]

subFakeCBICQuantileCalibration <- calibrateBrainChartsIDQuantilePenalty(subFakeCBIC, phenotype = "CT", largeSiteOutput = fullSampleCBICCalibration)
# subFakeCUNYQuantileCalibration <- calibrateBrainChartsIDQuantilePenalty(subFakeCUNY, phenotype = "CT", largeSiteOutput = fullSampleCUNYCalibration)
# subFakeSIQuantileCalibration <- calibrateBrainChartsIDQuantilePenalty(subFakeSI, phenotype = "CT", largeSiteOutput = fullSampleSICalibration)
# subFakeRUQuantileCalibration <- calibrateBrainChartsIDQuantilePenalty(subFakeRU, phenotype = "CT", largeSiteOutput = fullSampleRUCalibration)
