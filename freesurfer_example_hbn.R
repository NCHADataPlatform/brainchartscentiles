source("calibrate_braincharts_centiles.R")

library('pracma')

# read the Freesurfer stats tables, change the paths as required
asegDF <- read.csv(file.path('HBN', 'stats_aseg.csv'))
lhCortThicknessDF <- read.csv(file.path('HBN', 'stats_aparc_lh_thickness.csv'))
lhCortSurfaceAreaDF <- read.csv(file.path('HBN', 'stats_aparc_lh_area.csv'))
rhCortThicknessDF <- read.csv(file.path('HBN', 'stats_aparc_lh_thickness.csv')) ## wrong
rhCortSurfaceAreaDF <- read.csv(file.path('HBN', 'stats_aparc_rh_area.csv'))

# set the row names to participant IDs
row.names(asegDF) <- asegDF$Measure.volume
row.names(lhCortSurfaceAreaDF) <- lhCortSurfaceAreaDF$lh.aparc.area
row.names(rhCortSurfaceAreaDF) <- lhCortSurfaceAreaDF$rh.aparc.area
row.names(lhCortThicknessDF) <- lhCortThicknessDF$lh.aparc.area
row.names(rhCortThicknessDF) <- lhCortThicknessDF$rh.aparc.area

brainChartsDF <- data.frame(row.names=asegDF$Measure.volume)

# make the IDP columns
brainChartsDF$WMV <- asegDF$CerebralWhiteMatterVol
brainChartsDF$sGMV <- asegDF$SubCortGrayVol
brainChartsDF$GWMV <- asegDF$CortexVol
brainChartsDF$etiv <- asegDF$eTIV # wrong
brainChartsDF$Ventricles <- asegDF$VentricleChoroidVol
brainChartsDF$TCV <- asegDF$SupraTentorialVolNotVent
brainChartsDF$SA <- lhCortSurfaceAreaDF$lh_WhiteSurfArea_area + rhCortSurfaceAreaDF$rh_WhiteSurfArea_area
brainChartsDF$CT <- (lhCortThicknessDF$lh_MeanThickness_thickness * lhCortSurfaceAreaDF$lh_WhiteSurfArea_area + lhCortThicknessDF$lh_MeanThickness_thickness * rhCortSurfaceAreaDF$rh_WhiteSurfArea_area) / brainChartsDF$SA
brainChartsDF$BrainSegVol.to.eTIV <- asegDF$BrainSegVol.to.eTIV

# put in other columns, leave as is
brainChartsDF$fs_version <- 'Custom_T1T2'
brainChartsDF$run <- 1
brainChartsDF$participant <- row.names(brainChartsDF)
brainChartsDF$country <- 'Australia'

# diagnosis, CN for control (for calibration), notCN for patient (not for calibration)
brainChartsDF$dx <- 'CN'

# study, must be different than the studies in the training set, which are
# 3R-BRAIN, ABCD, abide1, abide2, ADHD200, ADNI, AIBL, AOBA, AOMIC_ID1000,
# AOMIC_PIOP1, AOMIC_PIOP2, ARWIBO, BabyConnectome, BGSP, BHRCS, BSNIP,
# Calgary, CALM, CamCAN, CAMFT, Cornell_C1, Cornell_C2, cVEDA, devCCNP, dHCP,
# DLBS, EDSD, EMBARC, FemaleASD, FinnBrain, GUSTO, HABS, Harvard_Fetal1, HBN, HCP,
# HCP_lifespanA, HCP_lifespanD, IBIS, ICBM, IMAGEN, IMAP, IXI, LA5c, MCIC, Narratives,
# NHGRI, NIH, NIHPD, NSPN, NTB_Yale, OASIS3, Oulu, PCDC, PING, PNC, POND, PREVENTAD,
# RDB, SALD, SLIM, UCSD, UKB, VETSA, WAYNE

brainChartsDF$study <- 'MyStudy'

# demographics needed
# brainChartsDF$age_days  age in days
# brainChartsDF$sex  age in days
# brainChartsDF$sex <- factor(brainChartsDF$sex, labels = c("Male", "Female"), levels = c("M", "F"))

# the remaining code gets the demographic information
# this is specific to the HBN dataset, and will likely not apply to other datasets

# HBN example
demoCBIC <- read.csv(file.path('HBN', 'participants-CBIC.csv'))
demoCUNY <- read.csv(file.path('HBN', 'participants-CUNY.csv'))
demoSI <- read.csv(file.path('HBN', 'participants-SI.csv'))
demoRU <- read.csv(file.path('HBN', 'participants-RU.csv'))

row.names(demoCBIC) <- demoCBIC$participant_id
row.names(demoCUNY) <- demoCUNY$participant_id
row.names(demoSI) <- demoSI$participant_id
row.names(demoRU) <- demoRU$participant_id

# setting studies for individual sites
demoCBIC$study <- 'HBNCBIC'
demoCUNY$study <- 'HBNCUNY'
demoSI$study <- 'HBNSI'
demoRU$study <- 'HBNRU'

# merge all
demoDF <- bind_rows(demoCBIC, demoCUNY, demoSI, demoRU)

# insert age and sex
validRowNames <- row.names(brainChartsDF) %in% row.names(demoDF)
validRowNames <- row.names(brainChartsDF)[validRowNames]
brainChartsDF[validRowNames, 'age_days'] <- demoDF[validRowNames, 'Age'] * 365.25
brainChartsDF[validRowNames, 'sex'] <- demoDF[validRowNames, 'Sex']
brainChartsDF[validRowNames, 'study'] <- demoDF[validRowNames, 'study']

# make sure the sex field is correct
brainChartsDF$sex[brainChartsDF$sex == 0] <- 'M'
brainChartsDF$sex[brainChartsDF$sex == 1] <- 'F'

# this needs to be done for your dataset
brainChartsDF$sex <- factor(brainChartsDF$sex, labels = c("Male", "Female"), levels = c("M", "F"))

# removing NA rows
brainChartsDF <- brainChartsDF[complete.cases(brainChartsDF),]

# preparation is now complete

# choose cortical thickness as phenotype
phenotype <- "CT"

# print site summary
#print(summary(brainChartsDF[,c('study', 'age_days')]))

brainChartsDF$study <- factor(brainChartsDF$study, levels=c("HBNCBIC", "HBNCUNY", "HBNRU", "HBNSI"))
print("Site numbers")
print(table(brainChartsDF[,c('study')]))
print(brainChartsDF %>% count(study, sex, sort = FALSE))

# do the full sample calibrations per site
fullSampleCBIC <- calibrateBrainCharts(filter(brainChartsDF, study == "HBNCBIC"), phenotype = phenotype)
fullSampleCUNY <- calibrateBrainCharts(filter(brainChartsDF, study == "HBNCUNY"), phenotype = phenotype)
fullSampleSI <- calibrateBrainCharts(filter(brainChartsDF, study == "HBNSI"), phenotype = phenotype)
fullSampleRU <- calibrateBrainCharts(filter(brainChartsDF, study == "HBNRU"), phenotype = phenotype)

# for each site, make a longitudinal second session with altered site parameters, generate new cortical thickness measures

# these wre variables are the sample-corrected distribution parameters for each measurement

brainChartsDF$mu.wre <- NA
brainChartsDF$sigma.wre <- NA
brainChartsDF$nu.wre <- NA

brainChartsDF[row.names(fullSampleCBIC$DATA.PRED2), 'mu.wre'] <- fullSampleCBIC$DATA.PRED2$mu.wre
brainChartsDF[row.names(fullSampleCBIC$DATA.PRED2), 'sigma.wre'] <- fullSampleCBIC$DATA.PRED2$sigma.wre
brainChartsDF[row.names(fullSampleCBIC$DATA.PRED2), 'nu.wre'] <- fullSampleCBIC$DATA.PRED2$nu.wre

brainChartsDF[row.names(fullSampleCUNY$DATA.PRED2), 'mu.wre'] <- fullSampleCUNY$DATA.PRED2$mu.wre
brainChartsDF[row.names(fullSampleCUNY$DATA.PRED2), 'sigma.wre'] <- fullSampleCUNY$DATA.PRED2$sigma.wre
brainChartsDF[row.names(fullSampleCUNY$DATA.PRED2), 'nu.wre'] <- fullSampleCUNY$DATA.PRED2$nu.wre

brainChartsDF[row.names(fullSampleRU$DATA.PRED2), 'mu.wre'] <- fullSampleRU$DATA.PRED2$mu.wre
brainChartsDF[row.names(fullSampleRU$DATA.PRED2), 'sigma.wre'] <- fullSampleRU$DATA.PRED2$sigma.wre
brainChartsDF[row.names(fullSampleRU$DATA.PRED2), 'nu.wre'] <- fullSampleRU$DATA.PRED2$nu.wre

brainChartsDF[row.names(fullSampleSI$DATA.PRED2), 'mu.wre'] <- fullSampleSI$DATA.PRED2$mu.wre
brainChartsDF[row.names(fullSampleSI$DATA.PRED2), 'sigma.wre'] <- fullSampleSI$DATA.PRED2$sigma.wre
brainChartsDF[row.names(fullSampleSI$DATA.PRED2), 'nu.wre'] <- fullSampleSI$DATA.PRED2$nu.wre

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
fakeSampleCBIC <- calibrateBrainCharts(filter(brainChartsDFFake, study == "HBNCBIC2"), phenotype = phenotype)
fakeSampleCUNY <- calibrateBrainCharts(filter(brainChartsDFFake, study == "HBNCUNY2"), phenotype = phenotype)
fakeSampleSI <- calibrateBrainCharts(filter(brainChartsDFFake, study == "HBNSI2"), phenotype = phenotype)
fakeSampleRU <- calibrateBrainCharts(filter(brainChartsDFFake, study == "HBNRU2"), phenotype = phenotype)

# print random effects for each site

# print("Random effects for each site")
# print(paste0("mu in the fake should offset by ", muNudge))
# print(paste0("sigma in the fake should offset by ", sigmaNudge))
# print("Random effects for CBIC, (real, fake)")
# print("mu")
# print(c(fullSampleCBIC$expanded$mu$ranef[["HBNCBIC"]], fakeSampleCBIC$expanded$mu$ranef[["HBNCBIC2"]]))
# print("sigma")
# print(c(fullSampleCBIC$expanded$sigma$ranef[["HBNCBIC"]], fakeSampleCBIC$expanded$sigma$ranef[["HBNCBIC2"]]))
# 
# print("Random effects for CUNY, (real, fake)")
# print("mu")
# print(c(fullSampleCUNY$expanded$mu$ranef[["HBNCUNY"]], fakeSampleCUNY$expanded$mu$ranef[["HBNCUNY2"]]))
# print("sigma")
# print(c(fullSampleCUNY$expanded$sigma$ranef[["HBNCUNY"]], fakeSampleCUNY$expanded$sigma$ranef[["HBNCUNY2"]]))
# 
# print("Random effects for RU, (real, fake)")
# print("mu")
# print(c(fullSampleRU$expanded$mu$ranef[["HBNRU"]], fakeSampleRU$expanded$mu$ranef[["HBNRU2"]]))
# print("sigma")
# print(c(fullSampleRU$expanded$sigma$ranef[["HBNRU"]], fakeSampleRU$expanded$sigma$ranef[["HBNRU2"]]))
# 
# print("Random effects for SI, (real, fake)")
# print("mu")
# print(c(fullSampleSI$expanded$mu$ranef[["HBNSI"]], fakeSampleSI$expanded$mu$ranef[["HBNSI2"]]))
# print("sigma")
# print(c(fullSampleSI$expanded$sigma$ranef[["HBNSI"]], fakeSampleSI$expanded$sigma$ranef[["HBNSI2"]]))

# now subsample the fake sites

# here, the subFake dataframes represent your "small" samples
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

# run the centile-informed calibration

subFakeCBICOutQuantile <- calibrateBrainChartsIDQuantilePenalty(subFakeCBIC, phenotype = "CT", largeSiteOutput = fullSampleCBIC)
subFakeCUNYOutQuantile <- calibrateBrainChartsIDQuantilePenalty(subFakeCUNY, phenotype = "CT", largeSiteOutput = fullSampleCUNY)
subFakeSIOutQuantile <- calibrateBrainChartsIDQuantilePenalty(subFakeSI, phenotype = "CT", largeSiteOutput = fullSampleSI)
subFakeRUOutQuantile <- calibrateBrainChartsIDQuantilePenalty(subFakeRU, phenotype = "CT", largeSiteOutput = fullSampleRU)

# accuracy evaluation

# run the normal Braincharts calibration on the small sites 
subFakeCBICOut <- calibrateBrainCharts(subFakeCBIC, phenotype = "CT")
subFakeCUNYOut <- calibrateBrainCharts(subFakeCUNY, phenotype = "CT")
subFakeRUOut <- calibrateBrainCharts(subFakeRU, phenotype = "CT")
subFakeSIOut <- calibrateBrainCharts(subFakeSI, phenotype = "CT")

# make plots 

fullSampleCBIC$DATA.PRED2$fitting <- 'large sample (cross-sectional method)'
subFakeCBICOutQuantile$DATA.PRED2$fitting <- 'small site (centile method)'
subFakeCBICOut$DATA.PRED2$fitting <- 'small site (cross-sectional method)'
fakeSampleCBIC$DATA.PRED2$fitting <- 'small site (ground truth)'

fullSampleCBIC$DATA.PRED2$sample <- "large"
fullSampleCBIC$DATA.PRED2
subFakeCBICOutQuantile$DATA.PRED2$sample <- "small"
subFakeCBICOut$DATA.PRED2$sample <- "small"
fakeSampleCBIC$DATA.PRED2$sample <- "small"

T <- bind_rows(fullSampleCBIC$DATA.PRED2, subFakeCBICOutQuantile$DATA.PRED2, subFakeCBICOut$DATA.PRED2, fakeSampleCBIC$DATA.PRED2)

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