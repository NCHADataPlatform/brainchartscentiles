# bccentiles

Software for improvement of small sample calibration using repeated measures from a large sample.

Please cite the publication ...

# Background

Braincharts [1] provides formulae that map predictor variables (age, sex, FreeSurfer version) to curves for distribution parameters that characterise brain IDPs cortical thickness, grey matter volume, cranial volume, cortical surface area, subcortical grey matter volume, cerebral white matter volume, lateral ventricle volumes. The curves are based on data from over 100,000 healthy individuals across the lifespan obrained from many samples. The distribution provided by Braincharts is the Generalized Gamma Distribution, which is a three-parameter distribution with location (mu) and shape (sigma, nu) variables. Calibration of non-training samples involves estimating offset factors for parameters mu and sigma that shift the predicted distribtions to match that of the novel sample. The Braincharts' method is known to be unstable for small (n < 100) samples. For a small sample with repeated measures in a large (n > 100) sample, such a sample has correction factors that are considered stable, this software allows researchers to leverage repeated samples to improve stability of the small sample correction factors. 

# Installation

This software runs as as addon to the Braincharts software, which was written in R. Obtain the software from GitHub as follows:

```
git clone https://github.com/NCHADataPlatform/brainchartscentiles.git
cd brainchartscentiles
git submodule init .
```

## R dependencies

Start an R session. Install the required packages for Braincharts:

```
install.packages(c('gamlss', 'tidyverse'))
```

Dependencies for the example scripts on the github repository:

```
install.packages(c('ggplot2', 'pracma'))
```

The previous two commands only need to be executed once. 

For each R session, execute the following to import the software:

```
source('calibrate_braincharts_centiles.R')
```

This imports all required scripts from the Braincharts software.

# Basic set up

Create a dataframe with the following columns:

- Predictor variables:
  -	*age_days*: age of subject in days
  -	*dx*: {“CN” control, “notCN” patient}, only CN subjects are used for calibration
  -	*sex*: factor column {“M”, “F”}
  - *study*: string representing the sample, must be different than the studies in the training set, which are*: 3R-BRAIN, ABCD, abide1, abide2, ADHD200, ADNI, AIBL, AOBA, AOMIC_ID1000, AOMIC_PIOP1, AOMIC_PIOP2, ARWIBO, BabyConnectome, BGSP, BHRCS, BSNIP, Calgary, CALM, CamCAN, CAMFT, Cornell_C1, Cornell_C2, cVEDA, devCCNP, dHCP, DLBS, EDSD, EMBARC, FemaleASD, FinnBrain, GUSTO, HABS, Harvard_Fetal1, HBN, HCP, HCP_lifespanA, HCP_lifespanD, IBIS, ICBM, IMAGEN, IMAP, IXI, LA5c, MCIC, Narratives, NHGRI, NIH, NIHPD, NSPN, NTB_Yale, OASIS3, Oulu, PCDC, PING, PNC, POND, PREVENTAD, RDB, SALD, SLIM, UCSD, UKB, VETSA, WAYNE
  - *participant*: participant ID, must be preserved for the repeated measures
  - *fs_version*: Freesurfer version used, supported values are {"FS6_T1", "FS6_T1T2", "Custom", "FS53", "Custom_T1T2", "FSInfant"}, ensure you choose a consistent version in your scripts, the particular choice doesn't matter too much, just ensure you always choose the same, Custom_T1T2 is used in the example scripts
-	IDP variables (one is chosen for analysis), only the analysis IDP is required:
    -	*CT*: mean cortical thickness, both hemispheres, mm
    -	*GMV*: total grey matter volume, mm3
    -	*TCV*: total cranial volume, mm3
    -	*SA*: total cortical surface area, both hemispheres, mm2
    - *sGMV*: subcortical grey matter volume, mm3 
    - *WMV*: cerebral white matter volume mm3
    - *Ventricles*: volume of lateral ventricles, mm3

# Worked example

The file freesurfer_example_hbn.R on the GitHub repo will now be described. This file can be used as a template for extracting user data from Freesurfer outputs.

Load the stats tables from your Freesurfer SUBJECTS_DIR, replace 'HBN', with yours:

```
asegDF <- read.csv(file.path('HBN', 'stats_aseg.csv'))
lhCortThicknessDF <- read.csv(file.path('HBN', 'stats_aparc_lh_thickness.csv'))
lhCortSurfaceAreaDF <- read.csv(file.path('HBN', 'stats_aparc_lh_area.csv'))
rhCortThicknessDF <- read.csv(file.path('HBN', 'stats_aparc_lh_thickness.csv'))
rhCortSurfaceAreaDF <- read.csv(file.path('HBN', 'stats_aparc_rh_area.csv'))
```

These csv files are created with the script `freesurfer_make_stats.sh` in th repository.

Set the row names to participant IDs

```
row.names(asegDF) <- asegDF$Measure.volume
row.names(lhCortSurfaceAreaDF) <- lhCortSurfaceAreaDF$lh.aparc.area
row.names(rhCortSurfaceAreaDF) <- lhCortSurfaceAreaDF$rh.aparc.area
row.names(lhCortThicknessDF) <- lhCortThicknessDF$lh.aparc.area
row.names(rhCortThicknessDF) <- lhCortThicknessDF$rh.aparc.area
```

Make the dataframe as described above. Note, this dataframe contains data from multiple samples.

```
brainChartsDF <- data.frame(row.names=asegDF$Measure.volume)
```

Make the IDP columns

```
brainChartsDF$WMV <- asegDF$CerebralWhiteMatterVol
brainChartsDF$sGMV <- asegDF$SubCortGrayVol
brainChartsDF$GWMV <- asegDF$CortexVol
brainChartsDF$etiv <- asegDF$eTIV
brainChartsDF$Ventricles <- asegDF$VentricleChoroidVol
brainChartsDF$TCV <- asegDF$SupraTentorialVolNotVent
brainChartsDF$SA <- lhCortSurfaceAreaDF$lh_WhiteSurfArea_area + rhCortSurfaceAreaDF$rh_WhiteSurfArea_area
brainChartsDF$CT <- (lhCortThicknessDF$lh_MeanThickness_thickness * lhCortSurfaceAreaDF$lh_WhiteSurfArea_area + lhCortThicknessDF$lh_MeanThickness_thickness * rhCortSurfaceAreaDF$rh_WhiteSurfArea_area) / brainChartsDF$SA
brainChartsDF$BrainSegVol.to.eTIV <- asegDF$BrainSegVol.to.eTIV
```

Put in other columns, leave as is

```
brainChartsDF$fs_version <- 'Custom_T1T2'
brainChartsDF$run <- 1
brainChartsDF$participant <- row.names(brainChartsDF)
brainChartsDF$country <- 'Australia'
```

Put in diagnosis, CN for control (for calibration), notCN for patient (not for calibration). Here, all subjects are controls.

```
brainChartsDF$dx <- 'CN'
```

Put in study, must be different than the studies in the training set, which are
 3R-BRAIN, ABCD, abide1, abide2, ADHD200, ADNI, AIBL, AOBA, AOMIC_ID1000,
 AOMIC_PIOP1, AOMIC_PIOP2, ARWIBO, BabyConnectome, BGSP, BHRCS, BSNIP,
 Calgary, CALM, CamCAN, CAMFT, Cornell_C1, Cornell_C2, cVEDA, devCCNP, dHCP,
 DLBS, EDSD, EMBARC, FemaleASD, FinnBrain, GUSTO, HABS, Harvard_Fetal1, HBN, HCP,
 HCP_lifespanA, HCP_lifespanD, IBIS, ICBM, IMAGEN, IMAP, IXI, LA5c, MCIC, Narratives,
 NHGRI, NIH, NIHPD, NSPN, NTB_Yale, OASIS3, Oulu, PCDC, PING, PNC, POND, PREVENTAD,
 RDB, SALD, SLIM, UCSD, UKB, VETSA, WAYNE

````
brainChartsDF$study <- 'MyStudy'
```

We will change this variable per site later. I will use the study names HBNCBIC, HBNCUNY, HBNSI, HBNRU for the sites.

Loading in the demographics will be specific to how that data is stored for your sample. The code for this will not be shown here. Refer to the source code if interested.

choose cortical thickness as phenotype

```
phenotype <- "CT"
```

Do the full sample calibrations per site

```
fullSampleCBIC <- calibrateBrainCharts(filter(brainChartsDF, study == "HBNCBIC"), phenotype = phenotype)
fullSampleCUNY <- calibrateBrainCharts(filter(brainChartsDF, study == "HBNCUNY"), phenotype = phenotype)
fullSampleSI <- calibrateBrainCharts(filter(brainChartsDF, study == "HBNSI"), phenotype = phenotype)
fullSampleRU <- calibrateBrainCharts(filter(brainChartsDF, study == "HBNRU"), phenotype = phenotype)
```

In this example, we generate fake data for each site, adding a year to each subject's age, simulating a second timepoint, and altering the distribution parameters. This would not be done for your sample but it is necessary here.
For each site, make a longitudinal second session with altered site parameters, generate new cortical thickness measures. The code won't be shown here but the dataframe brainChartsDFFake is made 

Calibrate the full fake samples

````
fakeSampleCBIC <- calibrateBrainCharts(filter(brainChartsDFFake, study == "HBNCBIC2"), phenotype = phenotype)
fakeSampleCUNY <- calibrateBrainCharts(filter(brainChartsDFFake, study == "HBNCUNY2"), phenotype = phenotype)
fakeSampleSI <- calibrateBrainCharts(filter(brainChartsDFFake, study == "HBNSI2"), phenotype = phenotype)
fakeSampleRU <- calibrateBrainCharts(filter(brainChartsDFFake, study == "HBNRU2"), phenotype = phenotype)
```

Note, this "full" version of the small sample would normally be unavailable, but is used here for "ground truth". To generate a "small" sample, we simply take 50 random subjects from the fake small site:

```
n <- 50
subFakeCBIC <- filter(brainChartsDFFake, study == "HBNCBIC2")
P <- randperm(1:nrow(subFakeCBIC))
subFakeCBIC <- subFakeCBIC[P[1:n],]
```

So that subFakeCBIC would be the measured small sample you wish to calibrate

Run the normal Braincharts calibration and the Quantile-informated calibration, using the full sample output as the large site reference

```
subFakeCBICOut <- calibrateBrainCharts(subFakeCBIC, phenotype = "CT")
subFakeCBICOutQuantile <- calibrateBrainChartsIDQuantilePenalty(subFakeCBIC, phenotype = "CT", largeSiteOutput = fullSampleCBIC)
```

Make plots 

```
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
```

This plot shows the predicted distribution quantiles. These should match the data.
