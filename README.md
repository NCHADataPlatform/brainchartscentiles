# Braincharts Centile informed calibration for repeated measures

Software for improvement of small sample calibration using repeated measures from a large sample.

Please cite the publication ...

# Background

Braincharts [[1]](#1). provides formulae that map predictor variables (age, sex, FreeSurfer version) to curves for distribution parameters that characterise brain IDPs cortical thickness, grey matter volume, cranial volume, cortical surface area, subcortical grey matter volume, cerebral white matter volume, lateral ventricle volumes. The curves are based on data from over 100,000 healthy individuals across the lifespan obrained from many samples. The distribution provided by Braincharts is the Generalized Gamma Distribution, which is a three-parameter distribution with location ($$\mu$$) and shape ($$\sigma$$, $$\nu$$) variables. Calibration of non-training samples involves estimating offset factors for parameters $$\mu$$ and $$\sigma$$ ($$\nu$$ is fixed) that shift the predicted distribtions to match that of the novel sample. The Braincharts' method is known to be unstable for small (n < 100) samples. For a small sample with repeated measures in a large (n > 100) sample, such a sample has correction factors that are considered stable, this software allows researchers to leverage repeated samples to improve stability of the small sample correction factors. 

# Quick start

This software runs as as addon to the Braincharts software, which was written in R. In a shell, obtain the software from GitHub as follows:

```
git clone https://github.com/NCHADataPlatform/brainchartscentiles.git
cd brainchartscentiles
git submodule init .
git submodule update
```

## R dependencies

In R, install the required packages for Braincharts:

```
install.packages(c('gamlss', 'tidyverse'))
```

Then install the dependencies for the example scripts on the github repository:

```
install.packages(c('ggplot2', 'pracma'))
```

The previous two commands only need to be executed once. 

For each R session, execute the following to import the software:

```
source('<install directory>/calibrate_braincharts_centiles.R')
```

This imports all required scripts from the Braincharts software.

# Basic set up

Users are required to have two collections of data prior to starting, predictor and IDP variables. These are as follows:

1. Predictor variables:
    -	*age_days*: age of subject in days
    -	*dx*: {“CN” control, “notCN” patient}, only CN subjects are used for calibration
    -	*sex*: factor column {“M”, “F”}
    - *study*: string representing the sample, must be different than the studies in the training set, which are*: 3R-BRAIN, ABCD, abide1, abide2, ADHD200, ADNI, AIBL, AOBA, AOMIC_ID1000, AOMIC_PIOP1, AOMIC_PIOP2, ARWIBO, BabyConnectome, BGSP, BHRCS, BSNIP, Calgary, CALM, CamCAN, CAMFT, Cornell_C1, Cornell_C2, cVEDA, devCCNP, dHCP, DLBS, EDSD, EMBARC, FemaleASD, FinnBrain, GUSTO, HABS, Harvard_Fetal1, HBN, HCP, HCP_lifespanA, HCP_lifespanD, IBIS, ICBM, IMAGEN, IMAP, IXI, LA5c, MCIC, Narratives, NHGRI, NIH, NIHPD, NSPN, NTB_Yale, OASIS3, Oulu, PCDC, PING, PNC, POND, PREVENTAD, RDB, SALD, SLIM, UCSD, UKB, VETSA, WAYNE
    - *participant*: participant ID, must be preserved for the repeated measures
    - *fs_version*: Freesurfer version used, supported values are {"FS6_T1", "FS6_T1T2", "Custom", "FS53", "Custom_T1T2", "FSInfant"}, ensure you choose a consistent version in your scripts, the particular choice doesn't matter too much, just ensure you always choose the same, Custom_T1T2 is used in the example scripts
1.	IDP variables (one is chosen for analysis), only the analysis IDP is required:
    -	*CT*: mean cortical thickness, both hemispheres, mm
    -	*GMV*: total grey matter volume, mm3
    -	*TCV*: total cranial volume, mm3
    -	*SA*: total cortical surface area, both hemispheres, mm2
    - *sGMV*: subcortical grey matter volume, mm3 
    - *WMV*: cerebral white matter volume mm3
    - *Ventricles*: volume of lateral ventricles, mm3

The predictor variables need to be acuquired by the investigators. The IDP variables can be derived from the Freesurfer stats output files, namely the aseg.stats for volumetric variables and ?h.aparc.stats for surface-based values. We provide scripts for automatic extraction of IDP variables from Freesurfer stats files. A worked example of a full analysis will now be described.

# Worked example

This worked example details how to extract user data from Freesurfer outputs that can be calibrated using the proposed method. To test the "small"-sample calibration procedures addressed by this method specifically, we generate synthetic ("fake") longitudinal data, compare the calibration results with the original (cross-sectional) method and visualize the results. This worked example is also given as the file `freesurfer_example_hbn.R` on the GitHub repo. Given that the sample sizes are "large" (n > 100) this example can perform accuracy evaluation of the centile-method in a manner performed in the paper. Specifically the example:

1. Performs calibration on the full samples (ground truth)
1. Makes fake longitudinal versions of the existing sites to simulate repeated measures
1. Make ground truth calibration on the "fake" full samples
1. Subsample the fake samples to make "small" samples (n < 100) and perform the original cross-sectional and centile informed (proposed method) calibration. We can then compare the original and proposed method calibration accuracy on the fake subsamples.

## Freesurfer data preparation

Assuming you have a Freesurfer directory with completed pipelines for your sample. Use the script `freesurfer_make_stats.sh` in the GitHub repo to extract all subjects' volumetric and surface measures into csv files as follows:

- *stats_aseg.csv*: volumetric stats
- *stats_aparc_lh_thickness.csv*, *stats_aparc_rh_thickness.csv*: cortical surface thickness files
- *stats_aparc_lh_area.csv*, *stats_aparc_rh_area.csv*: cortical surface area files

For the worked example on the repo, these files are in the directory `HBN`.

The idea of presenting the worked example is to demonstrate the steps of the software and, importantly, act as a template that can be modified for users' own data.

## Loading/formatting data in R

Load the stats tables from your Freesurfer SUBJECTS_DIR. This example loads the data in the `HBN` directory of the repository, make sure the working directory of R is the installation directory of the software for this to work. Code:

```
asegDF <- read.csv(file.path('HBN', 'stats_aseg.csv'))
lhCortThicknessDF <- read.csv(file.path('HBN', 'stats_aparc_lh_thickness.csv'))
lhCortSurfaceAreaDF <- read.csv(file.path('HBN', 'stats_aparc_lh_area.csv'))
rhCortThicknessDF <- read.csv(file.path('HBN', 'stats_aparc_lh_thickness.csv'))
rhCortSurfaceAreaDF <- read.csv(file.path('HBN', 'stats_aparc_rh_area.csv'))
```

These lines of code ensure the rows are indexed by participant IDs:

```
row.names(asegDF) <- asegDF$Measure.volume
row.names(lhCortSurfaceAreaDF) <- lhCortSurfaceAreaDF$lh.aparc.area
row.names(rhCortSurfaceAreaDF) <- lhCortSurfaceAreaDF$rh.aparc.area
row.names(lhCortThicknessDF) <- lhCortThicknessDF$lh.aparc.area
row.names(rhCortThicknessDF) <- lhCortThicknessDF$rh.aparc.area
```

This line creates the main dataframe (described above), with the row names indexed by subject ID:

```
brainChartsDF <- data.frame(row.names=asegDF$Measure.volume)
```

The following lines of code extract the relevant variables from the stats and creates the Braincharts' IDP values:

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

These lines of code create other required columns

```
brainChartsDF$fs_version <- 'Custom_T1T2'
brainChartsDF$run <- 1
brainChartsDF$participant <- row.names(brainChartsDF)
brainChartsDF$country <- 'Australia'
```

This line sets the diagnosis column for each subject. Braincharts requires labels *CN* for control subjects (for calibration), and *notCN* for patient or non-control subjects (not for calibration). Here, all subjects are controls so we set them as such:

```
brainChartsDF$dx <- 'CN'
```

Users can have non-control subjects in the dataframe. These subjects' data will be ignored for calibration but their estimated parameter distrubution values and centiles will be available in the outputs.

The following line sets the *study* column:

```
brainChartsDF$study <- 'MyStudy'
```

Here, *study* means sample. You may have multiple studies in the dataframe; each study will be calibrated independently. The study can be named anything except the study names of the Braincharts' training set, which are as follows: 3R-BRAIN, ABCD, abide1, abide2, ADHD200, ADNI, AIBL, AOBA, AOMIC_ID1000, AOMIC_PIOP1, AOMIC_PIOP2, ARWIBO, BabyConnectome, BGSP, BHRCS, BSNIP, Calgary, CALM, CamCAN, CAMFT, Cornell_C1, Cornell_C2, cVEDA, devCCNP, dHCP, DLBS, EDSD, EMBARC, FemaleASD, FinnBrain, GUSTO, HABS, Harvard_Fetal1, HBN, HCP, HCP_lifespanA, HCP_lifespanD, IBIS, ICBM, IMAGEN, IMAP, IXI, LA5c, MCIC, Narratives, NHGRI, NIH, NIHPD, NSPN, NTB_Yale, OASIS3, Oulu, PCDC, PING, PNC, POND, PREVENTAD,  RDB, SALD, SLIM, UCSD, UKB, VETSA, WAYNE. 

In our example, our study values will be: HBNCBIC, HBNCUNY, HBNSI, HBNRU for the sites.

Loading in the demographics will be specific to how that data is stored for your sample. The code for this will not be shown here. Refer to the source code if interested.

This line sets a phenotype variable, which indicates the IDP we want to analyse. Here, we use *CT* for *cortical thickness*:

```
phenotype <- "CT"
```

## Full sample calibrations

First, do the full-site calibrations using the original cross-sectional method on the full samples:

```
fullSampleCBIC <- calibrateBrainCharts(filter(brainChartsDF, study == "HBNCBIC"), phenotype = phenotype)
fullSampleCUNY <- calibrateBrainCharts(filter(brainChartsDF, study == "HBNCUNY"), phenotype = phenotype)
fullSampleSI <- calibrateBrainCharts(filter(brainChartsDF, study == "HBNSI"), phenotype = phenotype)
fullSampleRU <- calibrateBrainCharts(filter(brainChartsDF, study == "HBNRU"), phenotype = phenotype)
```

In a real example, these results would represent results from your "large" sample.

In this example, we generate fake data for each site, adding a year to each subject's age, simulating a second timepoint, and altering the distribution parameters. For each site, make a longitudinal second session with altered site parameters, generate new cortical thickness measures. The code won't be shown here but the dataframe brainChartsDFFake contains altered *age_days* and *CT* colums. As before, we calibrate the full fake samples, but the results of that are only used for accuracy evaluation and wouldn't be available in practice.

## Results description

Here, the results data structure will be described, this is purely for reference. The dataframe *fullSampleCBIC* is the output

- *fullSampleCBIC$expanded*: The model parameters 
    - *fullSampleCBIC$expanded$mu$ranef*: The $\mu$ site-effect estimates. The user-provided study name will be the estimated effect for the sample, i.e. fullSampleCBIC$expanded$mu$ranef[['HBNCBIC']]
    - *fullSampleCBIC$expanded$sigma$ranef*: The $\mu$ site-effect estimates. The user-provided study name will be the estimated effect for the sample, i.e. fullSampleCBIC$expanded$mu$ranef[['HBNCBIC']]
- *fullSampleCBIC$DATA.PRED*: The data with predicted values from the model *before* calibration
- *fullSampleCBIC$DATA.PRED2*: The data with predicted values for parameters *after* calibration
    - *fullSampleCBIC$DATA.PRED2$meanCT2transformed*



## Small sample calibrations

In the code, the data frame *subFakeCBIC* is a random subsample of 50 subjects from the synthetic longitudinal timepoint for the CBIC sample. In a real example, this dataframe would contain the data from your "small" sample.

We calibrate using the centile-informed repeated measures method, using the large site calibration output previously computed as follows:

```
subFakeCBICOutQuantile <- calibrateBrainChartsIDQuantilePenalty(subFakeCBIC, phenotype = "CT", largeSiteOutput = fullSampleCBIC)
```

## Accuracy evaluation

This part is only relevant if you want to evaluate the accuracy of the cross-sectional (original) and centile-informed method on the subsample calibration. Firstly, run the Braincharts calibration on the full fake sample:

```
subFakeCBICOut <- calibrateBrainCharts(subFakeCBIC, phenotype = "CT")

```

Then, make plots to show the predicted and ground truth centiles

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

# References

<a id="1">[1]</a> 
Bethlehem, R.A.I., Seidlitz, J., White, S.R. et al. Brain charts for the human lifespan. Nature 604, 525–533 (2022). https://doi.org/10.1038/s41586-022-04554-y
