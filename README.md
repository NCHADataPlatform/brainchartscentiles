# Braincharts Centile informed calibration for repeated measures

Software for improvement of small sample calibration using repeated measures from a large sample.

Please cite the publication ...

# Background

Braincharts [[1]](#1) provides formulae that map predictor variables (age, sex, FreeSurfer version) to curves for distribution parameters that characterise brain IDPs cortical thickness, grey matter volume, cranial volume, cortical surface area, subcortical grey matter volume, cerebral white matter volume, lateral ventricle volumes. The curves are based on data from over 100,000 healthy individuals across the lifespan obrained from many samples. The distribution provided by Braincharts is the Generalized Gamma Distribution, which is a three-parameter distribution with location ($$\mu$$) and shape ($$\sigma$$, $$\nu$$) variables. Calibration of non-training samples involves estimating offset factors for parameters $$\mu$$ and $$\sigma$$ ($$\nu$$ is fixed) that shift the predicted distribtions to match that of the novel sample. The Braincharts' method is known to be unstable for small (n < 100) samples. For a small sample with repeated measures in a large (n > 100) sample, such a sample has correction factors that are considered stable, this software allows researchers to leverage repeated samples to improve stability of the small sample correction factors. 

# Installation

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

# Tutorials

Please find tutorials on how to use the software at these links as follows:

- [Calibration example](https://nchadataplatform.github.io/brainchartscentiles/calibration_example.html): Data set up and calibration for a small site calibration using the Healthy Brain Network dataset [[2]](#2).
- [calibration accuracy evaluation](https://nchadataplatform.github.io/brainchartscentiles/calibration_accuracy.html): An evaluation of the accuracy of this centile-informed method on a small sample using the public Healthy Brain Network dataset [[2]](#2).

# References

<a id="1">[1]</a> 
Bethlehem, R.A.I., Seidlitz, J., White, S.R. et al. Brain charts for the human lifespan. Nature 604, 525–533 (2022). https://doi.org/10.1038/s41586-022-04554-y
<a id="2">[2]</a> Alexander, L. et al. An open resource for transdiagnostic research in pediatric mental health and learning disorders. Scientific Data 4, 170181 (2017).