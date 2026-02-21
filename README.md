# muggerCRB

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18467553.svg)](https://doi.org/10.5281/zenodo.18467553)

**Multi-scale habitat modeling to improve spatial prioritization for mugger crocodile conservation in riverine landscapes**

This repository contains the code used for the analyses presented in the manuscript:
> *Multi-scale habitat modeling to improve spatial prioritization for mugger crocodile conservation in riverine landscapes*
> 
> Journal: Ecological Applications

The repository implements a multi-scale habitat modeling workflow for evaluating habitat suitability of mugger crocodiles (*Crocodylus palustris*) in Cauvery River Basin.

---

## Repository contents

The repository is organized as follows:

- `01 Bias-Weighted Background Sampling.R`  
  Bias surface creation and background point generation

- `02 Background Sensitivity Test (Uniform vs Bias-weighted).R`  
  Sensitivity analysis of background sampling approaches

- `03 Correlation & Variable Selection.R`  
  Collinearity assessment and predictor selection

- `04 Spatial Block CV vs Naïve CV.R`  
  Comparison of spatial block cross-validation and naïve cross-validation

- `Main File.R`  
  Main script for model fitting, validation, prediction, response curves, and spatial prioritization

---

## Reproducibility

The scripts are intended to reproduce the analytical workflow and modeling framework described in the manuscript.  
Paths to data files may need to be modified by users to match their local directory structure.

---

## Software and versions

Analyses were conducted using:

- **R** (version ≥ 4.2)
- **QGIS** (version 3.40.5–Bratislava) for spatial data processing and visualisation
