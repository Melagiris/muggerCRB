#muggerCRB  -  Multi-scale habitat modeling to improve spatial prioritization for mugger crocodile conservation in riverine landscapes

R scripts used for the analysis
1) 01 Bias-Weighted Background Sampling.R   --> Bias surface creation & background point generation 
2) 02 Background Sensitivity Test (Uniform vs Bias-5k vs Bias-10k).R  --> Uniform vs bias-weighted sensitivity test
3) 03 Correlation & variable Selection.R
4) 04 Spatial Block CV vs Naïve CV.R
5) Main file.R --> stepAIC, Model fitting (GLM), Validation, Extracting final best model, Prediction, Response curves (Elith et al., 2005), Systematic Conservation Planning using 'prioritizr' (Hanson et al., 2025), Priority Index
