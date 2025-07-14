# FSJhelperEffect
This repository stores the raw data and scripts used in the submitted manuscript for Summers et al. 2025, studying the impact of helper sex and relatedness on the fitness effects of helping

All analysis scripts are located in the \scripts folder. These scripts include:
- HelperModels_20250714.R : a script that takes in the raw anonymized data (FitnessDataSet_20250714.rdata) to create regression models linking breeder survival, offspring survival, and nestling production to the number, sex, and relatedness of the helpers present on a territory. All models are fit using the glmmTMB package. This script includes tests using the DHARMa package to check model fit, overdispersion, and outliers. This script outputs HelperModels_20250714.rdata.

- allModelFigures20250714.R : a script that takes in the fit models (HelperModels_20250714.rdata) and raw anonymized data (FitnessDataSet_20250714.rdataa) to create tables and figures.


The raw anonymized data in FitnessDataSet_20250714.rdata includes the following data frames:
- Offspring.input: Offspring survival data frame with individual ID (USFWS), natal year (Year), natal nest ID (NatalNest), territory ID (Terr), territory-year joint ID (TerrYr), sex, julian hatch date (HatchDate), individual pedigree-based inbreeding coefficient (pedF), # of years parents have been paired (pairYear), survival outcome to next breeding season (surv), number of helpers (num.helpers), ratio of female to male helpers (sex.ratio), ratio of related...
