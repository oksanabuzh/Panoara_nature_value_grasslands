# Code for Janišová et al. (2025) [Biological Conservation](https://www.sciencedirect.com/journal/biological-conservation)

## Janišová M., Magnes M., Škodová I., Iuga A., Ivașcu A., Ivașcu C.M., Ďuricová V., Tarog A., Kromka M., Buzhdygan O.Y. (2025) Role of livestock and traditional management practices in maintaining high nature value grasslands


# Project Structure

This project is structured as follows:

```md
.gitignore
.Rproj.user/
data/
    Divers_LandUse_Soil_Variables.csv
    Community_composition_DungExperiment.csv
    Community_composition_VegetationPlots.csv
    NMDS_data.csv
    soil_NPK_PCA
analysis/
    0.1_NMDS.R
    0.2_PCA_soil.R
    04_Summary_Statistics.R
    05_GLMMs.R
    05.2_GLMMs_R2_plot.R
    06.1_SEM_SpRich.R
    06.2_SEM_NMDS.R
    06.3_SEM_R2_plot.R
    07_Seed_Exper.R
results/
Panoara.Rproj
```

## R Files

### Folder `analysis`

In this folder, the raw data is prepared for analysis and the analysis is performed


- `1.0_PCA_soil.R`: performs PCA for soil data
- `2.0_NMDS.R`: calculates community composition (NMDS) for each plot, performs the PERMANOVA analysis and ordination plots.
- `3.0_Summary_Statistics.R`: summarizes data and calculates summary statistics and correlations among the measures of plant community (using GLMMs)
- `05_GLMMs.R`: performs the GLMM analysis for plant species richness (field data).
- `06.1_SEM_SpRich.R`: performs the SEM analysis, calculates direct & indirect effects 
for species richness 
- `06.2_SEM_NMDS.R`: performs the SEM analysis, calculates direct & indirect effects
for community composition (NMDS)
- `06.3_SEM_NMDS.R`: 
- `07_Seed_Exper.R`: run GLMMs for the Dung Experiment Data.

## Data files

### Folder `data`

This folder contains the raw data files and the the processed data files.

Raw data files:
- `Community_composition_VegetationPlots.csv`: contains a list of plant species and their cover for each field plot
- `Community_composition_DungExperiment.csv`: contains a list of plant species and their abundances for the seed experiment 
- `Divers_LandUse_Soil_Variables.csv`: This file contains environmental variables and biodiversity measures for all plots

Processed files:
- `NMDS_data.csv`: contains NMDS scores, produced by analysis/`2.0_NMDS.R`
- `soil_NPK_PCA.csv`: contains PCA scores, produced by analysis/`1.0_PCA_soil.R` 



## Results files

### Folder `results`

`PERMANOVA_Field_Data.csv`- results of PERMANOVA analysis for field data, produced by analysis/`2.0_NMDS.R`.

`PERMANOVA_Experiment.csv`- results of PERMANOVA analysis for seed-experiment data, produced by analysis/`2.0_NMDS.R`.

`Sp.occurances_field.csv` - species occurrence (% parcels where species was found).
Produced by analysis/`3.0_Summary_Statistics.R`.



## Other Files

- `.gitignore`: This file specifies intentionally untracked files that Git should ignore.
- `Panoara.Rproj`: This is the R project file for this project.