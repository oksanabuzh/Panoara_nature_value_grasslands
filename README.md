# Code for Janišová et al. (2025) Biological Conservation

## Janišová M., Škodová I., Magnes M., Iuga A., Biro A.-S., Ivașcu C.M., Ďuricová V., Buzhdygang O. 
## Role of livestock and traditional management practices in maintaining high nature value grasslands [Biological Conservation](https://www.sciencedirect.com/journal/biological-conservation), 2025 
### https://doi.org/10.1016/j.biocon.2025.111301

* correspondence to oksana.buzh@fu-berlin.de

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
    1.0_PCA_soil.R
    2.0_NMDS.R
    3.0_Summary_Statistics.R
    4.1_GLMMs_Field.R
    4.2_GLMMs_Seed_Exper
    4.3_GLMMs_R2_plot
   results/
Panoara.Rproj
```

## R Files

### Folder [`analysis`](analysis)

In this folder, the raw data is prepared for analysis and the analysis is performed

- `1.0_PCA_soil.R`: performs PCA for soil data
- `2.0_NMDS.R`: calculates community composition (NMDS) for each plot, performs the PERMANOVA analysis and ordination plots.
- `3.0_Summary_Statistics.R`: summarizes data and calculates summary statistics and correlations among the measures of plant community (using GLMMs)
- `4.1_GLMMs_Field.R`: performs the GLMM analysis for plant species richness for the field data.
- `4.2_GLMMs_Seed_Exper.R`: run GLMMs for the Dung Experiment Data.
- `4.3_GLMMs_R2_plot.R`: create R2 plot for the GLMMs (Fig. 3).
- `5.1_SEM_SpRich.R`: performs the SEM analysis for species richness 
- `5.2_SEM_NMDS.R`: performs the SEM analysis for community composition (NMDS)
- `5.3_SEM_R2_plot.R`: creates R2 plot for the SEM (Fig. 5E)

## Data files

### Folder [`data`](data)

This folder contains the raw data files and the the processed data files. For details see  [data/README](data) 

Raw data files:
- `Community_composition_VegetationPlots.csv`: contains a list of plant species and their cover for each field plot
- `Community_composition_DungExperiment.csv`: contains a list of plant species and their abundances for the seed experiment 
- `Divers_LandUse_Soil_Variables.csv`: This file contains environmental variables and biodiversity measures for all plots

Processed files:
- `NMDS_data.csv`: contains NMDS scores, produced by analysis/`2.0_NMDS.R`
- `soil_NPK_PCA.csv`: contains PCA scores, produced by analysis/`1.0_PCA_soil.R` 


## Results files

### Folder [`results`](results)

- `PERMANOVA_Field_Data.csv`- results of PERMANOVA analysis for field data, produced by analysis/`2.0_NMDS.R`.
- `PERMANOVA_Experiment.csv`- results of PERMANOVA analysis for seed-experiment data, produced by analysis/`2.0_NMDS.R`.
- `Sp.occurances_field.csv` - species occurrence (% parcels where species was found).
Produced by analysis/`3.0_Summary_Statistics.R`.
- `Tables_Anova.csv`: contains GLMMs results, produced by analysis/`4.1_GLMMs_Field.R` 
- `Tables_Model_R2.csv`: contains model R2 for the GLMMs, produced by analysis/`4.1_GLMMs_Field.R` 
- `Tables_R2_partial.csv`: contains partial-R2 for fixed effects tested by the GLMMs, produced by analysis/`4.1_GLMMs_Field.R` 
- `SEM_SR_R2_partial.csv`: contains partial-R2 for the fixed effects for species richness, 
tested in the SEM, produced by analysis/`5.1_SEM_SpRich.R`
- `SEM_NMDS_R2_partial.csv`: contains partial-R2 for the fixed effects for species composition, 
tested in the SEM, produced by analysis/`5.2_SEM_NMDS.R`
- `SEM_SR_coefs.csv`: the estimates from the SEM, produced by analysis/`5.1_SEM_SpRich.R`
- `SEM_SR_coefs_Indirect.csv`: indirect effects from the SEM, produced by analysis/`5.1_SEM_SpRich.R`
- `SEM_NMDS_coefs.csv`: the estimates from the SEM, produced by analysis/`5.2_SEM_NMDS.R`
- `SEM_NMDS_coefs_Indirect.csv`: indirect effects from the SEM, produced by analysis/`5.2_SEM_NMDS.R`

## Other Files

- `.gitignore`: This file specifies intentionally untracked files that Git should ignore.
- `Panoara.Rproj`: This is the R project file for this project.
