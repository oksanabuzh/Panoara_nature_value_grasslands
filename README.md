# Code for Janišová et al. (2024) [Biological Conservation](https://www.sciencedirect.com/journal/biological-conservation)

## Janišová M., Magnes M., Škodová I., Iuga A., Ivașcu A., Ivașcu C.M., Ďuricová V., Tarog A., Kromka M., Buzhdygan O.Y. (2024) 
## Role of livestock in maintaining high nature value grasslands

Special issue in Biological Conservation ["Conservation of Palaearctic steppes and semi-natural grasslands: challenges and solutions"](https://www.sciencedirect.com/journal/biological-conservation/about/call-for-papers#conservation-of-palaearctic-steppes-and-semi-natural-grasslands-challenges-and-solutions)



# Project Structure

This project is structured as follows:

```md
.gitignore
.Rproj.user/
data/
    Variables_selected.csv
    Community_composition_DungExperiment.csv
    Community_composition_VegetationPlots.csv
    Diversity_&_NMDS_data.csv
    Panoara_Dat.csv
analysis/
    01_Calculate_diversity_&_composition.R
    02_PCA_variables.R
    03_Prepare_data.R
    04_Summary_Statistics.R
    05_GLMMs.R
    06_SEM.R
    old/
results/
Panoara.Rproj
```

## R Files

### Folder `analysis`

In this folder, the raw data is prepared for analysis and the analysis is performed

- `01_Calculate_diversity_&_composition.R`: calculates evenness and community composition (NMDS) for each plot, performs the PERMANOVA analysis and ordination plots.
- `02_PCA_variables.R`: performs PCA for soil data
- `03_Prepare_data.R`: prepares the data for analysis. It reads in the raw data, renames and filters the variables, join the data sets.
- `04_Summary_Statistics.R`: summarises data and calculates summary statistics.
- `05_GLMMs.R`: performs the GLMM analysis
- `06_SEM.R`: performs the SEM analysis
    

## Data files

### Folder `data`

This folder contains the raw data files and the the processed data files.

Raw data files:
- `Variables_selected`: This is the raw data file used in the `03_Prepare_data.R` script. It contains environmental variables and data from vegetation surveys and seed experiment.
- `Community_composition_VegetationPlots.csv`: contains a list of plant species and their cover for each field plot. 
- `Community_composition_DungExperiment.csv`: contains a list of plant species and their abundances for the seed experiment. 
- `headers.csv`: This file contains environmental variables for all plots.

Processed data files:
- `Panoara_Dat.csv`: This file is created by the `03_Prepare_data.R` script. It contains clean and joined dataset of the environmental variables, biodiversity measures, the NMDS and PCA scores.


## Other Files

- `.gitignore`: This file specifies intentionally untracked files that Git should ignore.
- `Panoara.Rproj`: This is the R project file for this project.



