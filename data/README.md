# Folder "data"

## About

Contains data that are used for analysis 

#### Structure

| Document                                   | Description                              |
| ------------------------------------------ |----------------------------------------- |
|`Community_composition_VegetationPlots.csv` | Community matrix for plant species composition at the field farm parcels     |
|`Community_composition_DungExperiment.csv`  | Community matrix for plant species composition germinated at  seed experiment     |
|`Diversity_&_NMDS_data.csv`                 | Diversity measures and community composition NMDS scores from the [diversity calculations and NMDS analysis](analysis/0.1_Calculate_diversity_&_composition.R)    |
|`soil_PC.csv`                               | scores from the [PCA analysis](analysis/0.2_PCA_soil.R)  |     
|`LandUse_soil_variables.csv`                | land use and soil variables | 



## Metadata

### -> `Community_composition_VegetationPlots.csv` 
Contains plant species composition at the field farm parcels

#### Variables:			
|Short name	| type |	Long name	| Description |
| ----------|------|------------| ------------|
| species	  |character | species name	    | species name                             |
| layer     |character | vegetation group	|VP - vascular, B - bryophyte, L - lichen	 |
| Farm parcels (as matrix columns)   | numeric  | 	plant cover in each parcel |cells represent cover (%) of each species at each parcel	 |


### ->  `Community_composition_DungExperiment.csv`
Contains plant species composition germinated at the seed experiment

#### Variables:			
|Short name	| type |	Long name	| Description |
| ----------|------|------------| ------------|
| species	  |character | species name	    | species name                             |
| layer     |character | vegetation group	|VP - vascular, B - bryophyte, L - lichen	 |
| Farm parcels (as matrix columns)   | numeric  | 	plant cover in each parcel |cells represent the cumulative cover (%) of each species at each parcel	 |



### ->  `Diversity_&_NMDS_data.csv`
Contains diversity measures and community composition NMDS scores from the [diversity calculations and NMDS analysis](analysis/0.1_Calculate_diversity_&_composition.R) 

#### Variables:			
|Short name	| type |	Long name	| Description |
| ----------|------|------------| ------------|
| Parcel_name	     | character  | parcel name	|parcel and farm ID |              
| CoverVP_field    | numeric    | plant cover  |cumulative cover (%) of plant community at each parcel, measured at the field	|
| SR_VP_field      | numeric    | species richness |species richness of plant community at each parcel, measured at the field	|
| EvennessVP_field | numeric    | evenness |evenness (inverse Simpson) of plant community at each parcel, measured at the field |
| ShannonVP_field  | numeric    | Shannon diversity |Shannon diversity measure of plant community at each parcel, measured at the field	|
| abundance_Exper  | numeric    | abundance | abundance of plant seedlings  germinated at  seed experiment	|
| SR_Exper  | numeric    | species richness | species richness of plant seedlings  germinated at  seed experiment	|
| Evenness_Exper  | numeric    | evenness | evenness (inverse Simpson) of plant seedlings  germinated at  seed experiment	|
| Shannon_Exper  | numeric    | Shannon diversity | Shannon diversity of plant seedlings  germinated at  seed experiment	|
| NMDS1_VP_field  | numeric    | NMDS1 score | scores for the 1st NMDS axis from the NMDS analysis of the species composition of plant community at each parcel, measured at the field	|
| NMDS2_VP_field  | numeric    | NMDS2 score | scores for the 2nd NMDS axis from the NMDS analysis of the species composition of plant community at each parcel, measured at the field	|
| NMDS1_exper  | numeric    | NMDS1 score | scores for the 1st NMDS axis from the NMDS analysis of the species composition of plant seedlings  germinated at seed experiment	|
| NMDS2_exper  | numeric    | NMDS2 score | scores for the 2nd NMDS axis from the NMDS analysis of the species composition of plant seedlings  germinated at seed experiment	|



### ->  `LandUse_soil_variables.csv`
Contains land use and soil variables

#### Variables:			
|Short name	| type |	Long name	| Description |
| ----------|------|------------| ------------|
| Parcel_name	     | character  | parcel name	|parcel and farm ID |     

                |  | 


### -> `alpha_beta_gamma_community_variabl.csv`
Combines all diversity measures and plant cover for each scale:
- alpha diversity and cover measures (SR, ENSPIE, and cover) include doubled 10 m2 plots
- gamma diversity and cover measures (SR, ENSPIE, and cover) include 100 m2 plots (i.e. the sample size is half of what we have for the alpha diversity)
- beta diversity measures (SR and ENSPIE) are calculated as gamma/alpha

SR - species richness

ENSPIE - evenness measure calculated as inverse Simpson using species cover

cover - cumulative cover of plant community

#### Variables:			
|Short name	| type |	Long name	| Description |
| ----------|------|------------| ------------|
| dataset	  |numeric | Dataset ID	| Dataset ID, sampled by different teams and years |
|series	| character	| Series ID	| 100 m2 plot that includes two 10 m2 plots that are nested within it |
|subplot	| character	| Subplot 	| one of two corners (i.e. 10 m 2 plots) nested within the 100 m2 plot (called series): NW - north west corner; SE - south east corner |
|plotID	| character	| PlotID	| Plot ID, combines information of both series and corner. It is a unique identification ID for the study plots|
|scale	| integer |	Spatial scale | Grain size of the sampled plots: 10 is 10 m2 plots; 100 is 100 m2 plots, NA - nonaplicable (for beta diversity)|
|type	| character	| Scale type as used in the paper	| Name of spatial scale: alpha - diversity and cover measures include doubled 10 m2 plots; gamma - diversity and cover measures include 100 m2 plots (i.e. the sample size is half of what we have for the alpha diversity); beta diversity measures are calculated as gamma/alpha |
|metric	| character	| Measure of plant community |	Measure of biodiversity and cover of plant community: SR - species richness; ENSPIE - evenness measure calculated as inverse Simpson using species cover; cover - cumulative cover of plant community |
|value	| numeric	| Value for the respective measure of plant community	| Value for the respective measure of plant community|


### -> `climate_PCA.csv`
Contains scores for the compound climate variable, derived from the [PCA analysis](analysis/PCA_climate.R), see `analysis/Readme.md` for details.



### -> `aggregation.csv`
Contains beta.BRAY.BAL - the proxy of intraspecific aggregation  for each plot. Spatial intraspecific aggregation was estimated by comparing dissimilarity in species covers between the two corners (i.e., two 10-m2 plots) within each 100-m2 plot. For this, we calculated  beta.BRAY.BAL - the balanced variation component of Bray–Curtis dissimilarity in species cover using ‘betapart’ package in R (Baselga & Orme, 2012). This measure is independent of total community abundance (total plant cover in our study) and measures the balanced variation in species abundance between two quadrats, i.e. when  cover increases for some species and decreases for others, maintaining similar total cover across quadrats, including also species turnover, where abundance of one species is replaced by other species (Baselga, 2017). Higher dissimilarity in covers of taxa between the two 10-m2 corners within the same 100-m2 plot implies higher intraspecific aggregation.

Baselga, A. (2017). Partitioning abundance-based multiple-site dissimilarity into components: balanced variation in abundance and abundance gradients. Methods in Ecology and Evolution, 8(7), 799–808. https://doi.org/10.1111/2041-210X.12693
Baselga, A., & Orme, C. D. L. (2012). Betapart: An R package for the study of beta diversity. Methods in Ecology and Evolution, 3(5), 808–812. https://doi.org/10.1111/j.2041-210X.2012.00224.x


#### Variables:			
|Short name	| type |	Long name	| Description |
| ----------|------|------------| ------------|
|series	| character	| Series ID	| 100 m2 plot that includes two 10 m2 plots that are nested within it |
|beta.BRAY.BAL	| numeric	| aggregation 	| Proxy of intraspecific aggregation of plant community |
