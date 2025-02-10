# Folder "data"

## About

Contains data that are used for analysis 

#### Structure

| Document                                   | Description                              |
| ------------------------------------------ |----------------------------------------- |
|`Community_composition_VegetationPlots.csv` | Community matrix for plant species composition at the field farm parcels     |
|`Community_composition_DungExperiment.csv`  | Community matrix for plant species composition germinated at  seed experiment     |
|`NMDS_data.csv`                 | Community composition NMDS scores from the [NMDS analysis](2.0_NMDS.R)    |
|`soil_NPK_PCA.csv`                               | scores from the [PCA analysis](Panoara_nature_value_grasslands/analysis/1.0_PCA_soil.R)  |     
|`Divers_LandUse_Soil_Variables.csv`                | biodiversity, land use and soil variables | 



## Metadata

### -> `Community_composition_VegetationPlots.csv` 
Contains plant species composition at the field farm parcels

#### Variables:			
|Short name	| type |	Long name	| Description |
| ----------|------|------------| ------------|
| species	  |character | species name	    | species name                             |
| Parcel_name | character | farm parcels|  farm parcels |
| cover     | numeric  | 	plant species cover in each parcel | cover (%) of each species at each parcel	 |


### ->  `Community_composition_DungExperiment.csv`
Contains plant species composition germinated at the seed experiment

#### Variables:			
|Short name	| type |	Long name	| Description |
| ----------|------|------------| ------------|
| species	  |character | species name	    | species name                             |
| Parcel_name | character | farm parcels|  farm parcels |
| abundance     | numeric  | 	plant species abundance in each parcel | abundance of each species at each parcel	 |




### ->  `NMDS_data.csv`
Contains NMDS scores from the [NMDS analysis](analysis/2.0_NMDS.R) 

#### Variables:			
|Short name	| type |	Long name	| Description |
| ----------|------|------------| ------------|
| Parcel_name	     | character  | parcel name	|parcel and farm ID |              
| NMDS1_field  | numeric    | NMDS1 score | scores for the 1st NMDS axis from the NMDS analysis of the species composition of plant community at each parcel, measured at the field	|
| NMDS2_field  | numeric    | NMDS2 score | scores for the 2nd NMDS axis from the NMDS analysis of the species composition of plant community at each parcel, measured at the field	|
| NMDS1_exper  | numeric    | NMDS1 score | scores for the 1st NMDS axis from the NMDS analysis of the species composition of plant seedlings  germinated at seed experiment	|
| NMDS2_exper  | numeric    | NMDS2 score | scores for the 2nd NMDS axis from the NMDS analysis of the species composition of plant seedlings  germinated at seed experiment	|



### ->  `Divers_LandUse_Soil_Variables.csv`
Contains environmental variables and biodiversity measures for all plots

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
