# Purpose:  PCA analysis for soil nutrient data


rm(list = ls())

# load packages ----

library(tidyverse)
library(ggplot2)
library(viridis)          
library(factoextra)

# data  ----

Dat <- read_csv("data/Biomass_nutrients.csv") %>% 
  mutate(Parcel_name=str_replace_all(Parcel_name, " ", "_")) %>% 
  mutate(plant_CN=C_Wassen/N_Wassen) %>% 
  mutate(plant_CN=case_when(plant_CN=="NaN"~0,
                            .default = plant_CN))
names(Dat)




plant_Dat <- Dat %>% 
  dplyr::select(P_Wassen, K_Wassen, plant_CN, Biomass_Wassen) %>% 
  rename("Plant P" = P_Wassen,
         "Plant C:N" = plant_CN,
         "Plant biomass" = Biomass_Wassen,
         "Plant K" = K_Wassen)

cor(plant_Dat)

  
# PCA analysis ----

plant_pca <- prcomp(Dat %>% 
                      dplyr::select(P_Wassen, K_Wassen) %>% 
                      rename("Plant P" = P_Wassen,
                            # "Plant C:N" = plant_CN,
                            # "Plant biomass" = Biomass_Wassen,
                             "Plant K" = K_Wassen),
                   scale = TRUE)

plant_pca

summary(plant_pca)

# Plots----

fviz_eig(plant_pca)

fviz_pca_var(plant_pca,
             axes = c(1, 2),
        #     col.var = "cos2", # Color by contributions to the PC1
         #   gradient.cols = c("#00AFBB", "yellow", "red", "blue"),
             repel = TRUE # Avoid text overlapping
)

plant_var <- get_pca_var(plant_pca)
plant_var$contrib
plant_var$cor
plant_ind <- get_pca_ind(plant_pca)
plant_ind
plant_ind$coord

round(plant_var$cor, 2)
head(plant_var$contrib)
head(plant_var$cor)


library("corrplot")

## Correlation with the PC axes ----

corrplot(plant_var$cor, #[,c("Dim.1", "Dim.2")], 
         is.corr=F,  
         cl.pos = 'r', cl.ratio = 0.8, tl.col="black") 

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))

corrplot(plant_var$cor,
         is.corr=F,  method="color", 
         col=col(200),   cl.ratio = 1,  
         #type="lower", # order="hclust", 
         addCoef.col = "black", number.cex = 0.7, # Add coefficient of correlation
         tl.col="black", tl.srt=45 #, #Text label color and rotation
         # Combine with significance
         # p.mat = p.mat, sig.level = 0.01, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         #  diag=FALSE 
)



## Contribution to the PC axes ----
corrplot(plant_var$contrib[,c("Dim.1", "Dim.2", "Dim.3")], is.corr=F,
         cl.ratio = 0.8, tl.col="black")

col <- colorRampPalette(c("#00AFBB", "#E7B800", "#FC4E07"))
corrplot(plant_var$contrib, #[,c("Dim.1", "Dim.2")], 
         is.corr=F,  method="color", 
         col=col(6),   cl.ratio = 1,
         # type="lower", # order="hclust", 
         addCoef.col = "black", number.cex = 0.7, # Add coefficient of correlation
         tl.col="black", tl.srt=45 #, #Text label color and rotation
         # Combine with significance
         # p.mat = p.mat, sig.level = 0.01, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         #  diag=FALSE 
)



# Extract the PC axis coordinates and bind them to the dataset ----

pca1_plant <- plant_ind$coord[, 1] # use only the 1st axis

plant_data <- bind_cols(Dat %>% dplyr::select(Parcel_name, Biomass_Wassen, plant_CN),
                       plant_PK_PC = pca1_plant) %>% 
  rename(plant_biomass= Biomass_Wassen)


write_csv(plant_data, "data/plant_biomass_CN_PK.PCA.csv")

