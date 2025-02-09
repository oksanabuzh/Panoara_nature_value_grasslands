# Purpose:  PCA analysis for soil nutrient data

# rm(list = ls())
# load packages ----

library(tidyverse)
library(ggplot2)
library(viridis)          
library(factoextra)

# data  ----

Dat <- read_csv("data/Divers_LandUse_Soil_Variables.csv") 

names(Dat)

soil_Dat <- Dat %>% 
  dplyr::select(soil_P_acces, soil_CN,  soil_K_acces) %>% 
  rename("Soil P" = soil_P_acces,
         "Soil C:N" = soil_CN,
         "Soil K" = soil_K_acces)
  
# PCA analysis ----

soil_pca <- prcomp(soil_Dat,
                   scale = TRUE)

soil_pca

summary(soil_pca)

# Plots----

fviz_eig(soil_pca)

fviz_pca_var(soil_pca,
             axes = c(1, 2),
        #     col.var = "cos2", # Color by contributions to the PC1
         #   gradient.cols = c("#00AFBB", "yellow", "red", "blue"),
             repel = TRUE # Avoid text overlapping
)

soil_var <- get_pca_var(soil_pca)
soil_var$contrib
soil_var$cor
soil_ind <- get_pca_ind(soil_pca)
soil_ind
soil_ind$coord

round(soil_var$cor, 2)
head(soil_var$contrib)
head(soil_var$cor)


library("corrplot")

## Correlation with the PC axes ----

corrplot(soil_var$cor, #[,c("Dim.1", "Dim.2")], 
         is.corr=F,  
         cl.pos = 'r', cl.ratio = 0.8, tl.col="black") 

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))

corrplot(soil_var$cor,
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
corrplot(soil_var$contrib[,c("Dim.1", "Dim.2", "Dim.3")], is.corr=F,
         cl.ratio = 0.8, tl.col="black")

col <- colorRampPalette(c("#00AFBB", "#E7B800", "#FC4E07"))
corrplot(soil_var$contrib, #[,c("Dim.1", "Dim.2")], 
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

pca1_soil <- soil_ind$coord[, 1] # use only the 1st axis

soil_text <- bind_cols(Dat %>% dplyr::select(Parcel_name),
                       soil_NPK_PC = pca1_soil)

write_csv(soil_text, "data/soil_NPK_PCA.csv")

