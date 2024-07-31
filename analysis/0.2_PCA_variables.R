# Purpose:  PCA for soil data


rm(list = ls())

library(tidyverse)
library(ggplot2)
library(viridis)          

library(factoextra)

# Read data  ----


Dat <- read_csv("data/Panoara_Dat.csv") 
names(Dat)
soil_dat <- Dat %>% 
  dplyr::select(humus,  soil_pH,
                soil_Ca_acces , soil_MG_acces,                                            
               # soil_K_acces, 
               soil_P_acces,  
                soil_CN, soil_N_tot)



# filter(!Parcel_name=="Brade_1") # extreme outlier

str(soil_dat)
names(soil_dat)



soil_pca <- prcomp(soil_dat, #  %>% dplyr::select(soil_P_acces, soil_CN),
                   scale = TRUE)

soil_pca
# round(soil_pca$rotation[, 1], 2)

summary(soil_pca)

# Plots
fviz_eig(soil_pca)

fviz_pca_var(soil_pca,
             axes = c(2, 1),
             # col.var = "cor", # Color by contributions to the PC
             # gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE # Avoid text overlapping
)

fviz_pca_var(soil_pca,
             axes = c(1, 2),
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "yellow", "red", "blue"),
             repel = TRUE # Avoid text overlapping
)

soil_text_var <- get_pca_var(soil_pca)
soil_text_var$contrib
soil_text_var$cor
soil_text_ind <- get_pca_ind(soil_pca)
soil_text_ind
soil_text_ind$coord

round(soil_text_var$cor, 2)
head(soil_text_var$contrib)
head(soil_text_var$cor)


library("corrplot")
corrplot(soil_text_var$cor, #[,c("Dim.1", "Dim.2")], 
         is.corr=F,  
         cl.pos = 'r', cl.ratio = 0.8, tl.col="black") 
corrplot(soil_text_var$contrib[,c("Dim.1", "Dim.2")], is.corr=F,
         cl.ratio = 0.8, tl.col="black")

col <- colorRampPalette(c("#00AFBB", "#E7B800", "#FC4E07"))

corrplot(soil_text_var$contrib, #[,c("Dim.1", "Dim.2")], 
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

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))

corrplot(soil_text_var$cor, #[,c("Dim.1", "Dim.2")], 
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



soil_text_var$cor




# Extract the PC axis coordinates and bind them to the dataset
pca1_soil_text <- soil_text_ind$coord[, 1]
pca2_soil_text <- soil_text_ind$coord[, 2]

soil_text <- bind_cols(Dat %>% dplyr::select(Parcel_name),
                       PC1_soil = pca1_soil_text,
                       PC2_soil = pca2_soil_text)


write_csv(soil_text, "data/soil_PC.csv")


