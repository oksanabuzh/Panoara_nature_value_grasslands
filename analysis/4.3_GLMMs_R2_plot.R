# Purpose: Create R2 plot for the GLMMs (Fig. 3)

library(tidyverse)

# data 
R2_partial <- read_csv("results/Tables_R2_partial.csv")%>% 
  mutate(predictor=fct_relevel(predictor,
                                  c("Management type",
                                    "Management stability",
                                    "Mowing frequency",
                                    "Mowing delay",
                                    "Grazing intensity",
                                    "Grazer type", 
                                    "Grazing type",
                                    "Grazing season", 
                                    "Grazing legacy",
                                    "Corralling",
                                    "Ploughing", 
                                    "Manuring frequency", 
                                    "Cow manure application",
                                    "Litter removal",
                                    "Ant/mole-hill leveling",
                                    "Shrub/tree removal",
                                    "Moss removal",
                                    "Burning", 
                                    "Soil NPK" ,
                                    "Humus"))) %>%  
  mutate(variab_group=recode_factor(predictor,
                                    "Management type" ="Management type/stability/legacy",
                                    "Management stability" = "Management type/stability/legacy",
                                    "Mowing frequency" = "Mowing",
                                    "Mowing delay" = "Mowing",
                                    "Grazing intensity" = "Grazing",
                                    "Grazer type" = "Grazing",
                                    "Grazing type" = "Grazing",
                                    "Grazing season" = "Grazing",
                                    "Grazing legacy" = "Management type/stability/legacy",  
                                    "Manuring frequency" = "Manuring",
                                    "Cow manure application" = "Manuring",
                                    "Litter removal" = "Cleaning",
                                    "Ant/mole-hill leveling" = "Cleaning",
                                    "Shrub/tree removal" = "Cleaning",
                                    "Moss removal" = "Cleaning",
                                    Burning = "Cleaning",
                                    Corralling = "Management type/stability/legacy",
                                    Ploughing = "Management type/stability/legacy",
                                    "Soil NPK" = "Soil properties",
                                    "Humus" = "Soil properties")) %>%  
  arrange(predictor) 


R2_partial





## > Plot R2 ----

# set theme for plots
set_theme(base = theme_bw(),
          axis.textsize.x = 0.9, axis.textsize.y = 0.9, axis.textcolor = "black",
          axis.title.color = "black", axis.title.size = 1.2,
          geom.linetype = 1) #legend.pos = "None", 


dodge_width <- 0.5

plot <- ggplot(R2_partial,
               aes(y =reorder(predictor, desc(predictor)), x = Rsq, 
                   color = variab_group)) +
  geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
  geom_point(position = position_dodge(width = dodge_width), size = 6) +
  geom_errorbarh(aes(xmin = 0, xmax = Rsq),
                 position = position_dodata:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABIAAAASCAYAAABWzo5XAAAAWElEQVR42mNgGPTAxsZmJsVqQApgmGw1yApwKcQiT7phRBuCzzCSDSHGMKINIeDNmWQlA2IigKJwIssQkHdINgxfmBBtGDEBS3KCxBc7pMQgMYE5c/AXPwAwSX4lV3pTWwAAAABJRU5ErkJggg==dge(width = dodge_width), height = 0.1) +
  theme(legend.key=element_blank()) +
  # theme_bw()+
  theme(axis.text.x = element_text(size=10), 
        axis.text.y = element_text(size=11),
        axis.title=element_text(size=13, face="bold"), 
        legend.text = element_text(size=10)) +
  #  theme(axis.title.x=element_text(vjust=-0.1), axis.title.y=element_text(vjust=2))+
  labs(y="Driver", color="Driver type",
       x= expression(paste("Variance explained, ", R^{2})),
       title = "Relative importance of biodiversity drivers") +
  scale_color_manual(values=c("cadetblue4",
                               "orange","#F8766D","#9BB018",
                              # "dodgerblue3","dodgerblue3",
                              "slateblue3",
                               "#BA7938"))

  

plot


# End ---------------


