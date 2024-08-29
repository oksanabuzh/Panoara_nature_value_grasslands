# Purpose: Calculate evenness and community composition (NMDS) for each plot


rm(list = ls()) # clears working environment

# dev.off() # shuts down the current graphical device

# load packages

library(tidyverse)
library(ggplot2)
library(vegan)

library(ggplot2)
library(ggrepel)
library(devtools)

# (1) Calculate diversity measures ----

## 1.1 Field data----

# read data
Community_VP_field <- read_csv("data/Community_composition_VegetationPlots.csv") %>% 
  rename(layer="layer (VP - vascular, B - bryophyte, L - lichen)") %>% 
  filter(layer == "VP") %>% # only for vascular plants 
  dplyr::select(-layer) %>% 
  pivot_longer(!Plot_Parcel, names_to = "Parcel_name", values_to = "cover") %>% 
  rename(Species="Plot_Parcel") %>% 
  mutate(Parcel_name=str_replace_all(Parcel_name, " ", "_")) %>% 
  filter(!is.na(cover))
  
Community_VP_field %>% pull(Parcel_name)%>% unique()

str(Community_VP_field)

# calculate diversity measures
VP_field <- Community_VP_field %>% 
  group_by(Parcel_name) %>% 
  summarise(
    CoverVP_field=sum(cover),
    SR_VP_field = n_distinct(Species),
    EvennessVP_field = vegan::diversity(cover, index = "invsimpson"),
    ShannonVP_field = vegan::diversity(cover, index = "shannon")) %>%
  ungroup()

VP_field


## 1.2 Experiment data----
# read data
Community_exper <- read_csv("data/Community_composition_DungExperiment.csv") %>% 
  rename(layer="layer (VP - vascular, B - bryophyte, L - lichen)") %>% 
  filter(layer == "VP") %>% # only for vascular plants 
  dplyr::select(-layer) %>% 
  pivot_longer(!Plot_Parcel, names_to = "Parcel_name", values_to = "abundance") %>% 
  rename(Species="Plot_Parcel") %>% 
  mutate(Parcel_name=str_replace_all(Parcel_name, " ", "_")) %>% 
  filter(!is.na(abundance)) 

Community_exper %>% pull(Parcel_name)%>% unique()

str(Community_exper)

# calculate diversity measures
VP_Exper <- Community_exper %>% 
  group_by(Parcel_name) %>% 
  summarise(
    abundance_Exper=sum(abundance),
    SR_Exper = n_distinct(Species),
    Evenness_Exper = vegan::diversity(abundance, index = "invsimpson"),
    Shannon_Exper = vegan::diversity(abundance, index = "shannon")) %>%
  ungroup()

VP_Exper


Diver_data <-  VP_field %>% 
  left_join(VP_Exper, by="Parcel_name")

#(2) NMDS analysis ----


## 2.1 Field data----
# convert plant list to matrix
compos_VP_field <- Community_VP_field %>% 
  pivot_wider(names_from = "Species", values_from = "cover") %>%
  replace(is.na(.), 0) 

# check if NA
which(is.na(compos_VP_field %>% 
              dplyr::select(-Parcel_name)))

#Check if any species with 0 cover 
which(colSums(compos_VP_field%>% 
                dplyr::select(-Parcel_name),
              na.rm=TRUE) %in% 0)



# assess how well dissimilarity metrics separate (assess) the data
# for this, read the environmental variables and join with the community data
Variables <- read_csv("data/Variables_clean.csv") %>% 
  dplyr::select(Parcel_name, Grazing_intensity_A, Mowing_frequency,
                Manuring_freq, humus)

Variables_dat <- compos_VP_field  %>% 
  left_join(Variables, by="Parcel_name")

vegan::rankindex(compos_VP_field %>% dplyr::select(-Parcel_name),
                 Variables_dat%>% dplyr::select(Grazing_intensity_A, Mowing_frequency,
                                                 Manuring_freq, humus))

### NMDS analysis

set.seed(10)

nmds1<-vegan::metaMDS(compos_VP_field%>% 
                        dplyr::select(-Parcel_name), 
                      distance="bray",k=2,trymax=100)

# Model fit (stress value)
nmds1 # the stress value shows how easy it was to condense multidimensional data into two dimensional space, below 0.2 is generally good

# Shepard plot and model fit (R2)
vegan::stressplot(nmds1, main = "Shepard plot")
nmds1$stress

# fast check plot
plot(nmds1)

vegan::scores(nmds1)

as.data.frame(vegan::scores(nmds1)$sites)
as.data.frame(vegan::scores(nmds1)$species)

#if any NA?
which(is.na(as.data.frame(vegan::scores(nmds1)$species)))

# add the NMDS scores to the data
NMDS_VP_field <- compos_VP_field %>% 
  dplyr::select(Parcel_name) %>% 
  mutate(NMDS1_VP_field = as.data.frame(vegan::scores(nmds1)$sites)$NMDS1,
         NMDS2_VP_field = as.data.frame(vegan::scores(nmds1)$sites)$NMDS2)

# add diversity measures to the data
Diver_NMDS_data <- Diver_data %>% 
  left_join(NMDS_VP_field, by="Parcel_name")

## 2.2 Experiment data----
# read data
compos_exper <- Community_exper %>% 
  pivot_wider(names_from = "Species", values_from = "abundance") %>%
  replace(is.na(.), 0) 

# check if NA
which(is.na(compos_exper %>% 
              dplyr::select(-Parcel_name)))

#Check if any species with 0 cover 
which(colSums(compos_exper%>% 
                dplyr::select(-Parcel_name),
              na.rm=TRUE) %in% 0)

# NMDS analysis
set.seed(11)
nmds2<-vegan::metaMDS(vegan::wisconsin(compos_exper%>% 
                        dplyr::select(-Parcel_name)), 
                      distance="bray",k=2,trymax=100)
# model fit (sterss)
nmds2 # the stress value shows how easy it was to condense multidimensional data into two dimensional space, below 0.2 is generally good
# Shepard plot and model fit (R2)
vegan::stressplot(nmds2, main = "Shepard plot")
nmds2$stress
eigenvals(nmds2)
summary(nmds2)
# fast plot
plot(nmds2)


vegan::scores(nmds2)

as.data.frame(vegan::scores(nmds2)$sites)
as.data.frame(vegan::scores(nmds2)$species)

# any NAs?
which(is.na(as.data.frame(vegan::scores(nmds2)$species)))

# add the NMDS scores to the data
NMDS_exper <- compos_exper %>% 
  dplyr::select(Parcel_name) %>% 
  mutate(NMDS1_exper = as.data.frame(vegan::scores(nmds2)$sites)$NMDS1,
         NMDS2_exper = as.data.frame(vegan::scores(nmds2)$sites)$NMDS2)

# add the biodiversity measures to the data
Diver_NMDS_data <- Diver_data %>% 
  left_join(NMDS_VP_field, by="Parcel_name") %>% 
  left_join(NMDS_exper, by="Parcel_name")


write_csv(Diver_NMDS_data, "data/Diversity_&_NMDS_data.csv")



# (3) PERMANOVA  analysis -----

# read data
Soil_PC <- read_csv("data/soil_PC.csv") %>% 
  mutate(PC1_soil_log= log(PC1_soil +10),
         PC2_soil_log= log(PC2_soil +10))

Soil_PC2 <- read_csv("data/soil_PC2.csv")


Data <- read_csv("data/Panoara_Dat.csv") %>%
  mutate(Grazer_type=str_replace_all(Grazer_type, "_", ", ")) %>% 
  mutate(Grazer_type=fct_relevel(Grazer_type,c("cow","sheep, goat","mixed"))) %>%
  mutate(Mowing_delay=fct_relevel(Mowing_delay,c("no mowing","June","July-August"))) %>%
  arrange(Mowing_delay) %>% 
  mutate(Dung_cover_log=log(Dung_cover+1)) %>% 
  mutate(Grazing_int_log = log1p(Grazing_intensity_A),
         Grazing_legacy =log1p(Grazing_intensity_B),
        # Corralling=case_when(Corralling=="no" ~ "0", 
        #                      Corralling=="yes" ~ "1"),
         Moss_removal=case_when(Moss_removal=="no" ~ 0, 
                                Moss_removal=="yes" ~ 1)) %>% 
  mutate(humus_log=log1p(humus)) %>% 
  mutate(Cow_dung=case_when(Cow_dung_applied=="present" ~ "1", 
                            Cow_dung_applied=="absent" ~ "0")) %>% 
  mutate(Ploughing=1/Last_ploughing)


# filter(!Parcel_name=="Brade_1") # extreme outlier


# add soil PCA to the data

Dat_field <- Data %>% 
  left_join(Soil_PC, by="Parcel_name") %>% 
  left_join(Soil_PC2, by="Parcel_name") %>% 
  mutate(SR_Exper=case_when(is.na(SR_Exper) ~ 0, .default=SR_Exper),
         Abund_D_E_exper=case_when(is.na(Abund_D_E_exper) ~ 0, .default=Abund_D_E_exper),
         abundance_Exper=case_when(is.na(abundance_Exper) ~ 0, .default=abundance_Exper),
         Evenness_Exper=case_when(is.na(Evenness_Exper) ~ 0, .default=Evenness_Exper),
         Shannon_Exper=case_when(is.na(Shannon_Exper) ~ 0, .default=Shannon_Exper)) %>% 
  mutate(Grazing_int_log = log1p(Grazing_intensity_A)) %>% 
  mutate(Mowing_delay=case_when(Mowing_delay=="no mowing" ~ 0,
                                Mowing_delay=="July-August" ~ 1,
                                Mowing_delay=="June" ~ 2)) %>% 
  mutate(Mowing_frq_scal=scale(Mowing_frequency, center=T, scale=T)) %>% 
  mutate(Mowing_frq_sqrd=Mowing_frq_scal^2,
         SR_Exper_log = log1p(SR_Exper+1))




Dat_exp <- Dat_field %>% filter(!Parcel_name=="Brade_1") %>%    # extrime outlyer
  filter(!is.na(SR_D_E_exper)) %>% 
  filter(!abundance_Exper==0)




## 3.1. PERMANOVA   field ----
summary(Dat_field)
names(Dat_field)

set.seed(10)


PERM_mod1 <- vegan::adonis2(compos_VP_field %>% dplyr::select(-Parcel_name) ~ 
                         Grazing_int_log +
                         poly(Mowing_frequency, 2)+
                         Manuring_freq + 
                         humus_log +  
                         PC1_soil_2 +
                     #  Grazer_type +
                       Grazing_season +
           #   Mowing_delay  + 
                       Ploughing + # Crops_planted + # correlates with Last_ploughing
                      # Shrub_tree_removal + # Cleaning +  # Shrub_tree_removal + # correlates with Cleaning
                       Litter_removal +
                      Moss_removal + 
                       Anthill_leveling + #  Molehill_leveling +# data are the same as Anthill_leveling
                     #  Burning + 
                       Corralling , #+
                     #  Cow_dung_applied , 
                      data=Dat_field,  strata=Dat_field$Farm,
                     permutations = 1000, method = "bray")

PERM_mod1

write.csv(PERM_mod1, "results/PERMANOVA_Field_Data.csv")


## 3.2. PERMANOVA   experiment ----

summary(Dat_exp)
names(Dat_exp)

set.seed(1)

PERM_mod2 <- vegan::adonis2(compos_exper %>% filter(!Parcel_name=="Brade_1") %>% dplyr::select(-Parcel_name) ~ 
                              Grazing_int_log +
                              Manuring_freq  + 
                              Grazing_season, 
                            data=Dat_exp,  strata=Dat_exp$Farm,
                            permutations = 1000, method = "bray")

PERM_mod2

write.csv(PERM_mod2, "results/PERMANOVA_Experiment.csv")


# (4) Plot NMDS ---- 


## 4.1. Seed Experiment ----
Dat_exp2 <- Dat_field %>% 
  filter(!is.na(SR_D_E_exper)) %>% 
  filter(!abundance_Exper==0)

# get coordinates
fit2 <- vegan::envfit(nmds2   ~  
                        Grazing_int_log +
                        Manuring_freq  + 
                        Grazing_season , 
                      data=Dat_exp2, perm=1000) #


fit2

# extract species scores
species.scores_exp <-  scores(nmds2, display = 'sp') %>% 
  as_tibble(rownames ='species') %>% 
  mutate(short_name = str_c(str_split_i(species, '\\s', 1) %>% str_sub(., 1, 4), # make short species names Genus species -> Gen.spe
                            str_split_i(species, '\\s', 2) %>% str_sub(., 1, 3), sep = '.'))

species.scores_exp


print(species.scores_exp, n=43)

sites.scores_exp <-  scores(nmds2, display = 'sites') %>% 
  as_tibble(rownames ='sites.scores') %>% 
  mutate(Grazer_type = factor(Dat_exp2$Grazer_type),
         Habitat = factor(Dat_exp2$habitat_corrected),
         Grazing_season = factor(Dat_exp2$Grazing_season))

sites.scores_exp

# calculate centroid for  Grazing_season
centr_exp <- sites.scores_exp %>% 
  group_by(Grazing_season) %>% 
  summarise(NMDS1c=mean(NMDS1),
            NMDS2c=mean(NMDS2)) %>% 
  ungroup()

# merge with site scores
sites.scores_exp <- sites.scores_exp %>% 
  left_join(centr_exp, by="Grazing_season")


p1  <- ggplot(data = sites.scores_exp, 
              aes(x = NMDS1, y = NMDS2)) +
  geom_vline(xintercept = 0, col="grey", linetype="dashed")+
  geom_hline(yintercept = 0, col="grey", linetype="dashed")+
  scale_shape_manual(values=c( 19, 1))+ 
  geom_point (data = species.scores_exp, size = 2, pch=19, colour= "grey") +
  geom_text_repel (data=species.scores_exp,aes(x=NMDS1,y=NMDS2,label=short_name), 
                   size = 3.6, colour= "black") + 
  theme(axis.title = element_text(size = 13, face = "bold", colour = "black"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, 
                                                                        colour = "grey30", linewidth = 0.5), 
        axis.ticks = element_blank(), 
        axis.text = element_text(size = 13, colour = "black"))

p1



# add ellipses showing system type 

p2  <- ggplot(data = sites.scores_exp, 
              aes(x = NMDS1, y = NMDS2)) +
  stat_ellipse(aes(fill=Habitat), alpha=.1,type='t',linewidth =1, geom="polygon")+
  geom_vline(xintercept = 0, col="grey", linetype="dashed")+
  geom_hline(yintercept = 0, col="grey", linetype="dashed")+
   geom_text_repel (data=species.scores_exp, aes(x=NMDS1,y=NMDS2,label=short_name), 
                   size = 3.6, 
                   colour= "forestgreen") + 
  theme(axis.title = element_text(size = 13, face = "bold", colour = "black"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, 
                                                                        colour = "grey30", size=0.5), 
        axis.ticks = element_blank(), 
        axis.text = element_text(size = 13, colour = "black"),
        legend.key = element_blank(), 
        legend.title = element_text(size = 11, face = "bold", colour = "black"), 
        legend.text = element_text(size = 11, colour = "black")) 


p2

# add spiders for grazing season

p3 <- ggplot(data = sites.scores_exp, 
             aes(x = NMDS1, y = NMDS2)) +
  stat_ellipse(aes(fill=Habitat), alpha=.1,type='t',linewidth =1, geom="polygon")+
  geom_vline(xintercept = 0, col="grey", linetype="dashed")+
  geom_hline(yintercept = 0, col="grey", linetype="dashed")+
  
  # add centroids:
  geom_segment(data = sites.scores_exp,        # spiders
               mapping = aes(xend = NMDS1c, yend = NMDS2c, color=Grazing_season),
               alpha=.5) + 
  geom_point(data = centr_exp %>%                   # centroids
               rename(NMDS1=NMDS1c, NMDS2=NMDS2c), size = 5,
             aes(color=Grazing_season),
             alpha=.5) +                        
  geom_point(data = sites.scores_exp, size = 2,      # sample scores
             aes(color=Grazing_season),
             alpha=.5) +                         
  coord_fixed()+
  # add species names
  scale_colour_manual(values = c("#F8766D", "#C77CFF", "#00BFC4"))+
  geom_text_repel (data=species.scores_exp,aes(x=NMDS1,y=NMDS2,label=species #short_name
                                               ), 
                   size = 3.6, 
                   colour= "forestgreen", fontface = 'bold', max.overlaps = Inf) + 
  theme(axis.title = element_text(size = 13, face = "bold", colour = "black"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, 
                                                                        colour = "grey30", size=0.5), 
        axis.ticks = element_blank(), 
        axis.text = element_text(size = 13, colour = "black"),
        legend.key = element_blank(), 
        legend.title = element_text(size = 11, face = "bold", colour = "black"), 
        legend.text = element_text(size = 11, colour = "black")) +
  labs(colour="Grazing season")




p3


#### Add significant predictors ----
names(fit2$vectors)
fit2$vectors$arrows 

# add standardized scores:
coord_cont_exp <- as.data.frame(scores(fit2, "vectors")) %>% 
  cbind(stand= fit2$vectors$arrows) %>% 
  mutate(Variables = rownames(scores(fit2, "vectors"))) %>% 
  mutate(Variables=recode_factor(Variables, 
                                 Grazing_int_log="Grazing",
                                 Manuring_freq="Manuring"))  

coord_cont_exp
# stand.NMDS1 and stand.NMDS2 are all the same lenght arrows for the posthoc plottig



# rescale  all arrows to fill an ordination plot, where fill =  shows proportion of plot to be filled by the arrows
coord_cont_exp_stnd  <- coord_cont_exp %>% 
  mutate(stand.NMDS1=stand.NMDS1 * ordiArrowMul(fit2, rescale=TRUE, fill = 0.4))%>% 
  mutate(stand.NMDS2=stand.NMDS2 * ordiArrowMul(fit2, rescale=TRUE, fill = 0.4))

# Create scores with the standard length for the posthoc plotting of arrows of the same size
# for this we would use raw vectors from envfit and not scores 
# as arrow lengths (scores) are created as NMDS1 vectors multiplied on sqrt(r2) from fit$vectors 


coord_cont_exp_stnd





p4 <- p3 + geom_segment(aes(x = 0, y = 0, xend = stand.NMDS1, yend = stand.NMDS2), 
                  data = coord_cont_exp_stnd, linewidth =0.7, alpha = 1, 
                  colour = "black", 
                  arrow = arrow(length = unit(0.031, "npc"))) +
  geom_text(data = coord_cont_exp_stnd, 
            aes(x =stand.NMDS1, y = stand.NMDS2, label = Variables), 
            colour = "black", fontface = "bold", size = 5 ,
            vjust=c(-1,1),
            hjust=0.5)

p4

ggsave("results/Fig_5.png", p4, width = 20, height = 20, 
       units = "cm")





## 4.2. Field data ----
