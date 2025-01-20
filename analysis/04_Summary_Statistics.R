# Purpose: Data wrangling and summary statistics
## Correlations among the measures of plant community, using GLMMs

rm(list = ls()) # clears working environment

# dev.off() # shuts down the current graphical device

# load packages

library(tidyverse)
library(ggplot2)
library(sjPlot)
library(patchwork)
library(performance)
library(lme4)
library(lmerTest)
library(car)
library(emmeans)
library(multcomp)
library(piecewiseSEM)
library(viridis)          


# set thepe for plots
set_theme(base = theme_bw(),
          axis.textsize.x = 1, axis.textsize.y = 1, axis.textcolor = "black",
          axis.title.color = "black", axis.title.size = 1.4,
          geom.linetype = 1) #legend.pos = "None", 




# (1) Summary for plant community ----

## field data ----
Community_VP_field <- read_csv("data/Community_composition_VegetationPlots.csv") %>% 
  rename(layer="layer (VP - vascular, B - bryophyte, L - lichen)") %>% 
  filter(layer == "VP") %>% # only for vascular plants 
  dplyr::select(-layer) %>% 
  pivot_longer(!species, names_to = "Parcel_name", values_to = "cover") %>% 
 # rename(Species="Plot_Parcel") %>% 
  mutate(Parcel_name=str_replace_all(Parcel_name, " ", "_")) %>% 
  filter(!is.na(cover))

Community_VP_field %>% pull(Parcel_name)%>% unique()

str(Community_VP_field)


Community_VP_field %>% 
  group_by(Parcel_name) %>% 
  count() # species number in each parcel


Community_VP_field %>% 
  dplyr::select(species) %>%
  distinct()
  
Sp.occur <- Community_VP_field %>% 
  group_by(species) %>% 
  count() %>% 
  mutate(fraction=round(n*100/44,0)) %>% # species occurrence (% parcels where species was found) %>% 
  arrange(desc(fraction)) %>% 
mutate(abbreviation = str_c(str_split_i(species, '\\s', 1) %>% str_sub(., 1, 4), # make short species names Genus species -> Gen.spe
                          str_split_i(species, '\\s', 2) %>% str_sub(., 1, 3), sep = '.')) %>% 
  relocate(abbreviation, .after = species) %>% 
  rename(parcel_number=n)

  
Sp.occur

write_csv(Sp.occur, "results/Sp.occurances_field.csv")


Sp.occur$Sp_ordered <- with(Sp.occur, reorder(species, fraction))


  ggplot(Sp.occur, aes(y =Sp_ordered , x = fraction)) +
  geom_bar(stat="identity", color="#64ABCE", fill="#64ABCE", alpha=.6, width=.4) +
    theme(axis.text.x = element_text(size=6), axis.text.y = element_text(size=4),
          axis.title=element_text(size=9)) +
 #  coord_flip() +
  xlab("") +
 # theme_bw() +
  labs(y=" ", color="", fill="",
       x="% parcels, where species occures")

## experiment data ----

Community_exper <- read_csv("data/Community_composition_DungExperiment.csv") %>% 
    rename(layer="layer (VP - vascular, B - bryophyte, L - lichen)") %>% 
    filter(layer == "VP") %>% # only for vascular plants 
    dplyr::select(-layer) %>% 
    pivot_longer(!species, names_to = "Parcel_name", values_to = "abundance") %>% 
    mutate(Parcel_name=str_replace_all(Parcel_name, " ", "_")) %>% 
    filter(!is.na(abundance)) 
  
Community_exper %>% pull(Parcel_name)%>% unique()
  

Sp.occur_exp <- Community_exper %>% 
  group_by(species) %>% 
  count() %>% 
  mutate(fraction=round(n*100/28,0)) %>% # species occurrence (% parcels where species was found) %>% 
  arrange(desc(fraction))

Sp.occur_exp

43*100/221 # % species dispersed through zoochory via excrements and farmyard dung

Sp.occur_exp$Sp_ordered <- with(Sp.occur_exp, reorder(species, fraction))


ggplot(Sp.occur_exp, aes(y =Sp_ordered , x = fraction)) +
  geom_bar(stat="identity", color="#64ABCE", fill="#64ABCE", alpha=.6, width=.4) +
  theme(axis.text.x = element_text(size=9), axis.text.y = element_text(size=9),
        axis.title=element_text(size=11)) +
  #  coord_flip() +
  xlab("") +
  # theme_bw() +
  labs(y="Seedling species", color="", fill="",
       x="% parcels, where species occured")

#  (2) Correlations among plant community measures  -----

Dat <- read_csv("data/Panoara_Dat.csv") %>%
  mutate(Grazer_type=str_replace_all(Grazer_type, "_", ", ")) %>% 
  mutate(Grazer_type=fct_relevel(Grazer_type,c("cow","sheep, goat","mixed"))) %>%
  mutate(Mowing_delay=fct_relevel(Mowing_delay,c("no mowing","June","July-August"))) %>%
  arrange(Mowing_delay) %>% 
  mutate(Dung_cover_log=log(Dung_cover+1)) %>% 
  mutate(SR_D_E_exper=case_when(is.na(SR_D_E_exper) ~ 0, .default=SR_D_E_exper),
         Abund_D_E_exper=case_when(is.na(Abund_D_E_exper) ~ 0, .default=Abund_D_E_exper)) %>% 
  mutate(Grazing_int_log = log1p(Grazing_intensity_A),
         Corralling=factor(Corralling))

str(Dat)
names(Dat)


Dat <- Dat %>% 
  mutate(zz = Plant_SR_vascular) # for predicting from the glmmPQL models

#


## Abundance ----

m1_abund <- glmer(Plant_SR_vascular ~  CoverVP_field + 
                    (1|Farm), family = "poisson", data = Dat) 

check_overdispersion(m1_abund) 
check_convergence(m1_abund) 

Anova(m1_abund)
summary(m1_abund)

# overdispersed data, use quasi

m2_abund <- glmmPQL(Plant_SR_vascular ~  CoverVP_field, random = ~ 1 | Farm,  data = Dat,
                    family = "quasipoisson") 

Anova(m2_abund)

model_data(m2_abund)

plot_model(m2_abund, type = "pred", terms = c("CoverVP_field"), show.data=T, dot.alpha=0.3, title="") 

m2_abund_pred <- get_model_data(m2_abund,type = "pred", terms="CoverVP_field[50:190, by=0.1]")

ggplot(m2_abund_pred, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1)+
  geom_point(data=Dat, aes(CoverVP_field, Plant_SR_vascular, fill = Mowing_delay), 
             fill ="#64ABCE", size=3, alpha=0.5, pch=21,
             position=position_jitter(width=0.05, height=0))+
  # scale_fill_manual(values=c("blue", "red", "green"))+
  theme(axis.title.x=element_text(vjust=-0.1), axis.title.y=element_text(vjust=2)) +
  labs(y="Species richness", x="Plant cover, %")+
  geom_line(linetype=1, linewidth=1) 


## Evenness ----

m1_Even <- glmer(Plant_SR_vascular ~  EvennessVP_field + 
                   (1|Farm), family = "poisson", data = Dat) 

check_overdispersion(m1_Even) 
check_convergence(m1_Even) 

Anova(m1_Even)
summary(m1_Even)


plot_model(m1_Even, type = "pred", show.data=T, dot.alpha=0.3, title="") 

m1_Even_pred <- get_model_data(m1_Even,type = "pred", terms="EvennessVP_field[2:30, by=0.1]")

ggplot(m1_Even_pred, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1)+
  geom_point(data=Dat, aes(EvennessVP_field, Plant_SR_vascular), 
             fill ="#64ABCE", size=3, alpha=0.5, pch=21,
             position=position_jitter(width=0.05, height=0))+
  # scale_fill_manual(values=c("blue", "red", "green"))+
  theme(axis.title.x=element_text(vjust=-0.1), axis.title.y=element_text(vjust=2)) +
  labs(y="Species richness", x="Evenness")+
  geom_line(linetype=1, linewidth=1) 




## Shannon  ----

m1_Shan <- glmer(Plant_SR_vascular ~  ShannonVP_field + 
                   (1|Farm), family = "poisson", data = Dat) 

check_overdispersion(m1_Shan) 
check_convergence(m1_Shan) 

Anova(m1_Shan)
summary(m1_Shan)
max(Dat$ShannonVP_field)

plot_model(m1_Shan, type = "pred", show.data=T, dot.alpha=0.3, title="") 

m1_Shan_pred <- get_model_data(m1_Shan,type = "pred", terms="ShannonVP_field[2:3.72, by=0.1]")

ggplot(m1_Shan_pred, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1)+
  geom_point(data=Dat, aes(ShannonVP_field, Plant_SR_vascular), 
             fill ="#64ABCE", size=3, alpha=0.5, pch=21,
             position=position_jitter(width=0.05, height=0))+
  # scale_fill_manual(values=c("blue", "red", "green"))+
  theme(axis.title.x=element_text(vjust=-0.1), axis.title.y=element_text(vjust=2)) +
  labs(y="Species richness", x="Shannon diversity")+
  geom_line(linetype=1, linewidth=1) 


## Biomass ----

m1_biom <- glmer(Plant_SR_vascular ~   poly(biom_log, 2) + 
                   (1|Farm), family = "poisson", data = Dat%>% 
                   mutate(biom_log=log1p(Biomass_dry_weight))) 

check_overdispersion(m1_biom) 
check_convergence(m1_biom) 

summary(m1_biom)
Anova(m1_biom)

# overdispersed data, use quasi/nb

m2_biom <- glmer.nb(Plant_SR_vascular ~  poly(biom_log, 2) +
                      (1|Farm),
                    data = Dat%>% 
                      mutate(biom_log=log1p(Biomass_dry_weight))) 


# m2_biom <- glmmPQL(Plant_SR_vascular ~ biom_log ,# poly(biom_log, 2) , 
#                   random = ~ 1 | Farm,  data = Dat %>% mutate(biom_log=log1p(Biomass_dry_weight)),
#                    family = "quasipoisson") 

car::Anova(m2_biom)


plot_model(m2_biom, type = "pred", show.data=T, dot.alpha=0.3, title="") 

m2_biom_pred <- get_model_data(m2_biom,type = "pred", terms="biom_log")


ggplot(m2_biom_pred, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1)+
  geom_point(data=Dat%>% 
               mutate(biom_log=log1p(Biomass_dry_weight)), aes(biom_log, Plant_SR_vascular), 
             fill ="#64ABCE", size=3, alpha=0.5, pch=21,
             position=position_jitter(width=0.05, height=0))+
  # scale_fill_manual(values=c("blue", "red", "green"))+
  theme(axis.title.x=element_text(vjust=-0.1), axis.title.y=element_text(vjust=2)) +
  labs(y="Species richness", x="Plant biomass (log)")+
  geom_line(linetype=1, linewidth=1) 





#####


### Seedlings

Data <- read_csv("data/Panoara_Dat.csv") %>%
  dplyr::select(Parcel_name, SR_D_E_exper, Abund_D_E_exper, D_E_exp_categ, D_E_exp_categ_2) %>% 
  filter(!D_E_exp_categ_2=="n") %>%  # remove parcels with no seedlings
  filter(!D_E_exp_categ_2=="c+s") %>%  # remove parcels with combinations of sheep and cow
  mutate(experiment=
           case_when(D_E_exp_categ_2=="c" & D_E_exp_categ=="dung" ~ "cow manure",
                     D_E_exp_categ_2=="c" & D_E_exp_categ=="excr" ~ "cow feces",
                     D_E_exp_categ_2=="s" & D_E_exp_categ=="excr" ~ "sheep feces")) %>% 
  filter(!Abund_D_E_exper>150) # remove outlire
  
 

ggplot(Data, aes(experiment, Abund_D_E_exper)) + 
  geom_boxplot(outlier.shape = NA, notch = F)+
  geom_point(alpha=0.7, pch=21, , fill ="#64ABCE",# size=4, #fill="gray", 
             position=position_jitter(width = 0.1, height = 0)) +
  # theme_bw()+
  labs(x =" ", y="Abundance of seedlings")+
  theme(axis.text = element_text(size=10), axis.title=element_text(size=12)) +
  theme(axis.title.x=element_text(vjust=-0.7), 
        axis.title.y=element_text(vjust=1.5),
        legend.key=element_blank()) 



ggplot(Data, aes(experiment, SR_D_E_exper)) + 
  geom_boxplot(outlier.shape = NA, notch = F)+
  geom_point(alpha=0.7, pch=21, , fill ="#64ABCE",# size=4, #fill="gray", 
             position=position_jitter(width = 0.1, height = 0)) +
  # theme_bw()+
  labs(x =" ", y="Species richness of seedlings")+
  theme(axis.text = element_text(size=10), axis.title=element_text(size=12)) +
  theme(axis.title.x=element_text(vjust=-0.7), 
        axis.title.y=element_text(vjust=1.5),
        legend.key=element_blank()) 
