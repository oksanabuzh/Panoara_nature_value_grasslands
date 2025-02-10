# Purpose: run GLMMs for Dung Experiment Data

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
library("viridis")          


# set theme for plots
set_theme(base = theme_bw(),
          axis.textsize.x = 1, axis.textsize.y = 1, axis.textcolor = "black",
          axis.title.color = "black", axis.title.size = 1.4,
          geom.linetype = 1) #legend.pos = "None", 


# Read data  ------------------------------------------------------------------

Dat <- read_csv("data/Divers_LandUse_Soil_Variables.csv") %>%
  mutate(Grazer_type=str_replace_all(Grazer_type, "_", ", ")) %>% 
  mutate(Mowing_delay=fct_relevel(Mowing_delay,c("no mowing","June","July-August"))) %>%
  arrange(Mowing_delay) %>% 
  mutate(Dung_cover_log=log(Dung_cover+1)) %>% 
  mutate(Plant_SR_exper=case_when(is.na(Plant_SR_exper) ~ 0, 
                                  .default=Plant_SR_exper),
         Plant_abundance_exper=case_when(is.na(Plant_abundance_exper) ~ 0, 
                                         .default=Plant_abundance_exper)) %>% 
   filter(!Parcel_name=="Farm_F_1") # extreme outlier

str(Dat)
names(Dat)
Dat %>% 
  pull(Plant_abundance_exper)


# GLMMs Species Richness (Dung Experiment Data) -------------------------------------

m1_SR_exp <- glmer (Plant_SR_exper ~   Plant_abundance_exper + (1|Farm), 
                    family = "poisson", 
                    data = Dat) 
ranef(m1_SR_exp) # random effects are zeros
summary(m1_SR_exp)
# Remove random effects

m2_SR_exp <- glm (Plant_SR_exper ~   Plant_abundance_exper, family = "poisson", 
                    data = Dat)
anova(m1_SR_exp, m2_SR_exp)


Anova(m2_SR_exp)
summary(m2_SR_exp)

plot_model(m2_SR_exp,type = "pred", terms = c("Plant_abundance_exper"), show.data=T, dot.alpha=0.3, title="")

SR_exp_pred <- get_model_data(m2_SR_exp,type = "pred", terms="Plant_abundance_exper[0:50, by=0.1]")

ggplot(SR_exp_pred, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1)+
  geom_point(data=Dat,
             aes(Plant_abundance_exper, Plant_SR_exper), fill = "#64ABCE",
             size=4, alpha=0.7, pch=21)+
  labs(y="Seedling species number", x='Seedling abundance')+
  geom_line(linetype=1, linewidth=1) 



# GLMMs Abundance (Dung Experiment Data) --------------------------------------------


m1_Abund_exp <- glmer (Plant_abundance_exper ~  log(Grazing_intensity_A+1) + 
                         Grazer_diversity +   
                         (1|Farm), family = "poisson", 
                       data = Dat %>% filter(!Parcel_name=="Farm_F_1")) # strong outlier in Plant_SR_exper and Plant_abundance_exper

  check_convergence(m1_Abund_exp)
ranef(m1_Abund_exp) # random effects are not zeros
hist(ranef(m1_Abund_exp)$`Farm`[,1])

check_collinearity(m1_Abund_exp)


m2_Abund_exp <- glm (Plant_abundance_exper ~  log(Grazing_intensity_A+1) + Grazer_diversity,  # Grazer_type + 
                     family = "poisson", 
                     data = Dat %>% filter(!Parcel_name=="Farm_F_1")) # strong outlier in Plant_SR_exper and Plant_abundance_exper
anova(m1_Abund_exp, m2_Abund_exp)



Anova(m1_Abund_exp)
summary(m1_Abund_exp)

# Grazing_intensity

plot_model(m1_Abund_exp,type = "pred", terms = c("Grazing_intensity_A"), show.data=T, dot.alpha=0.3, title="")
Abund_exp_pred1 <- get_model_data(m1_Abund_exp,type = "pred", terms="Grazing_intensity_A[0:1000, by=0.1]")

ggplot(Abund_exp_pred1, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1)+
  geom_point(data=Dat %>% filter(!Parcel_name=="Farm_F_1")%>% 
               mutate(Grazer_type=str_replace_all(Grazer_type_specifc, "_", ", ")),
             aes(Grazing_intensity_A, Plant_abundance_exper, fill = Grazer_type),
             size=3, alpha=0.7, pch=21)+
  scale_fill_viridis(discrete=TRUE, option = "D")+
  labs(y="Abundance", x='Grazing intensity', fill="Grazer type")+
  geom_line(linetype=1, linewidth=1)


# Grazer_diversity

plot_model(m1_Abund_exp,type = "pred", terms = c("Grazer_diversity"), show.data=T, dot.alpha=0.3, title="")
Abund_exp_pred2 <- get_model_data(m1_Abund_exp,type = "pred", terms="Grazer_diversity[1:4, by=0.1]")

ggplot(Abund_exp_pred2, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1)+
  geom_point(data=Dat %>% filter(!Parcel_name=="Farm_F_1")%>% 
               mutate(Grazer_type=str_replace_all(Grazer_type_specifc, "_", ", ")),
             aes(Grazer_diversity, Plant_abundance_exper, fill = Grazer_type),
             size=3.5, alpha=0.7, pch=21,
             position = position_jitter(width = 0.1,
                                        height=0))+
  scale_fill_viridis(discrete=TRUE, option = "D")+
  labs(y="Seedling abundance", x='Grazer diversity', fill="Combinations of grazer types")+
  geom_line(linetype=1, linewidth=1)

