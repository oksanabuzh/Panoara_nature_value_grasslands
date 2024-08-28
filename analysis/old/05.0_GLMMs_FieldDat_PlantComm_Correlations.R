# Purpose: Correlations among the measures of plant community, using GLMMs

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


# Read data  ----
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
  # filter(!Parcel_name=="Brade_1") # extreme outlier

str(Dat)
names(Dat)


Dat <- Dat %>% 
  mutate(zz = Plant_SR_vascular) # for predicting from the glmmPQL models



#  Correlations among plant community measures  -----

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






