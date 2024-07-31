# Purpose: to check and run the GLMMs

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




## 1) Species Richness ----
ggplot(Dat, aes(habitat_corrected, Plant_SR_vascular))+geom_boxplot()
ggplot(Dat, aes(Grazing_int_log, Plant_SR_vascular))+geom_point()
ggplot(Dat, aes(Grazer_type, Plant_SR_vascular))+geom_boxplot()
ggplot(Dat, aes(Grazer_diversity, Plant_SR_vascular))+geom_point()
ggplot(Dat, aes(Mowing_frequency, Plant_SR_vascular))+geom_point()
ggplot(Dat, aes(Mowing_delay, Plant_SR_vascular))+geom_boxplot()
ggplot(Dat, aes(Grazing_type, Plant_SR_vascular))+geom_boxplot()
ggplot(Dat, aes(Grazing_season, Plant_SR_vascular))+geom_boxplot()
ggplot(Dat, aes(Litter_removal, Plant_SR_vascular))+geom_boxplot()
ggplot(Dat, aes(Crops_planted, Plant_SR_vascular))+geom_boxplot()
ggplot(Dat, aes(Cleaning, Plant_SR_vascular))+geom_boxplot()
ggplot(Dat, aes(Shrub_tree_removal, Plant_SR_vascular))+geom_boxplot()

ggplot(Dat, aes(log(Last_ploughing), Plant_SR_vascular))+geom_point()

ggplot(Dat, aes(Corralling, Plant_SR_vascular))+geom_boxplot()



# habitat ----
m_SR_1 <- glmer(Plant_SR_vascular ~  habitat_corrected + 
                       (1|Farm), family = "poisson", data = Dat) 
summary(m_SR_1)
Anova(m_SR_1)

check_overdispersion(m_SR_1) # Overdispersion detected.

# Change family to quasipoisson

# library(MASS)
m_SR_1b <- glmmPQL(Plant_SR_vascular ~  habitat_corrected, random = ~ 1 | Farm,  data = Dat,
                family = "quasipoisson") 

PearsonResiduals <- resid(m_SR_1b, type = "pearson")
n_cases <- nrow(Dat) # extract number of samples
n_predictors <- length(fixef(m_SR_1b)) + 1 # extract number of estimates (plus intercept)
Overdispersion <- sum(PearsonResiduals^2) / (n_cases-n_predictors) # calculate overdispersion
Overdispersion # The overdispersion has decreased and is no longer an issue.

plot(m_SR_1b) #  the plot exhibits a slight funnel shape, but not drastically, and thus indicates heteroscedasticity.
qqnorm(resid(m_SR_1b))
qqline(resid(m_SR_1b))

Anova(m_SR_1b)

MuMIn::r.squaredGLMM(m_SR_1b)


emmeans_m1_habitat <- cld(emmeans(m_SR_1b, list(pairwise ~ habitat_corrected)), 
                         Letters = letters) %>% arrange(habitat_corrected)
emmeans_m1_habitat

ggplot(Dat, aes(habitat_corrected, Plant_SR_vascular)) + 
  geom_boxplot(outlier.shape = NA, notch = F)+
  geom_point(aes(fill = Grazer_type, size=Grazing_intensity_A),
    alpha=0.7, pch=21, # size=4, #fill="gray", 
             position=position_jitter(width = 0.1, height = 0)) +
  scale_fill_viridis(discrete=TRUE, option = "D")+
  theme_bw()+
  labs(x ="Management type", y="Species richess", fill="Grazer type", size="Grazing intensity")+
  theme(axis.text = element_text(size=10), axis.title=element_text(size=12)) +
  theme(axis.title.x=element_text(vjust=-0.7), axis.title.y=element_text(vjust=1.5)) +
  geom_text(data=emmeans_m1_habitat,aes(x=habitat_corrected, y=c(89, 87),
                                       label=emmeans_m1_habitat$.group),vjust=0.5, hjust=0, 
            size=4, col="black" , position=position_dodge(0))


# grazing, mowing----

m_SR_2 <- glmer(Plant_SR_vascular ~  
                Grazing_int_log +
                  #  poly(Grazing_int_log , 2)  +
                 poly(Mowing_frequency, 2)  +
                  (1|Farm), family = "poisson", data = Dat) 

plot(m_SR_2)
qqnorm(resid(m_SR_2))
qqline(resid(m_SR_2))

Anova(m_SR_2)
check_collinearity(m_SR_2)

plot_model(m_SR_2, type = "pred", terms = c("Grazing_int_log"), show.data=T, dot.alpha=0.3, title="") 


# check for hump-shape for grazing:
m_SR_2_ <- glmer(Plant_SR_vascular ~  
                   #  Grazing_int_log +
                   poly(Grazing_int_log , 2)  +
                  poly(Mowing_frequency, 2)  +
                  (1|Farm), family = "poisson", data = Dat) 

anova(m_SR_2, m_SR_2_) # selects against the hump-shape for grazing
Anova(m_SR_2)

check_overdispersion(m_SR_2) # Slight overdispersion detected

plot_model(m_SR_2, type = "pred", terms = c("Grazing_int_log"), show.data=T, dot.alpha=0.3, title="") 
plot_model(m_SR_2, type = "pred", terms = c("Mowing_frequency"), show.data=T, dot.alpha=0.3, title="") 

m_SR_2_mowing_pred <- get_model_data(m_SR_2,type = "pred", terms="Mowing_frequency[0:3, by=0.1]")

ggplot(m_SR_2_mowing_pred, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1)+
  geom_point(data=Dat, aes(Mowing_frequency, Plant_SR_vascular, fill = Mowing_delay), 
             size=3, alpha=0.5, pch=21,
             position=position_jitter(width=0.05, height=0))+
  scale_fill_manual(values=c("blue", "red", "green"))+
  theme(axis.title.x=element_text(vjust=-0.1), axis.title.y=element_text(vjust=2)) +
  labs(y="Species richness", x='Mowing frequency', fill="Mowing delay")+
  geom_line(linetype=1, linewidth=1) 


# Grazing intensity

m_SR_2_graz_pred <- get_model_data(m_SR_2,type = "pred", 
                                   terms="Grazing_int_log[2.39:6.9, by=0.001]")

ggplot(m_SR_2_graz_pred, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1)+
  geom_point(data=Dat,
             aes(Grazing_int_log, Plant_SR_vascular, fill = Grazer_type), 
             size=3, alpha=0.7, pch=21)+
  scale_fill_viridis(discrete=TRUE, option = "D")+
  labs(y="Species richness", x='Grazing intensity', fill="Grazer type")+
  geom_line(linetype=1, linewidth=1) 


# R2 for the entire model 
# R2m and R2c are marginal (for fixed predictors) and 
## conditional (for fixed and random predictors) coefficients of determination
MuMIn::r.squaredGLMM(m_SR_2)

# partial R2
R2=r2glmm::r2beta(m_SR_2,  partial = T, method = 'sgv', data=Dat)
R2
plot(R2)

#  grazer_type, mowing_delay ----

m_SR_3 <- glmer(Plant_SR_vascular ~  
                 #  Grazing_int_log +
                 Grazer_type +
                 # Grazing_type +
                  Grazing_season +
                 #   poly(Mowing_frequency, 2)  + # 
                  Mowing_delay + 
                  (1|Farm), family = "poisson", data = Dat) 


plot(m_SR_3)
qqnorm(resid(m_SR_3))
qqline(resid(m_SR_3))

Anova(m_SR_3)
summary(m_SR_3)
check_overdispersion(m_SR_3) #  overdispersion detected
check_collinearity(m_SR_3)



m_SR_3b <- glmer.nb(Plant_SR_vascular ~ 
                     # Grazing_int_log +
                     #
                       Grazer_type +
                      #Grazing_type +
                     # poly(Mowing_frequency, 2)  + # 
                     Mowing_delay + 
                      (1|Farm),
                        data = Dat) 


check_convergence(m_SR_3b)
check_collinearity(m_SR_3b)

Anova(m_SR_3b)
summary(m_SR_3b)
# R2 for the entire model 
# R2m and R2c are marginal (for fixed predictors) and 
## conditional (for fixed and random predictors) coefficients of determination
MuMIn::r.squaredGLMM(m_SR_3b)

# partial R2
r2glmm::r2beta(m_SR_3b,  partial = T, method = 'sgv', data=Dat)

## Grazer_type ----

emmeans_m_SR_3b_Grazer_type <- cld(emmeans(m_SR_3b, list(pairwise ~ Grazer_type)), 
                           Letters = letters, test="Tukey") %>% arrange(Grazer_type)
emmeans_m_SR_3b_Grazer_type

ggplot(Dat, aes(Grazer_type, Plant_SR_vascular)) + 
  geom_boxplot(outlier.shape = NA, notch = F)+
  geom_point(aes(fill = Grazer_type,size=Grazing_intensity_A), alpha=0.7, pch=21,
             position=position_jitter(width = 0.05, height = 0)) +
  scale_fill_viridis(discrete=TRUE, option = "D")+
  theme_bw()+
  labs(x ="Grazer type", y="Species richness", fill="Grazer type", size="Grazing intensity")+
  theme(axis.text = element_text(size=10), axis.title=element_text(size=15)) +
  theme(axis.title.x=element_text(vjust=-0.1), axis.title.y=element_text(vjust=2)) +
  
  geom_text(data=emmeans_m_SR_3b_Grazer_type,aes(x=Grazer_type, y=c(87, 76, 83),
                                         label=emmeans_m_SR_3b_Grazer_type$.group),vjust=0.5, hjust=0, 
            size=4, col="black" , position=position_dodge(0))


### mowing delay ----

emmeans_m_SR_3b_Mowing_delay <- cld(emmeans(m_SR_3b, list(pairwise ~ Mowing_delay)), 
                            Letters = letters, test="Tukey") %>% arrange(Mowing_delay)
emmeans_m_SR_3b_Mowing_delay

ggplot(Dat, aes(Mowing_delay, Plant_SR_vascular)) + 
  geom_boxplot(outlier.shape = NA, notch = F)+
  geom_point(aes(fill = Mowing_delay,size=Mowing_frequency), alpha=0.7, pch=21,
             position=position_jitter(width = 0.05, height = 0)) +
  scale_fill_manual(values=c("blue", "red", "green"))+
  theme_bw()+
  labs(x ="Mowing delay", y="Species richness", fill="Mowing frequency", size="Mowing frequency")+
  theme(axis.text = element_text(size=10), axis.title=element_text(size=15)) +
  theme(axis.title.x=element_text(vjust=-0.1), axis.title.y=element_text(vjust=2)) +
  
  geom_text(data=emmeans_m_SR_3b_Mowing_delay,aes(x=Mowing_delay, y=c(85, 78, 87),
                                          label=emmeans_m_SR_3b_Mowing_delay$.group),vjust=0.5, hjust=0, 
            size=4, col="black" , position=position_dodge(0))



# Grazing_season  ----


m_SR_4 <- glmer(Plant_SR_vascular ~  
                #  Grazing_int_log +
                 # Grazer_type +
                  #Grazing_type +
                  Grazing_season +
                #  poly(Mowing_frequency, 2)  + # 
                  # Mowing_delay + 
                  (1|Farm), family = "poisson", data = Dat) 


plot(m_SR_4)
qqnorm(resid(m_SR_4))
qqline(resid(m_SR_4))

Anova(m_SR_4)
summary(m_SR_4)
check_overdispersion(m_SR_4) #  overdispersion detected
check_collinearity(m_SR_4)



m_SR_4b <- glmer.nb(Plant_SR_vascular ~ 
                     # Grazing_int_log +
                     # Grazer_type +
                     # Grazing_type +
                      Grazing_season +
                     # poly(Mowing_frequency, 2)  +
                      (1|Farm),
                    data = Dat) 

Anova(m_SR_4b)

check_convergence(m_SR_4b)
check_collinearity(m_SR_4b)


emmeansm_SR_4b_Grazing_season<- cld(emmeans(m_SR_3b, list(pairwise ~ Grazing_season)), 
                                   Letters = letters, test="Tukey") %>% arrange(Grazing_season)
emmeans_m_SR_4b_Grazing_season


ggplot(Dat, aes( Grazing_season, Plant_SR_vascular)) + 
  geom_boxplot(outlier.shape = NA, notch = F)+
  geom_point(aes(fill = Grazer_type,size=Grazing_intensity_A), 
             alpha=0.7,  pch=21,
             position=position_jitter(width = 0.05, height = 0)) +
  scale_fill_viridis(discrete=TRUE, option = "D")+
  theme_bw()+
  labs(x ="Grazing season", y="Species richness", 
       fill="Grazer type", size="Grazing intensity") +
  theme(axis.text = element_text(size=10), axis.title=element_text(size=15)) +
  theme(axis.title.x=element_text(vjust=-0.1), axis.title.y=element_text(vjust=2)) 




ggplot(Dat, aes( Grazing_season, Plant_SR_vascular)) + 
  geom_boxplot(outlier.shape = NA, notch = F)+
  geom_point(aes(fill = Grazer_type,size=Grazing_intensity_A,
                 shape=Grazing_type), alpha=0.7, # pch=21,
             position=position_jitter(width = 0.05, height = 0)) +
  scale_fill_viridis(discrete=TRUE, option = "D")+
  scale_shape_manual(values = c(21, 22, 25)) +
  theme_bw()+
  labs(x ="Grazing season", y="Species richness", 
       fill="Grazer type", size="Grazing intensity" ,
     shape="Grazing type") +
  theme(axis.text = element_text(size=10), axis.title=element_text(size=15)) +
  theme(axis.title.x=element_text(vjust=-0.1), axis.title.y=element_text(vjust=2)) 



# Grazing_type  ----

ggplot(Dat, aes( Grazing_type, Plant_SR_vascular)) + 
  geom_boxplot(outlier.shape = NA, notch = F)+
  geom_point(aes(fill = Grazer_type,size=Grazing_intensity_A),  pch=21, alpha=0.7,
             position=position_jitter(width = 0.05, height = 0)) +
  scale_fill_viridis(discrete=TRUE, option = "D")+
  theme_bw()+
  labs(x ="Grazing type", y="Species richness", fill="Grazer type", size="Grazing intensity")+
  theme(axis.text = element_text(size=10), axis.title=element_text(size=15)) +
  theme(axis.title.x=element_text(vjust=-0.1), axis.title.y=element_text(vjust=2)) 



# Cleaning, ploughing, crops  ------

Dat %>% 
  group_by(Grazing_type, Corralling) %>% 
  count()

m_SR_5 <- glmer(Plant_SR_vascular ~  
                  Grazing_int_log +
                  poly(Mowing_frequency, 2)  +
                 log(Last_ploughing) + # Crops_planted +
                  Cleaning +  # Shrub_tree_removal +
                  Litter_removal +
                  Moss_removal + 
                  Anthill_leveling + # Molehill_leveling +  # (data same for both)
                  Burning + Corralling +
                  (1|Farm), family = "poisson", data = Dat) 

plot(m_SR_5)
qqnorm(resid(m_SR_5))
qqline(resid(m_SR_5))


Anova(m_SR_5)
check_collinearity(m_SR_5)
check_overdispersion(m_SR_5) 
check_convergence(m_SR_5) 

plot_model(m_SR_5, type = "pred", terms = c("Last_ploughing"), show.data=T, dot.alpha=0.3, title="") 


plot_model(m_SR_5, type = "pred", terms = c("Litter_removal"), show.data=T, dot.alpha=0.3, title="") 
plot_model(m_SR_5, type = "pred", terms = c("Corralling"), show.data=T, dot.alpha=0.3, title="") 

m_SR_2_mowing_pred <- get_model_data(m_SR_2,type = "pred", terms="Mowing_frequency[0:3, by=0.1]")

ggplot(m_SR_2_mowing_pred, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1)+
  geom_point(data=Dat, aes(Mowing_frequency, Plant_SR_vascular, fill = Mowing_delay), 
             size=3, alpha=0.5, pch=21,
             position=position_jitter(width=0.05, height=0))+
  scale_fill_manual(values=c("blue", "red", "green"))+
  theme(axis.title.x=element_text(vjust=-0.1), axis.title.y=element_text(vjust=2)) +
  labs(y="Species richness", x='Mowing frequency', fill="Mowing delay")+
  geom_line(linetype=1, linewidth=1) 

















m1_SR_field <- glmer(Plant_SR_vascular ~   
                       sqrt(Grazing_intensity_A +1)  +  Grazer_diversity +
                       Mowing_frequency +  # Mowing_delay + #
                       #Mowing_method + 
                       Excrements_cover +  
                       Manuring_freq +
                       (1|Farm),  family = "poisson", 
                data = Dat) 

check_convergence(m1_SR_field)
check_collinearity(m1_SR_field)

check_overdispersion(m1_SR_field)

Anova(m1_SR_field)
summary(m1_SR_field)



m2_SR_field <- glmer.nb(Plant_SR_vascular ~   
                       Abund_D_E_exper +   SR_D_E_exper +
                       sqrt(Grazing_intensity_A +1)  +  Grazer_diversity +
                       Mowing_frequency +  # Mowing_delay + #
                       #Mowing_method + 
                         Excrements_cover +
                         Manuring_freq + (1|Farm),
                     data = Dat) 

check_convergence(m2_SR_field)

Anova(m2_SR_field)
summary(m2_SR_field)


m3_SR_field <- glm.nb(Plant_SR_vascular ~   
                        Abund_D_E_exper + SR_D_E_exper +
                        sqrt(Grazing_intensity_A +1)  +  Grazer_diversity +
                        Mowing_frequency +  # Mowing_delay + #
                        #Mowing_method + 
                        Excrements_cover +
                        Manuring_freq, 
                     data = Dat)

anova(m2_SR_field, m3_SR_field)

summary(m3_SR_field)
Anova(m3_SR_field)

# check plots
plot_model(m2_SR_field,type = "pred", terms = c("Grazing_intensity_A"), show.data=T, dot.alpha=0.3, title="") 
plot_model(m2_SR_field,type = "pred", terms = c("SR_D_E_exper"), show.data=T, dot.alpha=0.3, title="") 
plot_model(m2_SR_field,type = "pred", terms = c("Abund_D_E_exper"), show.data=T, dot.alpha=0.3, title="")
plot_model(m2_SR_field,type = "pred", terms = c("Excrements_cover"), show.data=T, dot.alpha=0.3, title="")
plot_model(m2_SR_field,type = "pred", terms = c("Grazer_diversity"), show.data=T, dot.alpha=0.3, title="")
plot_model(m2_SR_field,type = "pred", terms = c("Manuring_freq"), show.data=T, dot.alpha=0.3, title="")
plot_model(m2_SR_field,type = "pred", terms = c("Mowing_frequency"), show.data=T, dot.alpha=0.3, title="")
plot_model(m2_SR_field,type = "pred", terms = c("Mowing_delay"), show.data=T, dot.alpha=0.3, title="")

### mowing ----

m1_SR_field_pred <- get_model_data(m2_SR_field,type = "pred", terms="Mowing_frequency[0:3, by=0.1]")

ggplot(m1_SR_field_pred, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1)+
  geom_point(data=Dat,
             aes(Mowing_frequency, Plant_SR_vascular, fill = Mowing_delay), 
             size=3, alpha=0.5, pch=21,
             position=position_jitter(width=0.05, height=0))+
  scale_fill_manual(values=c("blue", "red", "green"))+
  labs(y="Species richness", x='Mowing frequency', fill="Mowing delay")+
  geom_line(linetype=1, linewidth=1) 


# Grazing intensity

m1_graz_pred <- get_model_data(m2_SR_field,type = "pred", terms="Grazing_intensity_A[0:1000, by=0.1]")

ggplot(m1_graz_pred, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1)+
  geom_point(data=Dat %>% 
               mutate(Grazer_type=str_replace_all(Grazer_type, "_", ", ")),
             aes(Grazing_intensity_A, Plant_SR_vascular, fill = Grazer_type), 
             size=3, alpha=0.7, pch=21)+
  scale_fill_viridis(discrete=TRUE, option = "D")+
  labs(y="Species richness", x='Grazing intensity', fill="Grazer type")+
  geom_line(linetype=1, linewidth=1) 


m2_SR_field <- glmer.nb(Plant_SR_vascular ~   
                          Abund_D_E_exper +   SR_D_E_exper +
                        #  sqrt(Grazing_intensity_A +1)  +  Grazer_diversity +
                         # Mowing_frequency +  # Mowing_delay + #
                          #Mowing_method + 
                         # Excrements_cover +
                        #  Manuring_freq + 
                          (1|Farm),
                        data = Dat) 

check_convergence(m2_SR_field)

Anova(m2_SR_field)

# Abundance (experiment)
m1_abun_pred <- get_model_data(m2_SR_field,type = "pred", terms="Abund_D_E_exper[1:50, by=0.1]")

ggplot(m1_abun_pred, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1)+
  geom_point(data=Dat %>% 
               mutate(Grazer_type=str_replace_all(Grazer_type, "_", ", ")),
             aes(Abund_D_E_exper, Plant_SR_vascular, fill = Grazer_type, size=Excrements_cover), 
             alpha=0.7, pch=21)+
  scale_fill_viridis(discrete=TRUE, option = "D")+
  labs(y="Species richness", x='Abundance of seedlings', fill="Grazer type", size="Excrements cover")+
  geom_line(linetype=1, linewidth=1) 


# SR (experiment)
m1_sr_pred <- get_model_data(m2_SR_field,type = "pred", terms="SR_D_E_exper[0:9, by=0.1]")

ggplot(m1_sr_pred, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1)+
  geom_point(data=Dat %>% 
               mutate(Grazer_type=str_replace_all(Grazer_type, "_", ", ")),
             aes(SR_D_E_exper, Plant_SR_vascular, fill = Grazer_type, size=Excrements_cover), 
             alpha=0.7, pch=21)+
  scale_fill_viridis(discrete=TRUE, option = "D")+
  labs(y="Species richness", x='Species number of seedlings', fill="Grazer type", size="Excrements cover")+
  geom_line(linetype=1, linewidth=1) 


# Excrement cover
m1_Excrem_pred <- get_model_data(m2_SR_field,type = "pred", terms="Excrements_cover[0:4, by=0.1]")

ggplot(m1_Excrem_pred, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1)+
  geom_point(data=Dat %>% 
               mutate(Grazer_type=str_replace_all(Grazer_type, "_", ", ")),
             aes(Excrements_cover, Plant_SR_vascular, fill = Grazer_type), 
             size=3, alpha=0.7, pch=21)+
  scale_fill_viridis(discrete=TRUE, option = "D")+
  labs(y="Species richness", x='Excrement cover', fill="Grazer type")+
  geom_line(linetype=1, linewidth=1) 



# Manuring frequency

m1_Manur_pred <- get_model_data(m2_SR_field,type = "pred", terms="Manuring_freq[0:1, by=0.1]")

ggplot(m1_Manur_pred, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1)+
  geom_point(data=Dat,
             aes(Manuring_freq, Plant_SR_vascular, fill = Dung_cover_log), 
             size=3, alpha=0.7, pch=21,
             position = position_jitter(width = 0.03, height=0))+
  scale_fill_viridis(option = "C")+
  labs(y="Species richness", x='Manuring frequency', fill="Dung cover (log)")+
  geom_line(linetype=1, linewidth=1) 



# Grazer_diversity

m1_Grazer_div_pred <- get_model_data(m2_SR_field,type = "pred", terms="Grazer_diversity[1:4, by=0.1]")

ggplot(m1_Grazer_div_pred, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1)+
  geom_point(data=Dat,
             aes(Grazer_diversity, Plant_SR_vascular, 
             fill = Grazing_intensity_A, size=Grazing_intensity_A),
              alpha=0.7, pch=21,
             position = position_jitter(width = 0.03, height=0))+
  scale_fill_viridis(option = "C")+
 labs(y="Species richness", x='Livestock diversity', fill="Grazing intensity", size="Grazing intensity")+
  geom_line(linetype=1, linewidth=1) 

  
### Grazer_type ----
# fit new model with other variables
m5_SR_field <- glm.nb(Plant_SR_vascular ~   
                       Abund_D_E_exper + 
                         SR_D_E_exper +
                      # sqrt(Grazing_intensity_A +1)  +   
                     Grazer_type +
                       # Mowing_frequency + #  
                         Mowing_delay , #+ 
                       # Mowing_method + 
                      # Excrements_cover +  
                      # Manuring_freq ,
                     data = Dat) 

check_collinearity(m5_SR_field)

check_overdispersion(m5_SR_field)

Anova(m5_SR_field)
summary(m5_SR_field)



emmeans_m5_SR_field <- cld(emmeans(m5_SR_field, list(pairwise ~ Grazer_type)), 
                            Letters = letters, test="Tukey") %>% arrange(Grazer_type)
emmeans_m5_SR_field

ggplot(Dat, aes(Grazer_type, Plant_SR_vascular)) + 
  geom_boxplot(outlier.shape = NA, notch = F)+
  geom_point(aes(fill = Grazing_intensity_A,size=Grazing_intensity_A), alpha=0.7, pch=21,
             position=position_jitter(width = 0.05, height = 0)) +
  scale_fill_viridis(option = "C")+
  theme_minimal()+
  labs(x ="Grazer type used for grazing in the parcel", y="Species richness", fill="Grazing intensity", size="Grazing intensity")+
  theme(axis.text = element_text(size=10), axis.title=element_text(size=12)) +
  
  geom_text(data=emmeans_m5_SR_field,aes(x=Grazer_type, y=85,
                                          label=emmeans_m5_SR_field$.group),vjust=0.5, hjust=0, 
            size=4, col="black" , position=position_dodge(0))


# mowing delay

emmeans_m5_mow_field <- cld(emmeans(m5_SR_field, list(pairwise ~ Mowing_delay)), 
                           Letters = letters, test="Tukey") %>% arrange(Mowing_delay)
emmeans_m5_mow_field

ggplot(Dat, aes(Mowing_delay, Plant_SR_vascular)) + 
  geom_boxplot(outlier.shape = NA, notch = F)+
  geom_point(aes(fill = Mowing_frequency,size=Mowing_frequency), alpha=0.7, pch=21,
             position=position_jitter(width = 0.05, height = 0)) +
  scale_fill_viridis(option = "C")+
  theme_minimal()+
  labs(x ="Mowing delay", y="Species richness", fill="Mowing frequency", size="Mowing frequency")+
  theme(axis.text = element_text(size=10), axis.title=element_text(size=12)) +
  
  geom_text(data=emmeans_m5_mow_field,aes(x=Mowing_delay, y=86,
                                         label=emmeans_m5_mow_field$.group),vjust=0.5, hjust=0, 
            size=4, col="black" , position=position_dodge(0))


# mowing method

m6_SR_field <- glm.nb(Plant_SR_vascular ~   
                        Abund_D_E_exper + 
                        SR_D_E_exper +
                        # sqrt(Grazing_intensity_A +1)  +   
                        Grazer_type +
                        # Mowing_frequency + #  
                        Mowing_delay + 
                       Mowing_method , 
                      # Excrements_cover +  
                      # Manuring_freq ,
                      data = Dat) 

check_collinearity(m5_SR_field)

check_overdispersion(m5_SR_field)

Anova(m5_SR_field)
summary(m5_SR_field)

### check R2 ----
# R2 for the entire model 
# R2m and R2c are marginal (for fixed predictors) and 
## conditional (for fixed and random predictors) coefficients of determination
MuMIn::r.squaredGLMM(m_SR_field)

# partial R2
r2glmm::r2beta(m_SR_field,  partial = T, method = 'sgv', data=Dat %>% filter(!Parcel_name=="Brade_1", !is.na(Abund_D_E_exper)))

# check in piecewiseSEM
summary(psem(m_SR_field))

# grazing type:


  ggplot(Dat, aes(Grazing_type, Plant_SR_vascular)) + 
  geom_boxplot(outlier.shape = NA, notch = F) +
  geom_point(alpha=0.7, pch=21,
             position=position_jitter(width = 0.05, height = 0)) +
 # scale_fill_viridis(option = "C")+
  # theme_minimal()+
  labs(x ="Grazing type", y="Species richness")+
  theme(axis.text = element_text(size=10), axis.title=element_text(size=12)) 




 
# END --------------

## Replace NAs with "0" in experiment data ----
m3 <- glmer.nb (Plant_SR_vascular ~ Grazing_intensity_A + SR_D_E_exper +
                  Abund_D_E_exper + (1|Farm), 
                data = Dat %>% filter(!Parcel_name=="Brade_1") %>% 
                  mutate(SR_D_E_exper=case_when(is.na(SR_D_E_exper) ~ 0, .default=SR_D_E_exper),
                         Abund_D_E_exper=case_when(is.na(Abund_D_E_exper) ~ 0, .default=Abund_D_E_exper))
                )

# check plots
plot_model(m3,type = "pred", terms = c("Grazing_intensity_A"), show.data=T, dot.alpha=0.3, title="") +
plot_model(m3,type = "pred", terms = c("SR_D_E_exper"), show.data=T, dot.alpha=0.3, title="") +
plot_model(m3,type = "pred", terms = c("Abund_D_E_exper"), show.data=T, dot.alpha=0.3, title="")

Anova(m3)
summary(m3)

# Experiment data ----
# check correlation between SR and Abundance in Experiment data

m4 <- glmer (SR_D_E_exper ~ Abund_D_E_exper + (1|Farm), family = "poisson", data = Dat %>% 
               filter(!Parcel_name=="Brade_1")) # strong outlier in SR_D_E_exper and Abund_D_E_exper

check_convergence(m1)

# check model
plot(m4)
qqnorm(resid(m4))
qqline(resid(m4))

# check multicolinearity
# check_collinearity(m4)

# check overdispersion
sum(residuals(m4, type = "pearson")^2) / df.residual(m4)
# no overdispersion

# check plots
plot_model(m4,type = "pred", terms = c("Abund_D_E_exper"), show.data=T, dot.alpha=0.3, title="")

Anova(m4)

# Replace NAs with "0" in experiment data 
m4b <- glmer (SR_D_E_exper ~ Abund_D_E_exper + (1|Farm), family = "poisson", 
              data = Dat %>% 
               filter(!Parcel_name=="Brade_1") %>% 
                        mutate(SR_D_E_exper=case_when(is.na(SR_D_E_exper) ~ 0, .default=SR_D_E_exper),
                               Abund_D_E_exper=case_when(is.na(Abund_D_E_exper) ~ 0, .default=Abund_D_E_exper))
               )

# check plots
plot_model(m4b,type = "pred", terms = c("Abund_D_E_exper"), show.data=T, dot.alpha=0.3, title="")
              
Anova(m4b)





