# Purpose: GLMMs for Species Richness (SR)

rm(list = ls())

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


# set theme for plots
 set_theme(base = theme_bw(),
          axis.textsize.x = 0.9, axis.textsize.y = 0.9, axis.textcolor = "black",
          axis.title.color = "black", axis.title.size = 1.2,
          geom.linetype = 1) #legend.pos = "None", 



# data  ----

Soil_PC <- read_csv("data/soil_NPK_PCA.csv")  %>% 
  mutate(soil_NPK=-1*soil_NPK_PC)  # reverse the NMDS scores to make it positively correlated with the nutrients


Data <- read_csv("data/Panoara_Dat.csv") %>%
  mutate(Grazer_type=str_replace_all(Grazer_type, "_", ", ")) %>% 
  mutate(Grazer_type=fct_relevel(Grazer_type,c("cow","sheep, goat","mixed"))) %>%
  mutate(Mowing_delay=fct_relevel(Mowing_delay,c("no mowing","June","July-August"))) %>%
  arrange(Mowing_delay) %>% 
  mutate(Cleaning=fct_relevel(Cleaning,c("no","irregular","regular"))) %>%
  mutate(Dung_cover_log=log(Dung_cover+1)) %>% 
  mutate(SR_D_E_exper=case_when(is.na(SR_D_E_exper) ~ 0, .default=SR_D_E_exper),
         Abund_D_E_exper=case_when(is.na(Abund_D_E_exper) ~ 0, .default=Abund_D_E_exper)) %>% 
  mutate(Grazing_int_log = log1p(Grazing_intensity_A),
         Grazing_legacy =log1p(Grazing_intensity_B),
         Manuring_freq_log = log1p(Manuring_freq),
         Corralling=case_when(Corralling=="0" ~ "no", 
                              Corralling=="1" ~ "yes")) %>% 
  mutate(humus_log=log1p(humus)) %>% 
  mutate(Ploughing=1/Last_ploughing)
# filter(!Parcel_name=="Brade_1") # extreme outlier


Dat <- Data %>% 
  left_join(Soil_PC, by="Parcel_name")


str(Dat)
names(Dat)


Dat <- Dat %>% 
  mutate(zz = Plant_SR_vascular) # for predicting from the glmmPQL models





# GLMM models -----
##(1) Habitat, mang. stability , soil----
m_SR_1 <- glmer(Plant_SR_vascular ~  habitat_corrected + 
                  Management_stability +
                  humus_log + soil_NPK+
                (1|Farm), family = "poisson", data = Dat) 

check_convergence(m_SR_1)
check_collinearity(m_SR_1)
summary(m_SR_1)
Anova(m_SR_1)

check_overdispersion(m_SR_1) # Overdispersion detected.



m_SR_6 <- glmer(Plant_SR_vascular ~  
                  # Grazing_int_log +
                  Grazing_legacy + 
                  poly(Mowing_frequency, 2)  +
                  (1|Farm), family = "poisson", data = Dat) 

plot(m_SR_6)
qqnorm(resid(m_SR_6))
qqline(resid(m_SR_6))

Anova(m_SR_6)
check_collinearity(m_SR_6)
check_overdispersion(m_SR_6)


m_SR_6b <- glmer.nb(Plant_SR_vascular ~  
                      # Grazing_int_log +
                      Grazing_legacy + 
                      poly(Mowing_frequency, 2)  +
                      (1|Farm), data = Dat) 




# Change family to negative binomial
m_SR_1b <- glmer.nb(Plant_SR_vascular ~ 
            habitat_corrected +
             Management_stability +
              humus_log + soil_NPK+
              Grazing_legacy +
           (1|Farm), data = Dat) 

PearsonResiduals <- resid(m_SR_1b, type = "pearson")
n_cases <- nrow(Dat) # extract number of samples
n_predictors <- length(fixef(m_SR_1b)) + 1 # extract number of estimates (plus intercept)
Overdispersion <- sum(PearsonResiduals^2) / (n_cases-n_predictors) # calculate overdispersion
Overdispersion # The overdispersion has decreased and is no longer an issue.

plot(m_SR_1b) #  the plot exhibits a slight funnel shape, but not drastically, and thus indicates heteroscedasticity.
qqnorm(resid(m_SR_1b))
qqline(resid(m_SR_1b))

check_convergence(m_SR_1b)
check_collinearity(m_SR_1b)

Anova(m_SR_1b)

#### > habitat ----
emmeans_m1_habitat <- cld(emmeans(m_SR_1b, list(pairwise ~ habitat_corrected)), 
                          Letters = letters) %>% arrange(habitat_corrected)
emmeans_m1_habitat

ggplot(Dat, aes(habitat_corrected, Plant_SR_vascular)) + 
  geom_boxplot(outlier.shape = NA, notch = F)+
  geom_point(aes(fill = Grazer_type, size=Grazing_intensity_A),
             alpha=0.7, pch=21, # size=4, #fill="gray", 
             position=position_jitter(width = 0.1, height = 0)) +
  scale_fill_viridis(discrete=TRUE, option = "D")+
 # theme_bw()+
  labs(x ="Management type", y="Species richness", fill="Grazer type", 
       size="Grazing intensity")+
  theme(axis.text = element_text(size=10), axis.title=element_text(size=12)) +
  theme(axis.title.x=element_text(vjust=-0.7), 
        axis.title.y=element_text(vjust=1.5),
        legend.key=element_blank()) +
  geom_text(data=emmeans_m1_habitat,aes(x=habitat_corrected, y=c(89, 87),
                                        label=emmeans_m1_habitat$.group),vjust=0.5, hjust=0, 
            size=4, col="black" , position=position_dodge(0))



 ### > Management_stability ----
m_SR_2_mowing_pred <- get_model_data(m_SR_1b,type = "pred", terms="Management_stability[0:50, by=0.001]")

ggplot(m_SR_2_mowing_pred, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1)+
  geom_point(data=Dat, aes(Management_stability, Plant_SR_vascular, fill = habitat_corrected), 
             size=3, alpha=1, pch=21,
             position=position_jitter(width=0.05, height=0))+
  scale_fill_manual(values=c("forestgreen", "orange"))+
  theme(axis.title.x=element_text(vjust=-0.1), axis.title.y=element_text(vjust=2)) +
  labs(y="Species richness", x='Management stability', fill="Management type")+
  geom_line(linetype=1, linewidth=1) +
#  theme_bw()+
  theme(axis.text = element_text(size=10), 
        axis.title=element_text(size=12)) +
  theme(axis.title.x=element_text(vjust=-0.7), 
        axis.title.y=element_text(vjust=1.5),
       # legend.title = element_text(face="plain"),
        legend.key=element_blank()) 

max(Dat$humus_log)

### > Soil humus -----

m_SR_1b_humus_pred <- get_model_data(m_SR_1b,type = "pred", 
                                      terms="humus_log[1.87:3.63, by=0.001]")

ggplot(m_SR_1b_humus_pred, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1)+
  geom_point(data=Dat,
             aes(humus_log, Plant_SR_vascular), 
             size=3, alpha=0.7, pch=21, fill ="#64ABCE")+
  # scale_fill_viridis(discrete=TRUE, option = "D")+
  labs(y="Species richness", x='Humus', fill="Grazer type")+
#  theme_bw()+
  geom_line(linetype=1, linewidth=1) 


### > Soil NPK -----

m_SR_1b_soil_NPK_pred <- get_model_data(m_SR_1b,type = "pred", 
                                     terms="soil_NPK")

ggplot(m_SR_1b_soil_NPK_pred, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1)+
  geom_point(data=Dat,
             aes(soil_NPK, Plant_SR_vascular), 
             size=3, alpha=0.7, pch=21, fill ="#64ABCE")+
  # scale_fill_viridis(discrete=TRUE, option = "D")+
  labs(y="Species richness", x='Soil NPK', fill="Grazer type")+
 # theme_bw()+
  geom_line(linetype=1, linewidth=1) 



### > Legacy effect ----
m_SR_1_legac <- get_model_data(m_SR_1b, type = "pred", 
                                     terms="Grazing_legacy[0:4, by=0.001]")

ggplot(m_SR_1_legac, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1)+
  geom_point(data=Dat,
             aes(Grazing_legacy, Plant_SR_vascular, fill=Grazer_type), 
             size=3, alpha=0.7, pch=21)+
  scale_fill_viridis(discrete=TRUE, option = "D")+
  labs(y="Species richness", x='Grazing legacy', fill= "Grazer type")+
 # theme_bw()+
  geom_line(linetype=1, linewidth=1) +
  theme(legend.position="none")



### R2  ----
# R2 for the entire model
# R2m and R2c are marginal (for fixed predictors) and 
## conditional (for fixed and random predictors) coefficients of determination
R2_m_SR_1b=MuMIn::r.squaredGLMM(m_SR_1b)
R2_m_SR_1b

# partial R2 
R2_part_m_SR_1b=r2glmm::r2beta(m_SR_1b,  partial = T, method = 'sgv', data=Dat)
R2_part_m_SR_1b
plot(R2_part_m_SR_1b)

# (2) Grazing, Mowing----

m_SR_2 <- glmer(Plant_SR_vascular ~  
                  Grazing_int_log +
                  #  poly(Grazing_int_log , 2)  +
                  poly(Mowing_frequency, 2)  +
                 (1|Farm), family = "poisson", data = Dat) 

plot(c)
qqnorm(resid(m_SR_2))
qqline(resid(m_SR_2))

Anova(m_SR_2)
check_collinearity(m_SR_2)

plot_model(m_SR_2, type = "pred", terms = c("Grazing_int_log"), show.data=T, dot.alpha=0.3, title="") 
plot_model(m_SR_2, type = "pred", terms = c("Mowing_frequency"), show.data=T, dot.alpha=0.3, title="") 


# check for hump-shape for grazing:
m_SR_2b <- glmer(Plant_SR_vascular ~  
                   #  Grazing_int_log +
                   poly(Grazing_int_log , 2)  +
                   poly(Mowing_frequency, 2)  +
                  (1|Farm), family = "poisson", data = Dat) 

anova(m_SR_2, m_SR_2b) # selects m_SR_2, against the hump-shape for grazing
Anova(m_SR_2)

check_overdispersion(m_SR_2) # Slight overdispersion is fine


## -> Mowing_frequency -----
m_SR_2_mowing_pred <- get_model_data(m_SR_2,type = "pred", terms="Mowing_frequency[0:3, by=0.1]")

ggplot(m_SR_2_mowing_pred, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1)+
  geom_point(data=Dat, aes(Mowing_frequency, Plant_SR_vascular, fill = Mowing_delay), 
             size=3, alpha=0.5, pch=24,
             position=position_jitter(width=0.05, height=0))+
  scale_fill_manual(values=c("blue", "red", "green"))+
#  theme(axis.title.x=element_text(vjust=-0.1), axis.title.y=element_text(vjust=2)) +
  labs(y="Species richness", x='Mowing frequency', fill="Mowing delay")+
  theme(legend.key=element_blank(),
        legend.position = "none") +
  geom_line(linetype=1, linewidth=1) #+ theme(legend.position="none")



## -> Grazing intensity -----

m_SR_2_graz_pred <- get_model_data(m_SR_2,type = "pred", 
                                   terms="Grazing_int_log[2.39:6.9, by=0.001]")

ggplot(m_SR_2_graz_pred, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1)+
  geom_point(data=Dat,
             aes(Grazing_int_log, Plant_SR_vascular, fill = Grazer_type), 
             size=3, alpha=0.7, pch=21)+
  scale_fill_viridis(discrete=TRUE, option = "D")+
  labs(y="Species richness", x='Grazing intensity', fill="Grazer type")+
  theme(legend.key=element_blank(),
        legend.position = "none") +
  geom_line(linetype=1, linewidth=1) # +  theme(legend.position="none")


### R2  ----
# R2 for the entire model 
# R2m and R2c are marginal (for fixed predictors) and 
## conditional (for fixed and random predictors) coefficients of determination
R2_m_SR_2=MuMIn::r.squaredGLMM(m_SR_2)
R2_m_SR_2

# partial R2 
R2_part_m_SR_2=r2glmm::r2beta(m_SR_2,  partial = T, method = 'sgv', data=Dat)
R2_part_m_SR_2
plot(R2_part_m_SR_2)


#  (3) Grazer_type & Mowing_delay ----

m_SR_3 <- glmer(Plant_SR_vascular ~  
                  Grazer_type +
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
                      Grazer_type +
                      #  Grazing_type +
                      # poly(Mowing_frequency, 2)  + # 
                      Mowing_delay + 
                     # Mowing_method +
                      (1|Farm), data = Dat) 

check_convergence(m_SR_3b)
check_collinearity(m_SR_3b)

Anova(m_SR_3b)
summary(m_SR_3b)


# test Grazing_int_log as covariate

m_SR_3b_2 <- glmer.nb(Plant_SR_vascular ~ 
                  Grazing_int_log +
                      Grazer_type +
                 # poly(Mowing_frequency, 2) # The effect of mowing delay correlates strongly with the mowing frequency (not possible to fit both in the same model → high singularity) 
                      Mowing_delay + 
                      (1|Farm), data = Dat) 

check_convergence(m_SR_3b_2)
check_collinearity(m_SR_3b_2)

Anova(m_SR_3b_2)

## -> Grazer_type ----

emmeans_m_SR_3b_Grazer_type <- cld(emmeans(m_SR_3b, list(pairwise ~ Grazer_type)), 
                                   Letters = letters, test="Tukey") %>% arrange(Grazer_type)
emmeans_m_SR_3b_Grazer_type

ggplot(Dat, aes(Grazer_type, Plant_SR_vascular)) + 
  geom_boxplot(outlier.shape = NA, notch = F)+
  geom_point(aes(fill = Grazer_type,size=Grazing_intensity_A), alpha=0.7, pch=21,
             position=position_jitter(width = 0.05, height = 0)) +
  scale_fill_viridis(discrete=TRUE, option = "D")+
 # theme_bw()+
  labs(x ="Grazer type", y="Species richness", fill="Grazer type", size="Grazing intensity")+
  theme(axis.text.x = element_text(size=8.5), axis.text.y = element_text(size=9),
        axis.title=element_text(size=13)) +
  theme(axis.title.x=element_text(vjust=-0.1), axis.title.y=element_text(vjust=2),
        legend.key=element_blank(),
        legend.position = "none") +
    geom_text(data=emmeans_m_SR_3b_Grazer_type,aes(x=Grazer_type, y=c(87, 78, 83),
                                                 label=emmeans_m_SR_3b_Grazer_type$.group),vjust=0.5, hjust=0, 
          #  size=4, 
          col="black" , position=position_dodge(0))  # +   theme(legend.position="none")


## -> Mowing delay ----

emmeans_m_SR_3b_Mowing_delay <- cld(emmeans(m_SR_3b, list(pairwise ~ Mowing_delay)), 
                                    Letters = letters, test="Tukey") %>% arrange(Mowing_delay)
emmeans_m_SR_3b_Mowing_delay

ggplot(Dat, aes(Mowing_delay, Plant_SR_vascular)) + 
  geom_boxplot(outlier.shape = NA, notch = F)+
  geom_point(aes(fill = Mowing_delay,size=Mowing_frequency), alpha=0.5, pch=24,
             position=position_jitter(width = 0.05, height = 0)) +
  scale_fill_manual(values=c("blue", "red", "green"))+
#  theme_bw()+
  labs(x ="Mowing delay", y="Species richness", fill="Mowing frequency", size="Mowing frequency")+
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=9),
        axis.title=element_text(size=13)) +
  theme(axis.title.x=element_text(vjust=-0.1), axis.title.y=element_text(vjust=2),
        legend.key=element_blank()) +
  geom_text(data=emmeans_m_SR_3b_Mowing_delay,aes(x=Mowing_delay, y=c(85, 81, 89),
                                                  label=emmeans_m_SR_3b_Mowing_delay$.group),vjust=0.5, hjust=0, 
          col="black" , position=position_dodge(0)) # +  theme(legend.position="none")



### R2 ----
# R2 for the entire model 
# R2m and R2c are marginal (for fixed predictors) and 
## conditional (for fixed and random predictors) coefficients of determination
R2_m_SR_3b=MuMIn::r.squaredGLMM(m_SR_3b)
R2_m_SR_3b

# partial R2
R2_part_m_SR_3b= r2glmm::r2beta(m_SR_3b,  partial = T, method = 'sgv', data=Dat)
R2_part_m_SR_3b
plot(R2_part_m_SR_3b)

# (4) Grazing_type, Grazing_season  ----


m_SR_4 <- glmer(Plant_SR_vascular ~  
                  Grazing_int_log +
                  Grazing_type +
                  Grazing_season +
                 # poly(Mowing_frequency, 2)  + # 
                  (1|Farm), family = "poisson", data = Dat) 


plot(m_SR_4)
qqnorm(resid(m_SR_4))
qqline(resid(m_SR_4))

Anova(m_SR_4)
summary(m_SR_4)
check_overdispersion(m_SR_4) #  overdispersion detected
check_collinearity(m_SR_4)



m_SR_4b <- glmer.nb(Plant_SR_vascular ~ 
                     Grazing_int_log +
                       Grazing_type +
                       Grazing_season +
                    #   poly(Mowing_frequency, 2)  +
                      (1|Farm),
                    data = Dat) 

check_convergence(m_SR_4b)
check_collinearity(m_SR_4b)

Anova(m_SR_4b)



## > Grazing_season ----
emmeansm_SR_4b_Grazing_season <- cld(emmeans(m_SR_4b, list(pairwise ~ Grazing_season)), 
                                    Letters = letters, test="Tukey") %>% arrange(Grazing_season)
emmeansm_SR_4b_Grazing_season


ggplot(Dat, aes( Grazing_season, Plant_SR_vascular)) + 
  geom_boxplot(outlier.shape = NA, notch = F)+
  geom_point(aes(fill = Grazer_type,size=Grazing_intensity_A), 
             alpha=0.7,  pch=21,
             position=position_jitter(width = 0.05, height = 0)) +
  scale_fill_viridis(discrete=TRUE, option = "D")+
 # theme_bw()+
  labs(x ="Grazing season", y="Species richness", 
       fill="Grazer type", size="Grazing intensity") +
  theme(axis.text.x = element_text(size=11), axis.text.y = element_text(size=10), 
        axis.title=element_text(size=15)) +
  theme(axis.title.x=element_text(vjust=-0.1), 
        axis.title.y=element_text(vjust=2), legend.key=element_blank()) +
  geom_text(data=emmeansm_SR_4b_Grazing_season,aes(x=Grazing_season, y=c(87, 76, 85),
                                                 label=emmeansm_SR_4b_Grazing_season$.group),vjust=0.5, hjust=0, 
         #   size=4, 
         col="black" , position=position_dodge(0))



ggplot(Dat, aes( Grazing_season, Plant_SR_vascular)) + 
  geom_boxplot(outlier.shape = NA, notch = F)+
  geom_point(aes(fill = Grazer_type,size=Grazing_intensity_A,
                 shape=Grazing_type), alpha=0.7, # pch=21,
             position=position_jitter(width = 0.05, height = 0)) +
  scale_fill_viridis(discrete=TRUE, option = "D")+
  scale_shape_manual(values = c(21, 22, 25)) +
 # theme_bw()+
  labs(x ="Grazing season", y="Species richness", 
       fill="Grazer type", size="Grazing intensity" ,
       shape="Grazing type") +
  theme(axis.text = element_text(size=11), axis.title=element_text(size=15), 
        legend.key=element_blank()) +
  theme(axis.title.x=element_text(vjust=-0.1), axis.title.y=element_text(vjust=2)) 



## > Grazing_type  ----

emmeansm_SR_4b_Grazing_type <- cld(emmeans(m_SR_4b, list(pairwise ~ Grazing_type)), 
                                     Letters = letters, test="Tukey") %>% arrange(Grazing_type)
emmeansm_SR_4b_Grazing_type

ggplot(Dat, aes( Grazing_type, Plant_SR_vascular)) + 
  geom_boxplot(outlier.shape = NA, notch = F)+
  geom_point(aes(fill = Grazer_type,size=Grazing_intensity_A),  pch=21, alpha=0.7,
             position=position_jitter(width = 0.05, height = 0)) +
  scale_fill_viridis(discrete=TRUE, option = "D")+
 
  # theme_bw()+
  labs(x ="Grazing type", y="Species richness", fill="Grazer type", size="Grazing intensity")+
  theme(axis.text.x = element_text(size=11), axis.text.y = element_text(size=10), 
        axis.title=element_text(size=15)) +
  theme(axis.title.x=element_text(vjust=-0.1), 
        axis.title.y=element_text(vjust=2), legend.key=element_blank()) +
  geom_text(data=emmeansm_SR_4b_Grazing_type,aes(x=Grazing_type, y=c(86, 87, 77),
                                                   label=emmeansm_SR_4b_Grazing_type$.group),vjust=0.5, hjust=0, 
           # size=4, 
           col="black" , position=position_dodge(0))

### R2 ----
# R2 for the entire model 
# R2m and R2c are marginal (for fixed predictors) and 
## conditional (for fixed and random predictors) coefficients of determination
R2_m_SR_4b=MuMIn::r.squaredGLMM(m_SR_4b)
R2_m_SR_4b

# partial R2
R2_part_m_SR_4b= r2glmm::r2beta(m_SR_4b,  partial = T, method = 'sgv', data=Dat)
R2_part_m_SR_4b
plot(R2_part_m_SR_4b)


# (5) Cleaning, manuring  ------

m_SR_5.1.a <- glmer(Plant_SR_vascular ~  
                   Ploughing + # Crops_planted + # correlates with Last_ploughing
                   Shrub_tree_removal + # Cleaning +  # Shrub_tree_removal  correlates with Anthill_leveling
                   Litter_removal +
                   Moss_removal + 
                   #  Anthill_leveling + #  Molehill_leveling +# data are the same as Anthill_leveling
                   Burning + 
                   Corralling +
                   Manuring_freq +
                   Cow_dung_applied +
                  (1|Farm),  
                  family = "poisson", data = Dat)  


Anova(m_SR_5.1.a)



check_collinearity(m_SR_5.1.a)
check_overdispersion(m_SR_5.1.a)

m_SR_5b_1 <- glmer.nb(Plant_SR_vascular ~  
                         Ploughing + # Crops_planted + # correlates with Last_ploughing
                         Shrub_tree_removal + # Cleaning +  # Shrub_tree_removal  correlates with Anthill_leveling
                         Litter_removal +
                         Moss_removal + 
                         #  Anthill_leveling + #  Molehill_leveling +# data are the same as Anthill_leveling
                         Burning + 
                         Corralling +
                         Manuring_freq +
                         Cow_dung_applied +
                      (1|Farm), data = Dat) 

Anova(m_SR_5b_1)

# test for the Anthill_leveling instead of Shrub_tree_removal

m_SR_5b_2 <- glmer.nb(Plant_SR_vascular ~  
                        Ploughing + # Crops_planted + # correlates with Last_ploughing
                      #  Shrub_tree_removal + # Cleaning +  # Shrub_tree_removal  correlates with Anthill_leveling
                        Litter_removal +
                        Moss_removal + 
                        Anthill_leveling + #  Molehill_leveling +# data are the same as Anthill_leveling
                        Burning + 
                        Corralling +
                        Manuring_freq +
                        Cow_dung_applied +
                      (1|Farm), data = Dat) 


Anova(m_SR_5b_2)




ranef(m_SR_5b_2) # random effects are not zeros
hist(ranef(m_SR_5b_2)$`Farm`[,1])

check_convergence(m_SR_5b_2) 
check_collinearity(m_SR_5b_2)
check_convergence(m_SR_5b_2)
check_overdispersion(m_SR_5b_2)  # very low overdispersion




### R2 ----
# R2 for the entire model 
# R2m and R2c are marginal (for fixed predictors) and 
## conditional (for fixed and random predictors) coefficients of determination
R2_m_SR_5b_1=MuMIn::r.squaredGLMM(m_SR_5b_1)
R2_m_SR_5b_1

# partial R2
R2_part_m_SR_5b_1= r2glmm::r2beta(m_SR_5b_1,  partial = T, method = 'sgv', data=Dat)
plot(R2_part_m_SR_5b_1)

R2_part_m_SR_5b_2= r2glmm::r2beta(m_SR_5b_2,  partial = T, method = 'sgv', data=Dat)
R2_part_m_SR_5b_2
as_tibble(R2_part_m_SR_5b_2) %>% filter(Effect=="Anthill_levelingyes")

### > Litter_removal ----

emmeans_m_SR_5_Litter <- cld(emmeans(m_SR_5b_1, list(pairwise ~ Litter_removal)), 
                                    Letters = letters, test="Tukey") %>% arrange(Litter_removal)
emmeans_m_SR_5_Litter

ggplot(Dat, aes(Litter_removal, Plant_SR_vascular)) + 
  geom_boxplot(outlier.shape = NA, notch = F)+
  geom_point(alpha=0.6, pch=21, size=3,
             position=position_jitter(width = 0.05, height = 0),
             fill ="#64ABCE") +
 # scale_fill_manual(values=c("blue", "red", "green"))+
#  theme_bw()+
  labs(x ="Litter removal", y="Species richness")+
  theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=9),
        axis.title=element_text(size=13)) +
  theme(axis.title.x=element_text(vjust=-0.1), axis.title.y=element_text(vjust=2)) +
  
  geom_text(data=emmeans_m_SR_5_Litter,aes(x=Litter_removal, y=c(86, 90),
                                                  label=emmeans_m_SR_5_Litter$.group),
            vjust=0.5, hjust=0, #size=4, 
            col="black" , position=position_dodge(0))




### >  Corralling ----

emmeans_m_SR_5_Corralling <- cld(emmeans(m_SR_5b_1, list(pairwise ~ Corralling)), 
                                  Letters = letters, test="Tukey") %>% arrange(Corralling)
emmeans_m_SR_5_Corralling

ggplot(Dat, aes(Corralling, Plant_SR_vascular)) + 
  geom_boxplot(outlier.shape = NA, notch = F)+
  geom_point(alpha=0.6, pch=21, size=3,
             position=position_jitter(width = 0.05, height = 0),
             fill ="#64ABCE") +
  # scale_fill_manual(values=c("blue", "red", "green"))+
#  theme_bw()+
  labs(x ="Corraling", y="Species richness")+
  theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=9),
        axis.title=element_text(size=13)) +
  theme(axis.title.x=element_text(vjust=-0.1), axis.title.y=element_text(vjust=2)) +
  
  geom_text(data=emmeans_m_SR_5_Corralling, aes(x=Corralling, y=c(90, 70),
                                                label=c("a", "b")),vjust=0.5, hjust=0, # size=4, 
           col="black" , position=position_dodge(0))




### > Shrub_tree_removal ----
emmeans_m_SR_5_shrub <- cld(emmeans(m_SR_5b_1, list(pairwise ~ Shrub_tree_removal)), 
                                     Letters = letters, test="Tukey") %>% arrange(Shrub_tree_removal)
emmeans_m_SR_5_shrub

ggplot(Dat, aes(Shrub_tree_removal, Plant_SR_vascular)) + 
  geom_boxplot(outlier.shape = NA, notch = F)+
  geom_point(alpha=0.6, pch=21, size=3,
             position=position_jitter(width = 0.05, height = 0),
             fill ="#64ABCE") +
  # scale_fill_manual(values=c("blue", "red", "green"))+
#  theme_bw()+
  labs(x ="Shrub/tree removal", y="Species richness")+
  theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=9),
        axis.title=element_text(size=13)) +
  theme(axis.title.x=element_text(vjust=-0.1), axis.title.y=element_text(vjust=2))+
  geom_text(data=emmeans_m_SR_5_shrub,aes(x=Shrub_tree_removal, y=c(70, 90),
                                                   label=c("a", "b")),
            vjust=0.5, hjust=0, # size=4, 
            col="black", position=position_dodge(0)) 

### > Anthill_leveling  #  Molehill_leveling ----

emmeans_m_SR_5_Anthill <- cld(emmeans(m_SR_5b_2, list(pairwise ~ Anthill_leveling)), 
                              Letters = letters, test="Tukey") %>% arrange(Anthill_leveling)
emmeans_m_SR_5_Anthill

ggplot(Dat, aes(Anthill_leveling, Plant_SR_vascular)) + 
  geom_boxplot(outlier.shape = NA, notch = F)+
  geom_point(alpha=0.6, pch=21, size=3, fill ="#64ABCE",
             position=position_jitter(width = 0.05, height = 0))  +
  # scale_fill_manual(values=c("blue", "red", "green"))+
 # theme_bw()+
  labs(x ="Ant/mole-hill leveling", y="Species richness")+
  theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=9),
        axis.title=element_text(size=13)) +
  theme(axis.title.x=element_text(vjust=-0.1), axis.title.y=element_text(vjust=2)) +
  geom_text(data=emmeans_m_SR_5_Anthill,aes(x=Anthill_leveling, y=c(75, 90),
                                            label=emmeans_m_SR_5_Anthill$.group),
            vjust=0.5, hjust=0, # size=4, 
           col="black" , position=position_dodge(0)) 


### >  Ploughing ----
m_SR_5_plough_pred <- get_model_data(m_SR_5b_1,type = "pred", terms="Ploughing[0:0.2, by=0.01]")

min(Dat$Ploughing)

ggplot(m_SR_5_plough_pred, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1)+
  geom_point(data=Dat, aes(Ploughing, Plant_SR_vascular), 
             size=3, alpha=0.5, pch=21, fill ="#64ABCE",
             position=position_jitter(width=0.01, height=0))+
  # scale_fill_manual(values=c("blue", "red", "green"))+
  theme(axis.title.x=element_text(vjust=-0.1), axis.title.y=element_text(vjust=2)) +
  labs(y="Species richness", x='Ploughing intensity')+
  # scale_x_continuous(breaks=seq(0, 60, 10), 
  #                   labels=paste0(c("0", "10","20", "30", "40", "50", ">50"))) +
  geom_line(linetype=1, linewidth=1) 

### > Moss_removal ----

emmeans_m_SR_5_Moss <- cld(emmeans(m_SR_5b_1, list(pairwise ~ Moss_removal)), 
                                Letters = letters, test="Tukey") %>% arrange(Moss_removal)
emmeans_m_SR_5_Moss

ggplot(Dat, aes(Moss_removal, Plant_SR_vascular)) + 
  geom_boxplot(outlier.shape = NA, notch = F)+
  geom_point(alpha=0.6, pch=21, size=3,
             position=position_jitter(width = 0.05, height = 0),
             fill ="#64ABCE") +
  # scale_fill_manual(values=c("blue", "red", "green"))+
 # theme_bw()+
  labs(x ="Moss removal", y="Species richness")+
  theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=9),
        axis.title=element_text(size=13)) +
  theme(axis.title.x=element_text(vjust=-0.1), axis.title.y=element_text(vjust=2)) +
  geom_text(data=emmeans_m_SR_5_Moss,aes(x=Moss_removal, y=c(89, 87),
                                          label=emmeans_m_SR_5_Moss$.group),
            vjust=0.5, hjust=0, col="black", position=position_dodge(0)) 




### > Burning ----

emmeans_m_SR_5_Burning <- cld(emmeans(m_SR_5b_1, list(pairwise ~ Burning)), 
                                    Letters = letters, test="Tukey") %>% arrange(Burning)

emmeans_m_SR_5_Burning

ggplot(Dat, aes(Burning, Plant_SR_vascular)) + 
  geom_boxplot(outlier.shape = NA, notch = F)+
  geom_point(alpha=0.6, pch=21, size=3,
             position=position_jitter(width = 0.05, height = 0),
             fill ="#64ABCE") +
  # scale_fill_manual(values=c("blue", "red", "green"))+
#  theme_bw()+
  labs(x ="Burning", y="Species richness")+
  theme(axis.text.x = element_text(size=13), axis.text.y = element_text(size=9),
        axis.title=element_text(size=13)) +
  theme(axis.title.x=element_text(vjust=-0.1), axis.title.y=element_text(vjust=2))+
  geom_text(data=emmeans_m_SR_5_Burning,aes(x=Burning, y=c(86, 89),
                                         label=emmeans_m_SR_5_Burning$.group),
            vjust=0.5, hjust=0, col="black", position=position_dodge(0)) 




### > Manuring_freq   ----
m_SR_5_Manuring_pred <- get_model_data(m_SR_5b_1,type = "pred", terms="Manuring_freq[0:1, by=0.0001]")

ggplot(m_SR_5_Manuring_pred, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1)+
  geom_point(data=Dat, aes(Manuring_freq, Plant_SR_vascular), 
             size=3, alpha=0.5, pch=21, fill ="#64ABCE",
             position=position_jitter(width=0.01, height=0))+
  # scale_fill_manual(values=c("blue", "red", "green"))+
  theme(axis.title.x=element_text(vjust=-0.1), axis.title.y=element_text(vjust=2)) +
  labs(y="Species richness", x='Manuring frequency')+
 # theme_bw()+
  geom_line(linetype=1, linewidth=1) 




### Cow_dung_applied ----


emmeans_m_SR_5_Cow_dung  <- cld(emmeans(m_SR_5b_1, list(pairwise ~ Cow_dung_applied)), 
                               Letters = letters, test="Tukey") %>% arrange(Cow_dung_applied)
emmeans_m_SR_5_Cow_dung

ggplot(Dat, aes(Cow_dung_applied, Plant_SR_vascular)) + 
  geom_boxplot(outlier.shape = NA, notch = F)+
  geom_point(alpha=0.6, pch=21, size=3,
             position=position_jitter(width = 0.05, height = 0),
             fill ="#64ABCE") +
  # scale_fill_manual(values=c("blue", "red", "green"))+
#  theme_bw()+
  labs(x ="Cow manure applied", y="Species richness")+
  theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=9),
        axis.title=element_text(size=13)) +
  theme(axis.title.x=element_text(vjust=-0.1), axis.title.y=element_text(vjust=2)) +
  geom_text(data=emmeans_m_SR_5_Cow_dung,aes(x=Cow_dung_applied, y=c(89, 80),
                                            label=emmeans_m_SR_5_Cow_dung$.group),
            vjust=0.5, hjust=0,  col="black", position=position_dodge(0)) 


# Table results -----

## > Anova results----

Anova_table <- as_tibble(Anova(m_SR_1b), rownames ='variables') %>% 
  mutate(model="Mod_1") %>% 
  bind_rows(
    as_tibble(Anova(m_SR_2), rownames ='variables')%>% 
      mutate(model="Mod_2")) %>% 
  bind_rows(
    as_tibble(Anova(m_SR_3b), rownames ='variables')%>% 
      mutate(model="Mod_3"))%>% 
  bind_rows(
    as_tibble(Anova(m_SR_4b), rownames ='variables')%>% 
      mutate(model="Mod_4")) %>% 
  bind_rows(
    as_tibble(Anova(m_SR_5b_1), rownames ='variables')%>% 
      mutate(model="Mod_5"))  %>% 
  bind_rows(
    as_tibble(Anova(m_SR_5b_2), rownames ='variables')%>% 
      mutate(model="Mod_6") %>% filter(variables=="Anthill_leveling")) %>% 
  mutate(predictor=recode_factor(variables,
                                habitat_corrected ="Habitat type",
                                Management_stability = "Management stability",
                                Grazing_legacy = "Grazing legacy",  
                                "poly(Mowing_frequency, 2)" = "Mowing frequency",
                                Mowing_delay = "Mowing delay",
                                Grazing_int_log = "Grazing intensity",
                                Grazer_type = "Grazer type",
                                Grazing_type = "Grazing type",
                                Grazing_season = "Grazing season",
                                Corralling = "Corralling",
                                Manuring_freq = "Manuring frequency",
                                Cow_dung_applied = "Cow-dung application",
                                Litter_removal = "Litter removal",
                                Anthill_leveling = "Anthill leveling",
                                Shrub_tree_removal = "Shrub/tree removal",
                                Moss_removal = "Moss removal",
                                Burning = "Burning",
                                Ploughing = "Ploughing",
                                humus_log = "Humus",
                                soil_NPK = "Soil NPK") ) %>% 
           mutate(variab_group=recode_factor(variables,
                                         soil_NPK = "Soil",  humus_log = "Soil",  
                                         Grazing_legacy = "Management type/stability/legacy",  
                                         habitat_corrected ="Management type/stability/legacy",
                                         Management_stability = "Management type/stability/legacy",
                                         Grazing_int_log = "Grazing",
                                         "poly(Mowing_frequency, 2)" = "Mowing",
                                         Grazer_type = "Grazing",
                                         Grazing_type = "Grazing",
                                         Grazing_season = "Grazing",
                                         Mowing_delay = "Mowing",
                                         Litter_removal = "Cleaning",
                                         #   Corralling = "Corralling",
                                         # Ploughing = "Ploughing",
                                         Shrub_tree_removal = "Cleaning",
                                         Anthill_leveling = "Cleaning",
                                         Moss_removal = "Cleaning",
                                         # Burning = "Burning",
                                         Cow_dung_applied = "Manuring/dung ",
                                         Manuring_freq = "Manuring/dung "
                                         )) %>% 
  arrange(predictor) %>% 
  dplyr::select(variab_group, predictor, Chisq,	Df,	"Pr(>Chisq)")

Anova_table

write.csv(Anova_table, "results/Tables_Anova.csv")

## > Model R2 ----
R2_model <- as_tibble(R2_m_SR_1b,rownames ='type') %>% mutate(model="Mod_1") %>% 
  bind_rows(
    as_tibble(R2_m_SR_2,rownames ='type') %>% mutate(model="Mod_2")) %>% 
  bind_rows(
    as_tibble(R2_m_SR_3b,rownames ='type') %>% mutate(model="Mod_3")) %>%
  bind_rows(
    as_tibble(R2_m_SR_4b,rownames ='type') %>% mutate(model="Mod_4")) %>%
  bind_rows(
    as_tibble(R2_m_SR_4b,rownames ='type') %>% mutate(model="Mod_5")) %>%
  filter(type=="delta") %>% dplyr::select(-type)

R2_model

write.csv(R2_model, "results/Tables_Model_R2.csv")


## > R2 partial ----

R2_partial <- as_tibble(R2_part_m_SR_1b) %>% mutate(model="Mod_1") %>% 
  bind_rows(
    as_tibble(R2_part_m_SR_2) %>% mutate(model="Mod_2")) %>% 
  bind_rows(
    as_tibble(R2_part_m_SR_3b) %>% mutate(model="Mod_3")) %>% 
  bind_rows(
    as_tibble(R2_part_m_SR_4b) %>% filter(!Effect=="Grazing_int_log") %>% 
                                    mutate(model="Mod_4")) %>% 
  bind_rows(
    as_tibble(R2_part_m_SR_5b_1) %>% mutate(model="Mod_5")) %>% 
  bind_rows(
    as_tibble(R2_part_m_SR_5b_2) %>% filter(Effect=="Anthill_levelingyes") %>% 
                                    mutate(model="Mod_6")) %>% 
  dplyr::select(Effect, model, Rsq, upper.CL, lower.CL) %>% 
  filter(!Effect=="Model") %>% # , !model=="Mod_4"
  mutate(predictor=recode_factor(Effect,
                                habitat_correctedpasture ="Management type",
                                Management_stability = "Management stability",
                                Grazing_legacy = "Grazing legacy", 
                                "poly(Mowing_frequency, 2)2" = "Mowing frequency",
                                "poly(Mowing_frequency, 2)1" = "Mowing frequency",
                                "Mowing_delayJuly-August" = "Mowing delay",
                                "Mowing_delayJune" = "Mowing delay",
                                Grazing_int_log = "Grazing intensity",
                                  "Grazer_typesheep, goat" = "Grazer type",
                                "Grazer_typemixed" = "Grazer type",
                                "Grazing_typefree" = "Grazing type",
                                "Grazing_typeherding" = "Grazing type",
                                "Grazing_seasonWhole season" = "Grazing season",
                                "Grazing_seasonSummer" = "Grazing season",
                                 Corrallingyes = "Corralling",
                                 Manuring_freq = "Manuring frequency",
                                Cow_dung_appliedpresent = "Cow manure application",
                                 Litter_removalyes = "Litter removal",
                                 Anthill_levelingyes = "Ant/mole-hill leveling",
                                 Shrub_tree_removalyes = "Shrub/tree removal",
                                 Moss_removalyes = "Moss removal",
                                 Burningyes = "Burning",
                                 Ploughing = "Ploughing",
                                 soil_NPK = "Soil NPK",  
                                 humus_log = "Humus"))%>% 
  arrange(predictor) %>% 
  group_by(predictor) %>% 
  summarise(Rsq=sum(Rsq)) %>% 
  ungroup() %>% 
  mutate(variab_group=recode_factor(predictor,
                                       "Management type" ="Management type/stability/legacy",
                                    "Management stability" = "Management type/stability/legacy",
                                    "Grazing legacy" = "Management type/stability/legacy",  
                                     "Mowing frequency" = "Mowing",
                                     "Mowing delay" = "Mowing",
                                    "Grazing intensity" = "Grazing",
                                    "Grazer type" = "Grazing",
                                    "Grazing type" = "Grazing",
                                    "Grazing season" = "Grazing",
                                    Corralling = "Corralling",
                                    "Manuring frequency" = "Manuring/dung ",
                                    "Cow manure application" = "Manuring/dung ",
                                    "Litter removal" = "Cleaning",
                                    "Ant/mole-hill leveling" = "Cleaning",
                                    "Shrub/tree removal" = "Cleaning",
                                    "Moss removal" = "Cleaning",
                                     Burning = "Burning",
                                     Ploughing = "Ploughing",
                                    "Soil NPK" = "Soil properties",
                                    "Humus" = "Soil properties")) %>%  
  arrange(predictor) 
         
         
R2_partial

write.csv(R2_partial, "results/Tables_R2_partial.csv")


# End ---------------


