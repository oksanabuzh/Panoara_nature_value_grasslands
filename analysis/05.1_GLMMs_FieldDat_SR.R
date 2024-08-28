# Purpose: GLMMs for Species Richness (SR)
rm(list = ls())

# load packages

library(tidyverse)
library(ggplot2)
# library(devtools)
# devtools::install_github("strengejacke/sjPlot")
# library(insight)
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
          axis.textsize.x = 1, axis.textsize.y = 1, axis.textcolor = "black",
          axis.title.color = "black", axis.title.size = 1.4,
          geom.linetype = 1) #legend.pos = "None", 



# Read data  ----

Soil_PC <- read_csv("data/soil_PC.csv")

Soil_PC2 <- read_csv("data/soil_PC2.csv")

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
  left_join(Soil_PC, by="Parcel_name")%>% 
  left_join(Soil_PC2, by="Parcel_name")


str(Dat)
names(Dat)


Dat <- Dat %>% 
  mutate(zz = Plant_SR_vascular) # for predicting from the glmmPQL models




# Data Explorations ----
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
ggplot(Dat, aes(Shrub_tree_removal, Plant_SR_vascular))+geom_boxplot() + geom_point()
ggplot(Dat, aes(log(Last_ploughing), Plant_SR_vascular))+geom_point()
ggplot(Dat, aes(Corralling, Plant_SR_vascular))+geom_boxplot()
ggplot(Dat, aes(soil_CN, Plant_SR_vascular))+ geom_point()
ggplot(Dat, aes(soil_pH, Plant_SR_vascular))+ geom_point()
ggplot(Dat, aes(log(humus), Plant_SR_vascular))+ geom_point()
ggplot(Dat, aes(soil_P_acces, Plant_SR_vascular))+ geom_point()
ggplot(Dat, aes(soil_Ca_acces, Plant_SR_vascular))+ geom_point()
ggplot(Dat, aes(soil_MG_acces, Plant_SR_vascular))+ geom_point()
ggplot(Dat, aes(log(soil_K_acces), Plant_SR_vascular))+ geom_point()


# GLMM models -----
## Driver effects: -----

## habitat, mang. stability ----
m_SR_1 <- glmer(Plant_SR_vascular ~  habitat_corrected + 
                  Management_stability +
                (1|Farm), family = "poisson", data = Dat) 

check_convergence(m_SR_1)
check_collinearity(m_SR_1)
summary(m_SR_1)
Anova(m_SR_1)

check_overdispersion(m_SR_1) # Overdispersion detected.

# Change family to negative binomial
m_SR_1b <- glmer.nb(Plant_SR_vascular ~ 
            habitat_corrected +
             Management_stability +
           (1|Farm), data = Dat) 

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

#### habitat ----
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



 ### Management_stability ----
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
  theme(axis.text = element_text(size=10), 
        axis.title=element_text(size=12)) +
  theme(axis.title.x=element_text(vjust=-0.7), 
        axis.title.y=element_text(vjust=1.5),
        legend.title = element_text(face="plain")) 



# Grazing, Mowing----

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
plot_model(m_SR_2, type = "pred", terms = c("Mowing_frequency"), show.data=T, dot.alpha=0.3, title="") 


# check for hump-shape for grazing:
m_SR_2_ <- glmer(Plant_SR_vascular ~  
                   #  Grazing_int_log +
                   poly(Grazing_int_log , 2)  +
                   poly(Mowing_frequency, 2)  +
                  (1|Farm), family = "poisson", data = Dat) 

anova(m_SR_2, m_SR_2_) # selects against the hump-shape for grazing
Anova(m_SR_2)

check_overdispersion(m_SR_2) # Slight overdispersion detected


## -> Mowing_frequency -----
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
  geom_line(linetype=1, linewidth=1) 


# R2 for the entire model 
# R2m and R2c are marginal (for fixed predictors) and 
## conditional (for fixed and random predictors) coefficients of determination
MuMIn::r.squaredGLMM(m_SR_2)

# partial R2
R2=r2glmm::r2beta(m_SR_2,  partial = T, method = 'sgv', data=Dat)
R2
plot(R2)



#  Grazer_type & Mowing_delay ----

m_SR_3 <- glmer(Plant_SR_vascular ~  
                #   Grazing_int_log +
                  Grazer_type +
                 #  Grazing_type +
                 # Grazing_season +
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
# R2 for the entire model 
# R2m and R2c are marginal (for fixed predictors) and 
## conditional (for fixed and random predictors) coefficients of determination
MuMIn::r.squaredGLMM(m_SR_3b)

# partial R2
r2glmm::r2beta(m_SR_3b,  partial = T, method = 'sgv', data=Dat)

## -> Grazer_type ----

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


## -> Mowing delay ----

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
                  Grazing_int_log +
                 # Grazer_type +
                  Grazing_type +
                  Grazing_season +
                 # poly(Mowing_frequency, 2)  + # 
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
                     Grazing_int_log +
                    #   Grazer_type +
                       Grazing_type +
                       Grazing_season +
                    #  Mowing_delay +
                    #   poly(Mowing_frequency, 2)  +
                      (1|Farm),
                    data = Dat) 

Anova(m_SR_4b)

check_convergence(m_SR_4b)
check_collinearity(m_SR_4b)


emmeansm_SR_4b_Grazing_season <- cld(emmeans(m_SR_4b, list(pairwise ~ Grazing_season)), 
                                    Letters = letters, test="Tukey") %>% arrange(Grazing_season)
emmeansm_SR_4b_Grazing_season


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

emmeansm_SR_4b_Grazing_type <- cld(emmeans(m_SR_4b, list(pairwise ~ Grazing_type)), 
                                     Letters = letters, test="Tukey") %>% arrange(Grazing_type)
emmeansm_SR_4b_Grazing_type

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


Dat%>% pull(Moss_removal)

m_SR_5a <- glm(Plant_SR_vascular ~  
                    # Management_stability +
                    # habitat +
                    Grazing_int_log +
                  poly(Mowing_frequency, 2)  +
                 #  Mowing_frequency +
                 # Last_ploughing + # 
                  Ploughing + 
                  # Shrub_tree_removal + # Cleaning +  # Shrub_tree_removal + # correlates with Cleaning
                   Litter_removal +
                 #  Moss_removal + 
                    Anthill_leveling + #  Molehill_leveling +# data are the same as Anthill_leveling
                  # Burning + 
                    Corralling +
                 Manuring_freq +
                Cow_dung_applied +
                 humus_log + # 
             # PC1_soil +  
               PC1_soil_2 , 
                  family = "poisson", data = Dat)  


Anova(m_SR_5a)
check_collinearity(m_SR_5a)
check_overdispersion(m_SR_5a)


m_SR_5b <- glmer(Plant_SR_vascular ~  
                   # Management_stability +
                   # habitat +
                   Grazing_int_log +
                   poly(Mowing_frequency, 2)  +
                   #  Mowing_frequency +
                  #  poly(Ploughing) + # Crops_planted + # correlates with Last_ploughing
                    Last_ploughing +
                   # Shrub_tree_removal + # Cleaning +  # Shrub_tree_removal + # correlates with Cleaning
                   Litter_removal +
                   #  Moss_removal + 
                   Anthill_leveling + #  Molehill_leveling +# data are the same as Anthill_leveling
                   # Burning + 
                   Corralling +
                   Manuring_freq +
                   #  Cow_dung_applied  +
                   PC1_soil+ PC2_soil +
                  (1|Farm), family = "poisson", data = Dat)  

Anova(m_SR_5b)


anova(m_SR_5b, m_SR_5a)


ranef(m_SR_5b) # random effects are not zeros
hist(ranef(m_SR_5b)$`Farm`[,1])

check_convergence(m_SR_5b) 
check_collinearity(m_SR_5b)
check_convergence(m_SR_5b)
check_overdispersion(m_SR_5b)  # very low overdispersion


plot_model(m_SR_5b, type = "pred", terms = c("Last_ploughing"), show.data=T, dot.alpha=0.3, title="") 
plot_model(m_SR_5a, type = "pred", terms = c("Litter_removal"), show.data=T, dot.alpha=0.3, title="") 
plot_model(m_SR_5a, type = "pred", terms = c("Corralling"), show.data=T, dot.alpha=0.3, title="") 
plot_model(m_SR_5a, type = "pred", terms = c("Manuring_freq"), show.data=T, dot.alpha=0.3, title="") 
plot_model(m_SR_5a, type = "pred", terms = c("Cow_dung_applied"), show.data=T, dot.alpha=0.3, title="") 
plot_model(m_SR_5a, type = "pred", terms = c("Management_stability"), show.data=T, dot.alpha=0.3, title="") 
plot_model(m_SR_5a, type = "pred", terms = c("PC1_soil"), show.data=T, dot.alpha=0.3, title="") 

plot_model(m_SR_5a, type = "pred", terms = c("humus"), show.data=T, dot.alpha=0.3, title="") 
plot_model(m_SR_5a, type = "pred", terms = c("soil_P_acces"), show.data=T, dot.alpha=0.3, title="") 

plot_model(m_SR_5a, type = "pred", terms = c("soil_pH"), show.data=T, dot.alpha=0.3, title="") 
plot_model(m_SR_5a, type = "pred", terms = c("soil_CN"), show.data=T, dot.alpha=0.3, title="") 

### Last_ploughing ----
m_SR_5a_plough_pred <- get_model_data(m_SR_5a,type = "pred", terms="Ploughing[0:0.2, by=0.01]")

min(Dat$Ploughing)

ggplot(m_SR_5a_plough_pred, aes(x, predicted)) +
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



m_SR_5a_plough_pred <- get_model_data(m_SR_5a,type = "pred", terms="Last_ploughing[0:60, by=0.01]")

ggplot(m_SR_5a_plough_pred, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1)+
  geom_point(data=Dat, aes(Last_ploughing, Plant_SR_vascular), 
             size=3, alpha=0.5, pch=21, fill ="#64ABCE",
             position=position_jitter(width=0.08, height=0))+
  # scale_fill_manual(values=c("blue", "red", "green"))+
  theme(axis.title.x=element_text(vjust=-0.1), axis.title.y=element_text(vjust=2)) +
  labs(y="Species richness", x='Last ploughing (yrs. ago)')+
  scale_x_continuous(breaks=seq(0, 60, 10), 
                     labels=paste0(c("0", "10","20", "30", "40", "50", ">50"))) +
  geom_line(linetype=1, linewidth=1) 


### Litter_removal ----


emmeans_m_SR_5a_Litt_remov <- cld(emmeans(m_SR_5a, list(pairwise ~ Litter_removal)), 
                                    Letters = letters, test="Tukey") %>% arrange(Litter_removal)
emmeans_m_SR_5a_Litt_remov

ggplot(Dat, aes(Litter_removal, Plant_SR_vascular)) + 
  geom_boxplot(outlier.shape = NA, notch = F)+
  geom_point(alpha=0.6, pch=21, size=3,
             position=position_jitter(width = 0.05, height = 0),
             fill ="#64ABCE") +
 # scale_fill_manual(values=c("blue", "red", "green"))+
  theme_bw()+
  labs(x ="Litter removal", y="Species richness")+
  theme(axis.text = element_text(size=10), axis.title=element_text(size=15)) +
  theme(axis.title.x=element_text(vjust=-0.1), axis.title.y=element_text(vjust=2)) +
  
  geom_text(data=emmeans_m_SR_5a_Litt_remov,aes(x=Litter_removal, y=c(86, 90),
                                                  label=emmeans_m_SR_5a_Litt_remov$.group),vjust=0.5, hjust=0, 
            size=4, col="black" , position=position_dodge(0))




### Corralling ----




emmeans_m_SR_5a_Corralling <- cld(emmeans(m_SR_5a, list(pairwise ~ Corralling)), 
                                  Letters = letters, test="Tukey") %>% arrange(Corralling)
emmeans_m_SR_5a_Corralling

ggplot(Dat, aes(Corralling, Plant_SR_vascular)) + 
  geom_boxplot(outlier.shape = NA, notch = F)+
  geom_point(alpha=0.6, pch=21, size=3,
             position=position_jitter(width = 0.05, height = 0),
             fill ="#64ABCE") +
  # scale_fill_manual(values=c("blue", "red", "green"))+
  theme_bw()+
  labs(x ="Corraling", y="Species richness")+
  theme(axis.text = element_text(size=10), axis.title=element_text(size=15)) +
  theme(axis.title.x=element_text(vjust=-0.1), axis.title.y=element_text(vjust=2)) +
  
  geom_text(data=emmeans_m_SR_5a_Corralling,aes(x=Corralling, y=c(90, 70),
                                                label=c("a", "b")),vjust=0.5, hjust=0, 
            size=4, col="black" , position=position_dodge(0))


### Anthill_leveling  #  Molehill_leveling ----

emmeans_m_SR_5a_Anthill_removal <- cld(emmeans(m_SR_5a, list(pairwise ~ Anthill_leveling)), 
                                       Letters = letters, test="Tukey") %>% arrange(Anthill_leveling)
emmeans_m_SR_5a_Anthill_removal

ggplot(Dat, aes(Anthill_leveling, Plant_SR_vascular)) + 
  geom_boxplot(outlier.shape = NA, notch = F)+
  geom_point(alpha=0.6, pch=21, size=3, fill ="#64ABCE",
             position=position_jitter(width = 0.05, height = 0))  +
  # scale_fill_manual(values=c("blue", "red", "green"))+
  theme_bw()+
  labs(x ="Anthill / molehill leveling", y="Species richness")+
  theme(axis.text = element_text(size=10), axis.title=element_text(size=15)) +
  theme(axis.title.x=element_text(vjust=-0.1), axis.title.y=element_text(vjust=2)) +
  geom_text(data=emmeans_m_SR_5a_Anthill_removal,aes(x=Anthill_leveling, y=c(75, 90),
                                                     label=emmeans_m_SR_5a_Anthill_removal$.group),
            vjust=0.5, hjust=0, 
            size=4, col="black" , position=position_dodge(0)) 


###Shrub_tree_removal ----


m_SR_5c <- glm(Plant_SR_vascular ~  
                 Grazing_int_log +
                 poly(Mowing_frequency, 2)  +
                 poly(Last_ploughing, 2) + # Crops_planted + # correlates with Last_ploughing
                 Shrub_tree_removal + # Cleaning +  # Shrub_tree_removal  correlates with Anthill_leveling
                 Litter_removal +
                 #  Moss_removal + 
                 #  Anthill_leveling + #  Molehill_leveling +# data are the same as Anthill_leveling
                 # Burning + 
                 Corralling +
                Manuring_freq +
               #  Cow_dung_applied  +
                 PC1_soil+ PC2_soil, 
               family = "poisson", data = Dat)  


Anova(m_SR_5c)
summary(m_SR_5c)

emmeans_m_SR_5c_shrub_removal <- cld(emmeans(m_SR_5c, list(pairwise ~ Shrub_tree_removal)), 
                                     Letters = letters, test="Tukey") %>% arrange(Shrub_tree_removal)
emmeans_m_SR_5c_shrub_removal

ggplot(Dat, aes(Shrub_tree_removal, Plant_SR_vascular)) + 
  geom_boxplot(outlier.shape = NA, notch = F)+
  geom_point(alpha=0.6, pch=21, size=3,
             position=position_jitter(width = 0.05, height = 0),
             fill ="#64ABCE") +
  # scale_fill_manual(values=c("blue", "red", "green"))+
  theme_bw()+
  labs(x ="Shrub / tree removal", y="Species richness")+
  theme(axis.text = element_text(size=10), axis.title=element_text(size=15)) +
  theme(axis.title.x=element_text(vjust=-0.1), axis.title.y=element_text(vjust=2))+
  geom_text(data=emmeans_m_SR_5c_shrub_removal,aes(x=Shrub_tree_removal, y=c(70, 90),
                                                   label=emmeans_m_SR_5c_shrub_removal$.group),
            vjust=0.5, hjust=0, 
            size=4, col="black" , position=position_dodge(0)) 




### Moss_removal ----

emmeans_m_SR_5a_Moss_removal <- cld(emmeans(m_SR_5a, list(pairwise ~ Moss_removal)), 
                                Letters = letters, test="Tukey") %>% arrange(Moss_removal)
emmeans_m_SR_5a_Moss_removal

ggplot(Dat, aes(Moss_removal, Plant_SR_vascular)) + 
  geom_boxplot(outlier.shape = NA, notch = F)+
  geom_point(alpha=0.6, pch=21, size=3,
             position=position_jitter(width = 0.05, height = 0),
             fill ="#64ABCE") +
  # scale_fill_manual(values=c("blue", "red", "green"))+
  theme_bw()+
  labs(x ="Moss removal", y="Species richness")+
  theme(axis.text = element_text(size=10), axis.title=element_text(size=15)) +
  theme(axis.title.x=element_text(vjust=-0.1), axis.title.y=element_text(vjust=2)) 




### Burning ----

emmeans_m_SR_5a_Burning <- cld(emmeans(m_SR_5a, list(pairwise ~ Burning)), 
                                    Letters = letters, test="Tukey") %>% arrange(Burning)

emmeans_m_SR_5a_Burning

ggplot(Dat, aes(Moss_removal, Plant_SR_vascular)) + 
  geom_boxplot(outlier.shape = NA, notch = F)+
  geom_point(alpha=0.6, pch=21, size=3,
             position=position_jitter(width = 0.05, height = 0),
             fill ="#64ABCE") +
  # scale_fill_manual(values=c("blue", "red", "green"))+
  theme_bw()+
  labs(x ="Burning", y="Species richness")+
  theme(axis.text = element_text(size=10), axis.title=element_text(size=15)) +
  theme(axis.title.x=element_text(vjust=-0.1), axis.title.y=element_text(vjust=2))





### Manuring_freq   ----
m_SR_5a_Manuring_freq_pred <- get_model_data(m_SR_5a,type = "pred", terms="Manuring_freq[0:1, by=0.0001]")

ggplot(m_SR_5a_Manuring_freq_pred, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1)+
  geom_point(data=Dat, aes(Manuring_freq, Plant_SR_vascular), 
             size=3, alpha=0.5, pch=21, fill ="#64ABCE",
             position=position_jitter(width=0.01, height=0))+
  # scale_fill_manual(values=c("blue", "red", "green"))+
  theme(axis.title.x=element_text(vjust=-0.1), axis.title.y=element_text(vjust=2)) +
  labs(y="Species richness", x='Manuring frequency')+
  
  geom_line(linetype=1, linewidth=1) 




### Cow_dung_applied ----


emmeans_m_SR_5a_Cow_dung  <- cld(emmeans(m_SR_5a, list(pairwise ~ Cow_dung_applied)), 
                               Letters = letters, test="Tukey") %>% arrange(Cow_dung_applied)
emmeans_m_SR_5a_Cow_dung

ggplot(Dat, aes(Cow_dung_applied, Plant_SR_vascular)) + 
  geom_boxplot(outlier.shape = NA, notch = F)+
  geom_point(alpha=0.6, pch=21, size=3,
             position=position_jitter(width = 0.05, height = 0),
             fill ="#64ABCE") +
  # scale_fill_manual(values=c("blue", "red", "green"))+
  theme_bw()+
  labs(x ="Cow dung applied", y="Species richness")+
  theme(axis.text = element_text(size=10), axis.title=element_text(size=15)) +
  theme(axis.title.x=element_text(vjust=-0.1), axis.title.y=element_text(vjust=2)) 


### Soil PC -----

m_SR_2_soilPC1_pred <- get_model_data(m_SR_5_glm,type = "pred", 
                                      terms="PC1_soil[-3:8, by=0.001]")

ggplot(m_SR_2_soilPC1_pred, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1)+
  geom_point(data=Dat,
             aes(PC1_soil, Plant_SR_vascular), 
             size=3, alpha=0.7, pch=21, fill ="#64ABCE")+
  # scale_fill_viridis(discrete=TRUE, option = "D")+
  labs(y="Species richness", x='Soil PC1', fill="Grazer type")+
  geom_line(linetype=1, linewidth=1) 




# Grazing Legacy effect ------

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


Anova(m_SR_6b)
check_collinearity(m_SR_6b)


plot_model(m_SR_6b, type = "pred", terms = c("Grazing_legacy"), show.data=T, dot.alpha=0.3, title="") 
plot_model(m_SR_6b, type = "pred", terms = c("Mowing_frequency"), show.data=T, dot.alpha=0.3, title="") 


m_SR_6b_legac_pred <- get_model_data(m_SR_6b, type = "pred", 
                                   terms="Grazing_legacy[0:4, by=0.001]")

ggplot(m_SR_6b_legac_pred, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1)+
  geom_point(data=Dat,
             aes(Grazing_legacy, Plant_SR_vascular, fill=Grazer_type), 
             size=3, alpha=0.7, pch=21)+
  scale_fill_viridis(discrete=TRUE, option = "D")+
  labs(y="Species richness", x='Cumulative grazing intensity', fill= "Grazer type")+
  geom_line(linetype=1, linewidth=1) 

# End ---------------


