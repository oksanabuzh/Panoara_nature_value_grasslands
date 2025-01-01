# Purpose:  check relationship with plant biomass and nutrients

# load packages ----
library(tidyverse)
library(ggplot2)
library(performance)
library(lme4)
library(lmerTest)
library(car)
library(emmeans)
library(multcomp)
library(piecewiseSEM)
library(sjPlot)


# data ----

## soil PC data

Soil_PC <- read_csv("data/soil_NPK_PCA.csv") %>% 
  mutate(soil_NPK=-1*soil_NPK_PC)  # reverse the NMDS scores to make it positively correlated with the nutrients

Plant_NPK_mass <- read_csv("data/plant_biomass_CN_PK.PCA.csv") %>% 
  mutate(plant_PK_PC=-1*plant_PK_PC) # reverse the NMDS scores to make it positively correlated with the nutrients


# Biodiversity, NMDS and envoronmental data
Data <- read_csv("data/Panoara_Dat.csv") %>%
  mutate(Grazer_type=str_replace_all(Grazer_type, "_", ", ")) %>% 
  mutate(Grazer_type=fct_relevel(Grazer_type,c("cow","sheep, goat","mixed"))) %>%
  mutate(Mowing_delay=fct_relevel(Mowing_delay,c("no mowing","June","July-August"))) %>%
  arrange(Mowing_delay) %>% 
  mutate(Dung_cover_log=log(Dung_cover+1)) %>% 
  mutate(Grazing_int_log = log1p(Grazing_intensity_A),
         Grazing_legacy =log1p(Grazing_intensity_B),
         Corralling=case_when(Corralling=="0" ~ "no", 
                              Corralling=="1" ~ "yes")) %>% 
  mutate(humus_log=log1p(humus)) %>% 
  mutate(Cow_dung=case_when(Cow_dung_applied=="present" ~ "1", 
                            Cow_dung_applied=="absent" ~ "0"))

# join with soil PC
Dat <- Data %>% 
  left_join(Soil_PC, by="Parcel_name") %>% 
  left_join(Plant_NPK_mass, by="Parcel_name") %>% 
  mutate(Grazing_int_log = log1p(Grazing_intensity_A)) %>% 
  mutate(Mowing_delay=case_when(Mowing_delay=="no mowing" ~ 0,
                                Mowing_delay=="July-August" ~ 1,
                                Mowing_delay=="June" ~ 2)) %>% 
  #  filter(!Mowing_frequency==3) %>% 
  mutate(Mowing_frq_scal=scale(Mowing_frequency, center=T, scale=T)) %>% 
  mutate(Mowing_frq_sqrd=Mowing_frq_scal^2,
         SR_Exper_log = log1p(SR_Exper))





m1 <- glmer(Plant_SR_vascular ~   
                    #   abundance_Exper +  
                    #   SR_Exper_log + 
                    #   Grazing_int_log  +   
                     #  Mowing_frequency +
                    #   Mowing_frq_sqrd +
                     #  Manuring_freq + 
                       humus_log + 
                       soil_NPK +
                       plant_biomass + plant_CN + plant_PK_PC +
                       (1|Farm),  family = "poisson", 
                     data = Dat) 



check_convergence(m1)
check_collinearity(m1)
check_overdispersion(m1) 


Anova(m1)
summary(m1)



m2 <- glmer.nb(Plant_SR_vascular ~   
                           #   abundance_Exper +  
                           #   SR_Exper_log + 
                           #   Grazing_int_log  +   
                           #  Mowing_frequency +
                           #   Mowing_frq_sqrd +
                           #  Manuring_freq + 
                           humus_log + 
                           soil_NPK +
                           plant_biomass + plant_CN + plant_PK_PC +
                           (1|Farm),
                         data=Dat) 


Anova(m2)

check_convergence(m2)
check_collinearity(m2)
check_overdispersion(m2)

summary(m2)


# plant_biomass -----

m3 <- lmer(plant_biomass ~   
               Grazing_int_log  +   
              Mowing_frequency +
              Manuring_freq + 
         #    Cow_dung_applied +
              humus_log + 
              soil_NPK +
             # plant_biomass + plant_CN + plant_PK_PC +
              (1|Farm),   
            data = Dat %>% 
             filter(!humus>25)) 



check_convergence(m3)
check_collinearity(m3)


Anova(m3)

plot_model(m3, type = "pred", terms = c("Manuring_freq"), 
           show.data=T, dot.alpha=0.3, title="") 





# plant_CN -----
m4 <- lmer(plant_CN~   
             Grazing_int_log  +   
             Mowing_frequency +
             Manuring_freq + 
         #    Cow_dung_applied +
             humus_log + 
             soil_NPK +
             # plant_biomass + plant_CN + plant_PK_PC +
             (1|Farm),   
           data = Dat %>% 
             filter(!humus>25)) 



check_convergence(m4)
check_collinearity(m4)


Anova(m4)

plot_model(m4, type = "pred", terms = c("Grazing_int_log"), 
           show.data=T, dot.alpha=0.3, title="") 


# plant_PK_PC -----
m5 <- lmer(plant_PK_PC ~   
             Grazing_int_log  +   
             Mowing_frequency +
             Manuring_freq + 
        #     Cow_dung_applied +
             humus_log + 
             soil_NPK +
             # plant_biomass + plant_CN + plant_PK_PC +
             (1|Farm),   
           data = Dat %>% filter(!humus>25)) 



check_convergence(m5)
check_collinearity(m5)


Anova(m5)

plot_model(m5, type = "pred", terms = c("humus_log"), 
           show.data=T, dot.alpha=0.3, title="") 
