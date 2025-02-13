# Purpose:  Perform Structural Equation Modelling for species composition (NMDS as main response)



# load packages ----------------------------------------------------------------
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


# data -------------------------------------------------------------------------

## soil PC data

Soil_PC <- read_csv("data/soil_NPK_PCA.csv")  %>% 
  mutate(soil_NPK=-1*soil_NPK_PC)  # reverse the NMDS scores to make it positively correlated with the nutrients

## NMDS data

NMDS_dat <- read_csv("data/NMDS_data.csv")
NMDS_dat

## Environmental data

Dat <- read_csv("data/Divers_LandUse_Soil_Variables.csv") %>%
  left_join(NMDS_dat, by="Parcel_name") %>% 
  left_join(Soil_PC, by="Parcel_name") %>% 
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


# Data wrangling
SEM.dat <- Dat %>% filter(!Parcel_name=="Farm_F_1") %>%  # extreme outlyer
  filter(!Dung_for_experiment==0) %>% # remove cases when dung was not found on the field, thus no experiment was done
  mutate(Grazing_int_log = log1p(Grazing_intensity_A)) %>% 
  mutate(Mowing_delay=case_when(Mowing_delay=="no mowing" ~ 0,
                                Mowing_delay=="July-August" ~ 1,
                                Mowing_delay=="June" ~ 2)) %>% 
  #  filter(!Mowing_frequency==3) %>% 
  mutate(Mowing_frq_scal=scale(Mowing_frequency, center=T, scale=T)) %>% 
  mutate(Mowing_frq_sqrd=Mowing_frq_scal^2)




SEM.dat
summary(SEM.dat)
names(SEM.dat)

# analysis----------------------------------------------------------------------

# (1) Individual models for the SEM ----

## mod 1: NMDS 1 (field) ----

m1_NMDS1_field <- lmer(NMDS1_field ~   
                         NMDS1_exper +                                                                   
                         NMDS2_exper +
                         Grazing_int_log  +  
                         Mowing_frequency +
                         Mowing_frq_sqrd +
                         Manuring_freq + 
                         humus_log + 
                         soil_NPK +
                         (1|Farm), data = SEM.dat) 


#ranef(m1_NMDS1_field) # random effects are not zeros
hist(ranef(m1_SR_field)$`Farm`[,1])

plot(m1_NMDS1_field)
qqnorm(resid(m1_NMDS1_field))
qqline(resid(m1_NMDS1_field))

check_convergence(m1_NMDS1_field)
check_collinearity(m1_NMDS1_field)

Anova(m1_NMDS1_field)
summary(m1_NMDS1_field)



### Extract R2 ----
R2_m1_NMDS1_field=MuMIn::r.squaredGLMM(m1_NMDS1_field)
R2_m1_NMDS1_field

# partial R2 
R2_part_m1_NMDS1_field=r2glmm::r2beta(m1_NMDS1_field,  partial = T, method = 'sgv')
R2_part_m1_NMDS1_field
plot(R2_part_m1_NMDS1_field %>% filter(!Effect=="Model"))

R2_NMDS1_partial <- as_tibble(R2_part_m1_NMDS1_field) %>% mutate(responce="NMDS1")


## mod 2: NMDS 2 (field) ----

m1_NMDS2_field <- lmer(NMDS2_field ~   
                         NMDS1_exper +                                                                   
                         NMDS2_exper +
                         Grazing_int_log  +  
                         Mowing_frequency +
                         Mowing_frq_sqrd +
                         Manuring_freq + 
                         humus_log + 
                         soil_NPK +
                         (1|Farm), data = SEM.dat) 


# ranef(m1_NMDS2_field) # random effects are not zeros
hist(ranef(m1_NMDS2_field)$`Farm`[,1])

plot(m1_NMDS2_field)
qqnorm(resid(m1_NMDS2_field))
qqline(resid(m1_NMDS2_field))

check_convergence(m1_NMDS2_field)
check_collinearity(m1_NMDS2_field)

Anova(m1_NMDS2_field)
summary(m1_NMDS2_field)


### Extract R2 ----
R2_m1_NMDS2_field=MuMIn::r.squaredGLMM(m1_NMDS2_field)
R2_m1_NMDS2_field

# partial R2 
R2_part_m1_NMDS2_field=r2glmm::r2beta(m1_NMDS2_field,  partial = T, method = 'sgv')
R2_part_m1_NMDS2_field
plot(R2_part_m1_NMDS2_field %>% filter(!Effect=="Model"))

R2_NMDS2_partial <- as_tibble(R2_part_m1_NMDS2_field) %>% mutate(responce="NMDS2")


R2_NMDS1_partial

write_csv(R2_NMDS1_partial%>% 
            bind_rows(R2_NMDS2_partial),
          "results/SEM_NMDS_R2_partial.csv")

## mod 3: NMDS 1 (experiment data) ----

m1_NMDS1_exp <- lmer (NMDS1_exper ~   
                        Grazing_int_log  +  
                        Manuring_freq + 
                    #  humus_log +
                     #  soil_NPK + 
                        (1|Farm), data = SEM.dat) 


ranef(m1_NMDS1_exp) # random effects are not zeros
hist(ranef(m1_NMDS1_exp)$`Farm`[,1])

m2_NMDS1_exp <- lm (NMDS1_exper ~   
                      Grazing_int_log  +  
                      Manuring_freq ,
                   #   humus_log +
                   #   soil_NPK,
                    data = SEM.dat)

check_collinearity(m2_NMDS1_exp)

Anova(m2_NMDS1_exp)

par(mfrow = c(2, 2))
plot(m2_NMDS1_exp)
par(mfrow = c(1, 1))


## mod 4: NMDS 2 (experiment data) ----

m1_NMDS2_exp <- lmer (NMDS2_exper ~ 
                        Grazing_int_log  +  
                        Manuring_freq + 
                    #    humus_log +
                     #  soil_NPK +
                        (1|Farm), data = SEM.dat) 

ranef(m1_NMDS2_exp) # random effects are not zeros
hist(ranef(m1_NMDS2_exp)$`Farm`[,1])


m2_NMDS2_exp <- lm(NMDS2_exper ~   
                     Grazing_int_log  +  
                     Manuring_freq , 
                  #   humus_log +
                   #  soil_NPK,
                   data = SEM.dat)

check_collinearity(m2_NMDS2_exp)

Anova(m2_NMDS2_exp)

par(mfrow = c(2, 2))
plot(m2_NMDS2_exp)
par(mfrow = c(1, 1))

## mod 5: Soil humus ----
m1_humus <- lmer(humus_log ~   
                  Grazing_int_log  +  
                  Mowing_frequency +
                   Manuring_freq +  
                   (1|Farm), data = SEM.dat) 

check_collinearity(m1_humus)

ranef(m1_humus) # random effects are zeros
hist(ranef(m1_humus)$`Farm`[,1])

# Remove random effects
m2_humus <- lm(humus_log ~    
                 Grazing_int_log  +  
                 Mowing_frequency +
                 Manuring_freq, data = SEM.dat)

Anova(m2_humus)

par(mfrow = c(2, 2))
plot(m2_humus)
par(mfrow = c(1, 1))

summary(m2_humus)


## mod 6:  Soil PC ----

m1_Soil_PC2 <- lmer(soil_NPK ~   
                      Grazing_int_log  +  
                      Mowing_frequency +
                      Manuring_freq +  
                      (1|Farm), data = SEM.dat) 

check_collinearity(m1_Soil_PC2)


plot(m1_Soil_PC2)
qqnorm(resid(m1_Soil_PC2))
qqline(resid(m1_Soil_PC2))


ranef(m1_Soil_PC2) # random effects 
hist(ranef(m1_Soil_PC2)$`Farm`[,1])

m2_Soil_PC2 <- lm(soil_NPK ~   
                    Grazing_int_log  +  
                    Mowing_frequency +
                    Manuring_freq, data = SEM.dat) 



Anova(m2_Soil_PC2)



# (2) Structural Equation modellind, PSEM ----

psem_model <- psem (m1_NMDS1_field, m1_NMDS2_field,
                    m2_NMDS1_exp, m2_NMDS2_exp,
                    m2_humus,
                    m2_Soil_PC2,
                    Mowing_frequency %~~% Mowing_frq_sqrd,
                    humus_log %~~%  soil_NPK)


# summary(psem_model, .progressBar =F,  conserve = TRUE)


coefic <- coefs(psem_model) 
coefic


#Fisher C statistic as global goodness of model fit:
fisherC(psem_model, .progressBar =F,  conserve = TRUE)

write_csv(coefic, "results/SEM_NMDS_coefs.csv")


# Rough draft of the plot (dirty plot)

plot(psem_model, digits=2, layout = "dot", 
     node_attrs = list(shape = "rectangle", 
                       color = "black",fillcolor = "white", width=1.4, 
                       distortion=5)) 


# (3) Indirect effects ------

Coefs <- coefic[1:26,]  %>%   
  dplyr::select(Response, Predictor, Std.Estimate) %>% 
  mutate(coefs=c("a1.1", "a1.2", "a1.3", "a1.4", "a1.5", "a1.6", "a1.7", "a1.8", 
                 "a2.1", "a2.2", "a2.3", "a2.4", "a2.5", "a2.6", "a2.7", "a2.8", 
                 "b1", "b2", # "b3", "b4",
                 "c1", "c2", # "c3", "c4",
                 "d1", "d2", "d3", 
                 "e1", "e2", "e3")) %>% 
  relocate(coefs)


Coefs

# Assign multiple values: 
library(zeallot)

Coefs$Std.Estimate
# assign letters to the coefficients
c(a1.1,a1.2,a1.3,a1.4,a1.5,a1.6,a1.7,a1.8, 
  a2.1,a2.2,a2.3,a2.4,a2.5,a2.6,a2.7,a2.8, 
  b1,b2, # b3,b4,
  c1,c2,# c3,c4,
  d1,d2,d3, 
  e1,e2,e3) %<-%   
  abs(Coefs$Std.Estimate)

a1
b1



## Manuring ----
# direct effect
manur_Direct_NMDS1 <- a1.6
manur_Direct_NMDS2 <- a2.6

# indirect
manur_Indir1_NMDS1 <- b2 * a1.1 # through NMDS1
manur_Indir2_NMDS1 <- c2 * a1.2 # through NMDS2
manur_Indir3_NMDS1 <- d3 * a1.7 # through humus
manur_Indir4_NMDS1 <- e3 * a1.8 # through soil nutrients

manur_Indir_seed_NMDS1 <- manur_Indir1_NMDS1 + manur_Indir2_NMDS1 
manur_Indir_soil_NMDS1 <- manur_Indir3_NMDS1 + manur_Indir4_NMDS1

manur_Indir1_NMDS2 <- b2 * a2.1 # through NMDS1
manur_Indir2_NMDS2 <- c2 * a2.2 # through NMDS2
manur_Indir3_NMDS2 <- d3 * a2.7 # through humus
manur_Indir4_NMDS2 <- e3 * a2.8 # through soil nutrients

manur_Indir_seed_NMDS2 <- manur_Indir1_NMDS2 + manur_Indir2_NMDS2 
manur_Indir_soil_NMDS2 <- manur_Indir3_NMDS2 + manur_Indir4_NMDS2

## Grazing ----
# direct effect
grazing_Direct_NMDS1 <- a1.3
grazing_Direct_NMDS2 <- a2.3

# indirect

grazing_Indir1_NMDS1 <- b1 * a1.1 # through NMDS 1
grazing_Indir2_NMDS1 <- c1 * a1.2 # through NMDS 2
grazing_Indir3_NMDS1 <- d1 * a1.7 # through humus
grazing_Indir4_NMDS1 <- e1 * a1.8 # through soil nutrients

grazing_Indir_seed_NMDS1 <- grazing_Indir1_NMDS1 + grazing_Indir2_NMDS1 
grazing_Indir_soil_NMDS1 <- grazing_Indir3_NMDS1 + grazing_Indir4_NMDS1

grazing_Indir1_NMDS2 <- b1 * a2.1 # through NMDS 1
grazing_Indir2_NMDS2 <- c1 * a2.2 # through NMDS 2
grazing_Indir3_NMDS2 <- d1 * a2.7 # through humus
grazing_Indir4_NMDS2 <- e1 * a2.8 # through soil nutrients

grazing_Indir_seed_NMDS2 <- grazing_Indir1_NMDS2 + grazing_Indir2_NMDS2 
grazing_Indir_soil_NMDS2 <- grazing_Indir3_NMDS2 + grazing_Indir4_NMDS2

## Mowing ----
# direct effect
mowing_Direct_NMDS1 <- a1.4+a1.5
mowing_Direct_NMDS2 <- a2.4+a2.5


# indirect

mowing_Indir3_NMDS1 <- d2 * a1.7 # through humus
mowing_Indir4_NMDS1 <- e2 * a1.8 # through soil nutrients

mowing_Indir_soil_NMDS1 <- mowing_Indir3_NMDS1 + mowing_Indir4_NMDS1

mowing_Indir3_NMDS2 <- d2 * a2.7 # through humus
mowing_Indir4_NMDS2 <- e2 * a2.8 # through soil nutrients

mowing_Indir_soil_NMDS2 <- mowing_Indir3_NMDS2 + mowing_Indir4_NMDS2


## seed dispersal ----
# direct effect
seed_Direct_NMDS1 <- a1.1+a1.2
seed_Direct_NMDS2 <- a2.1+a2.2
# indirect
# NA

Direct = c(grazing_Direct_NMDS1, grazing_Direct_NMDS2,
           manur_Direct_NMDS1, manur_Direct_NMDS2, 
           mowing_Direct_NMDS1, mowing_Direct_NMDS2,
           seed_Direct_NMDS1, seed_Direct_NMDS2)


Indirect_soil = c(grazing_Indir_soil_NMDS1, grazing_Indir_soil_NMDS2, 
                  manur_Indir_soil_NMDS1, manur_Indir_soil_NMDS2, 
                  mowing_Indir_soil_NMDS1,  mowing_Indir_soil_NMDS2,
                  0, 0)

Indirect_seed = c(grazing_Indir_seed_NMDS1, grazing_Indir_seed_NMDS2,
                  manur_Indir_seed_NMDS1, manur_Indir_seed_NMDS2,
                  0, 0,
                  0, 0)

Variable = c("Grazing", "Grazing", 
              "Manuring","Manuring",
              "Mowing",  "Mowing",
             "Dispersal", "Dispersal")

NMDS = c("NMDS1", "NMDS2", 
         "NMDS1", "NMDS2",
         "NMDS1", "NMDS2",
         "NMDS1", "NMDS2")

Coefs_summar <- tibble(Variable, NMDS, Direct, Indirect_soil, Indirect_seed) %>% 
  mutate(Total=Direct + Indirect_soil + Indirect_seed)

Coefs_summar

Coefs_summar_sum <- Coefs_summar %>% 
  dplyr::select(-NMDS) %>% 
  group_by(Variable) %>% 
 summarise_all(sum, na.rm=T) %>% 
  ungroup() %>% 
  mutate(Variable = factor(Variable, levels=c("Grazing", "Manuring", "Mowing",
                                            "Dispersal"))) %>% 
  mutate(Variable = fct_recode(Variable, "Grazing intensity" = "Grazing", 
                               "Manuring frequency" = "Manuring", 
                               "Mowing frequency" ="Mowing",
                               "Seed dispersal" = "Dispersal")) %>% 
  arrange(Variable)%>% 
  pivot_longer(-Variable, names_to = "Effect_type", values_to = "Std.Estimate") %>% 
  mutate(Std.Est=case_when (Std.Estimate==0 ~ NA,
                            .default=Std.Estimate)) %>%
  mutate(Effect_type=fct_recode(Effect_type,
                                "Indirect through soil properties" ="Indirect_soil",
                                "Indirect through seed dispersal" ="Indirect_seed",
                                "Direct effect" = "Direct",
                                "Total effect" = "Total")) %>% 
  mutate(Effect_type=fct_relevel(Effect_type, 
                              c("Indirect through soil properties", 
                              "Indirect through seed dispersal" ,
                              "Direct effect" ,
                              "Total effect")))


Coefs_summar_sum$Effect_type

write_csv(Coefs_summar_sum, "results/SEM_NMDS_coefs_Indirect.csv")

# Plots ----
Coefs_summar_sum  %>% 
  mutate(Effect_type=fct_relevel(Effect_type, 
                                 c("Indirect through soil properties", 
                                   "Indirect through seed dispersal" ,
                                   "Direct effect" ,
                                   "Total effect")))



dodge_width <- 0.5



plot <- ggplot(Coefs_summar_sum, aes(y =Effect_type , x = Std.Est, 
                                     color = Effect_type)) +
  geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
  geom_point(position = position_dodge(width = dodge_width), size = 4) +
  geom_errorbarh(aes(xmin = 0, xmax = Std.Estimate),
                 position = position_dodge(width = dodge_width), height = 0.1
  ) +
 facet_wrap(~Variable) +
  scale_color_manual(values = c("#BC7B3A", "#A021FF","#7CAE00","black"))+
#  MetBrewer::scale_color_met_d("Kandinsky") +
  theme_bw()+
  labs(y="Effect type", color="Effect type",
       x="Standardised effect size (absolute value)", title="Effects on community composition")

plot

plot2 <- ggplot(Coefs_summar_sum%>% 
                  filter(!Variable=="Mowing frequency"),
                aes(y =Effect_type , x = Std.Est, 
                                  color = Effect_type, fill = Effect_type)) +
  geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
  geom_col(position = position_dodge(width = 0.8), width = 0.8,size=0.5) +
  #  geom_point(position = position_dodge(width = dodge_width), size = 4) +
  # geom_errorbarh(aes(xmin = 0, xmax = Std.Estimate),
  #                position = position_dodge(width = dodge_width), height = 0.1
  # ) +
  facet_wrap(~Variable #, scales = "free_x"
             ) +
  #  MetBrewer::scale_color_met_d("Kandinsky") +
  scale_color_manual(values = c("#BC7B3A", "#A021FF","#719E00","gray22"))+
  scale_fill_manual(values = c("#DFBA95", "#E6C5FF","#C4E0B2","gray77"))+
  theme_bw()+
  theme(legend.position = "None",
        axis.text.x = element_text(size=7, color="black"), 
        axis.text.y = element_text(size=9, color="black", face="bold"),
        axis.title=element_text(size=10, face="bold")) +
  labs(y=" ", color=" ", fill=" ",
       x="Standardised effect size (absolute value)" )#, title="Effects on community composition")

plot2

