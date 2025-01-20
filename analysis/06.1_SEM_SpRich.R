# Purpose:  Perform Structural Equation Modelling for Species richness (as main response)
#dev.off()
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

Diver_NMDS_dat <- read_csv("data/Diversity_&_NMDS_data.csv")


  
# Biodiversity, NMDS and envoronmental data
Data <- read_csv("data/LandUse_soil_variables.csv") %>%
  left_join(Soil_PC, by="Parcel_name") %>% 
  left_join(Diver_NMDS_dat, by="Parcel_name") %>% 
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
SEM.dat <- Data %>% filter(!Parcel_name=="Farm_F_1") %>%  # extrime outlyer
  filter(!Dung_for_experiment==0) %>% 
  mutate(Grazing_int_log = log1p(Grazing_intensity_A)) %>% 
  mutate(SR_Exper=case_when(is.na(SR_Exper) ~ 0, .default=SR_Exper),
        # Abund_D_E_exper=case_when(is.na(Abund_D_E_exper) ~ 0, .default=Abund_D_E_exper),
         abundance_Exper=case_when(is.na(abundance_Exper) ~ 0, .default=abundance_Exper)) %>% 
  mutate(Mowing_delay=case_when(Mowing_delay=="no mowing" ~ 0,
                                Mowing_delay=="July-August" ~ 1,
                                Mowing_delay=="June" ~ 2)) %>% 
  #  filter(!Mowing_frequency==3) %>% 
  mutate(Mowing_frq_scal=scale(Mowing_frequency, center=T, scale=T)) %>% 
  mutate(Mowing_frq_sqrd=Mowing_frq_scal^2,
         SR_Exper_log = log1p(SR_Exper))


Data$Plant_SR_vascular
Data$SR_D_E_exper
Data$SR_Exper
Data$No_dung_experiment
SEM.dat$SR_VP_field

SEM.dat
summary(SEM.dat)
names(SEM.dat)

# analysis----

# (1) Individual models for the SEM ----

## mod 1: Species Richness (field data)----

m1_SR_field <- glmer(SR_VP_field ~ #Plant_SR_vascular ~   
                       abundance_Exper +  
                       SR_Exper_log + 
                       Grazing_int_log  +   
                       Mowing_frequency +
                       Mowing_frq_sqrd +
                       Manuring_freq + 
                       humus_log + 
                       soil_NPK +
                       (1|Farm),  family = "poisson", 
                     data = SEM.dat) 



check_convergence(m1_SR_field)
check_collinearity(m1_SR_field)
check_overdispersion(m1_SR_field) # low overdispersion


Anova(m1_SR_field)
summary(m1_SR_field)

# random effects
ranef(m1_SR_field) 
hist(ranef(m1_SR_field)$`Farm`[,1])

### Extract R2 ----
R2_m1_SR_field=MuMIn::r.squaredGLMM(m1_SR_field)
R2_m1_SR_field

# partial R2 
R2_part_m1_SR_field=r2glmm::r2beta(m1_SR_field,  partial = T, method = 'sgv', data=SEM.dat)
R2_part_m1_SR_field
plot(R2_part_m1_SR_field %>% filter(!Effect=="Model"))

R2_partial <- as_tibble(R2_part_m1_SR_field) %>% mutate(responce="SR")

write_csv(R2_partial, "results/SEM_SR_R2_partial.csv")

## mod 2: Species Richness (experiment data) ----

m1_SR_exp <- lmer (SR_Exper_log ~   
                     #  abundance_Exper +
                     Grazing_int_log  +  
                     Manuring_freq +                      
                  #   soil_NPK +
                  #   humus_log +
                     (1|Farm), data = SEM.dat) 



Anova(m1_SR_exp)


plot(m1_SR_exp)
qqnorm(resid(m1_SR_exp))
qqline(resid(m1_SR_exp))

check_convergence(m1_SR_exp)
check_collinearity(m1_SR_exp)

ranef(m1_SR_exp)
hist(ranef(m1_SR_exp)$`Farm`[,1])


plot_model(m1_SR_exp, type = "pred", terms = c("Grazing_int_log"), 
           show.data=T, dot.alpha=0.3, title="") 

## mod 3: Abundance (Experiment data) ----
m1_Abund_exp <- glmer (abundance_Exper ~ 
                         Grazing_int_log  + 
                         Manuring_freq + 
                       #  soil_NPK +
                      #   humus_log + 
                         (1|Farm), family = "poisson", 
                       data = SEM.dat) 

check_convergence(m1_Abund_exp)
check_collinearity(m1_Abund_exp)

Anova(m1_Abund_exp)
summary(m1_Abund_exp)

ranef(m1_Abund_exp) 
hist(ranef(m1_Abund_exp)$`Farm`[,1])

check_overdispersion(m1_Abund_exp)

# high overdispersion, use negative binomial
m2_Abund_exp <- glmer.nb(abundance_Exper ~
                           Grazing_int_log  +   
                           Manuring_freq + 
                      #   humus_log +
                      #     soil_NPK +
                           (1|Farm),
                         data=SEM.dat) 


Anova(m2_Abund_exp)

check_convergence(m2_Abund_exp)
check_collinearity(m2_Abund_exp)
check_overdispersion(m2_Abund_exp)

summary(m2_Abund_exp)


## mod 4: Soil humus ----
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


## mod 5:  Soil PC ----

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

psem_model <- psem (m1_SR_field, 
                    m1_SR_exp, 
                    m2_Abund_exp,
                    m1_humus,
                    m1_Soil_PC2,
                    Mowing_frequency %~~% Mowing_frq_sqrd,
                    abundance_Exper %~~% SR_Exper_log,
                    humus_log %~~%  soil_NPK)




# summary(psem_model, .progressBar =F,  conserve = TRUE)

coefic <- coefs(psem_model) 
coefic

#Fisher C statistic as global goodness of model fit:
fisherC(psem_model, .progressBar =F,  conserve = TRUE)

write_csv(coefic, "results/SEM_SR_coefs.csv")





plot(psem_model)

plot(psem_model, digits=2, layout = "dot", 
     node_attrs = list(shape = "rectangle", 
                       color = "black",fillcolor = "white", width=1.4, 
                       distortion=5)) 


# (3) Indirect effects ------

Coefs <- coefic[1:18,]  %>%   
  dplyr::select(Response, Predictor, Std.Estimate) %>% 
  mutate(coefs=c("a1", "a2", "a3", "a4", "a5", "a6", "a7", "a8", 
                  "b2", "b3", # "b4", "b5",
                 "c1", "c2", # "c3", 
                 "d1", "d2", "d3", 
                 "e1", "e2", "e3")) %>% 
  relocate(coefs)


Coefs

# Assign multiple values: 
library(zeallot)

Coefs$Std.Estimate
# assign letters to the coefficients
c(a1,a2,a3,a4,a5,a6,a7,a8, 
  b2,b3, # b4,b5,
  c1,c2,# c3,
  d1,d2,d3, 
  e1,e2,e3) %<-%   
  Coefs$Std.Estimate

a1




## Manuring ----
# direct effect
manur_Direct <- a6
# indirect
manur_Indir1 <- c2 * a1 # through SeedAbund
manur_Indir2 <- b3 * a2 # through SR_Exper
manur_Indir3 <- d3 * a7 # through humus
manur_Indir4 <- e3 * a8 # through soil nutrients

manur_Indir_seed <- manur_Indir1 + manur_Indir2 
manur_Indir_soil <- manur_Indir3 + manur_Indir4

## Grazing ----
# direct effect
grazing_Direct <- a3
# indirect
grazing_Indir1 <- c1 * a1 # through SeedAbund
grazing_Indir2 <- b2 * a2 # through SR_Exper
grazing_Indir3 <- d1 * a7 # through humus
grazing_Indir4 <- e1 * a8 # through soil nutrients

grazing_Indir_seed <- grazing_Indir1 + grazing_Indir2 
grazing_Indir_soil <- grazing_Indir3 + grazing_Indir4



## Mowing ----
# direct effect
mowing_Direct <- a4
# indirect

mowing_Indir3 <- d2 * a7 # through humus
mowing_Indir4 <- e2 * a8 # through soil nutrients

mowing_Indir_soil <- mowing_Indir3 + mowing_Indir4

## seed dispersal ----
# direct effect
seed_Direct <- a1#+a2
# indirect

Direct = c(grazing_Direct, manur_Direct, mowing_Direct,  seed_Direct)
Indirect_soil = c(grazing_Indir_soil, manur_Indir_soil, mowing_Indir_soil,  0)
Indirect_seed = c(grazing_Indir_seed,manur_Indir_seed, 0,  0)
Variable = c("Grazing", "Manuring", "Mowing",  "Dispersal")

Coefs_summar <- tibble(Variable, Direct, Indirect_soil, Indirect_seed) %>% 
  mutate(Total=Direct + Indirect_soil + Indirect_seed) %>% 
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
                                "Total effcet" = "Total")) %>% 
  mutate(Effect_type=fct_relevel(Effect_type, 
                                 c("Indirect through soil properties", 
                                   "Indirect through seed dispersal" ,
                                   "Direct effect" ,
                                   "Total effcet")))


Coefs_summar

write_csv(Coefs_summar, "results/SEM_SR_coefs_Indirect.csv")

## > Plot effects----

Coefs_summar <- read_csv("results/SEM_SR_coefs_Indirect.csv") %>% 
  mutate(Effect_type=fct_relevel(Effect_type, 
                                 c("Indirect through soil properties", 
                                   "Indirect through seed dispersal" ,
                                   "Direct effect" ,
                                   "Total effect")))

dodge_width <- 0.5

plot <- ggplot(Coefs_summar, aes(y =Effect_type , x = Std.Est, 
                                     color = Effect_type)) +
  geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
  geom_point(position = position_dodge(width = dodge_width), size = 4) +
  geom_errorbarh(aes(xmin = 0, xmax = Std.Estimate),
                 position = position_dodge(width = dodge_width), height = 0.1
  ) +
  facet_wrap(~Variable) +
  #  MetBrewer::scale_color_met_d("Kandinsky") +
  scale_color_manual(values = c("#BC7B3A", "#C77CFF","#00BFC4","gray37"))+
  theme_bw()+
  labs(y="Effect type", color="Effect type",
       x="Standardised effect strength and direction", title="Effects on plant species richness")

plot



plot2 <- ggplot(Coefs_summar %>% 
                  filter(!Variable=="Mowing frequency"), 
                aes(y =Effect_type , x = Std.Est, 
                                 color = Effect_type, fill = Effect_type)) +
  geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
  geom_col(position = position_dodge(width = 0.8), width = 0.8,size=0.5) +
#  geom_point(position = position_dodge(width = dodge_width), size = 4) +
 # geom_errorbarh(aes(xmin = 0, xmax = Std.Estimate),
 #                position = position_dodge(width = dodge_width), height = 0.1
 # ) +
  facet_wrap(~Variable) + #, scales = "free_x") +
  #  MetBrewer::scale_color_met_d("Kandinsky") +
  scale_color_manual(values = c("#BC7B3A", "#A021FF","#719E00","gray22"))+
  scale_fill_manual(values = c("#DFBA95", "#E6C5FF","#C4E0B2","gray77"))+
  theme_bw()+
  theme(legend.position = "None",
        axis.text.x = element_text(size=7, color="black"), 
        axis.text.y = element_text(size=9, color="black", face="bold"),
        axis.title=element_text(size=10, face="bold")) +
  labs(y=" ", color="", fill="",
       x="Standardised effect size and direction") #, title="Effects on plant species richness")

plot2
