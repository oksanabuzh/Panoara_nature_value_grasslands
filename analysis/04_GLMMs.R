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
library("viridis")          


# set thepe for plots
set_theme(base = theme_bw(),
          axis.textsize.x = 1, axis.textsize.y = 1, axis.textcolor = "black",
          axis.title.color = "black", axis.title.size = 1.4,
          geom.linetype = 1) #legend.pos = "None", 


# Read data  ----
Dat <- read_csv("data/Panoara_Dat.csv")

str(Dat)
names(Dat)


 
# Dung Experiment Data-----

## species richness
m1_SR_exp <- glmer (SR_D_E_exper ~  # Abund_D_E_exper + 
                     # D_E_exp_type +
                       D_E_exp_categ +
                      # D_E_exp_categ_2 +
                    #  log(Excrements_cover+1) + log(Dung_cover+1) +
                      # Manuring_freq + # Cow_dung_applied +
                      log(Grazing_intensity_A) + Grazer_diversity + # Grazer_type + 
               # Mowing_delay + Mowing_frequency + Mowing_method +
               (1|Farm), family = "poisson", data = Dat %>% 
               filter(!Parcel_name=="Brade 1")) # strong outlier in SR_D_E_exper and Abund_D_E_exper

check_collinearity(m1_SR_exp)

Anova(m1_SR_exp)
summary(m1_SR_exp)
plot_model(m1_SR_exp,type = "pred", terms = c("D_E_exp_categ"), show.data=T, dot.alpha=0.3, title="")
plot_model(m1_SR_exp,type = "pred", terms = c("D_E_exp_categ_2"), show.data=T, dot.alpha=0.3, title="")
plot_model(m1_SR_exp,type = "pred", terms = c("D_E_exp_type"), show.data=T, dot.alpha=0.3, title="")
plot_model(m1_SR_exp,type = "pred", terms = c("Excrements_cover"), show.data=T, dot.alpha=0.3, title="")
plot_model(m1_SR_exp,type = "pred", terms = c("Dung_cover"), show.data=T, dot.alpha=0.3, title="")
plot_model(m1_SR_exp,type = "pred", terms = c("Manuring_freq"), show.data=T, dot.alpha=0.3, title="")
plot_model(m1_SR_exp,type = "pred", terms = c("Cow_dung_applied"), show.data=T, dot.alpha=0.3, title="")
plot_model(m1_SR_exp,type = "pred", terms = c("Grazing_intensity_A"), show.data=T, dot.alpha=0.3, title="")
plot_model(m1_SR_exp,type = "pred", terms = c("Grazer_type"), show.data=T, dot.alpha=0.3, title="")
plot_model(m1_SR_exp,type = "pred", terms = c("Grazer_diversity"), show.data=T, dot.alpha=0.3, title="")


# abundance
m1_Abund_exp <- glmer (Abund_D_E_exper ~  # SR_D_E_exper+
                         # D_E_exp_type +
                         D_E_exp_categ_2 +
                        # Excrements_cover + Dung_cover +
                         # Manuring_freq +  # Cow_dung_applied + 
                         log(Grazing_intensity_A) + Grazer_diversity + # Grazer_type +
                      # Mowing_delay +  Mowing_frequency + Mowing_method +
                      (1|Farm), family = "poisson", data = Dat %>% 
                      filter(!Parcel_name=="Brade 1")) # strong outlier in SR_D_E_exper and Abund_D_E_exper

check_convergence(m1_Abund_exp)
Anova(m1_Abund_exp)
summary(m1_Abund_exp)

plot_model(m1_Abund_exp,type = "pred", terms = c("D_E_exp_categ"), show.data=T, dot.alpha=0.3, title="")

plot_model(m1_Abund_exp,type = "pred", terms = c("D_E_exp_type"), show.data=T, dot.alpha=0.3, title="")
plot_model(m1_Abund_exp,type = "pred", terms = c("Manuring_freq"), show.data=T, dot.alpha=0.3, title="")
plot_model(m1_Abund_exp,type = "pred", terms = c("Excrements_cover"), show.data=T, dot.alpha=0.3, title="")
plot_model(m1_Abund_exp,type = "pred", terms = c("Dung_cover"), show.data=T, dot.alpha=0.3, title="")
plot_model(m1_Abund_exp,type = "pred", terms = c("Cow_dung_applied"), show.data=T, dot.alpha=0.3, title="")

plot_model(m1_Abund_exp,type = "pred", terms = c("Grazing_intensity_A"), show.data=T, dot.alpha=0.3, title="")
plot_model(m1_Abund_exp,type = "pred", terms = c("Grazer_type"), show.data=T, dot.alpha=0.3, title="")
plot_model(m1_Abund_exp,type = "pred", terms = c("Grazer_diversity"), show.data=T, dot.alpha=0.3, title="")




# (2) Field Data -----


## 1) Dung cover ----

m1_dung <- lmer (log(Dung_cover+1) ~  Manuring_freq +  # Cow_dung_applied +
                   (1|Farm), 
                 data = Dat # %>%  filter(!Parcel_name=="Brade 1")
)

plot(m1_dung) 
qqnorm(resid(m1_dung))
qqline(resid(m1_dung))

Anova(m1_dung)

# check plots
plot_model(m1_dung,type = "pred", terms = c("Manuring_freq"), show.data=T, dot.alpha=0.3, title="")


m1_pred_Manuring_freq <- get_model_data(m1_dung,type = "pred", terms="Manuring_freq[0:1,by=0.001]")

ggplot(m1_pred_Manuring_freq, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1)+
  geom_point(data=Dat, aes(Manuring_freq, Dung_cover), 
             position = position_jitter(w = 0.05, h = 0), 
             size=3, alpha=0.7, pch=21, fill="#64ABCE")+
  labs(y="Dung cover", x='Manuring frequency')+
  geom_line(linetype=1, linewidth=1) 



## 2) Excrement cover ----

m1_Excrem <- lmer(log(Excrements_cover+1) ~                      
                    Grazing_intensity_A +  # Grazer_diversity +  
                    Grazer_type + 
                    (1|Farm), 
                  data = Dat )

plot(m1_Excrem) 
qqnorm(resid(m1_Excrem))
qqline(resid(m1_Excrem))

Anova(m1_Excrem)

# check plots
plot_model(m1_Excrem,type = "pred", terms = c("Grazing_intensity_A"), show.data=T, dot.alpha=0.3, title="")
plot_model(m1_Excrem,type = "pred", terms = c("Grazer_diversity"), show.data=T, dot.alpha=0.3, title="")
plot_model(m1_Excrem,type = "pred", terms = c("Grazer_type"), show.data=T, dot.alpha=0.3, title="")

# Grazing intensity
m1_Excrem_pred <- get_model_data(m1_Excrem,type = "pred", terms="Grazing_intensity_A[0:1000, by=0.1]")

ggplot(m1_Excrem_pred, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1)+
  geom_point(data=Dat %>% 
               mutate(Grazer_type=str_replace_all(Grazer_type, "_", ", ")),
             aes(Grazing_intensity_A, Excrements_cover, fill = Grazer_type), 
             size=3, alpha=0.7, pch=21)+
  scale_fill_viridis(discrete=TRUE, option = "D")+
  labs(y="Excrement cover", x='Grazing intensity', fill="Grazer type")+
  geom_line(linetype=1, linewidth=1) 


## 3) Species Richness ----
m_SR_field <- glmer.nb (Plant_SR_vascular ~ Grazing_intensity_A + 
                      # Experiment data
                          Abund_D_E_exper +  SR_D_E_exper +
                      # Management
                    Mowing_frequency + # Mowing_method +  Mowing_delay +
                   # Manuring_freq +
                    (1|Farm), 
                data = Dat %>% filter(!Parcel_name=="Brade 1")) 

check_convergence(m_SR_field)

check_collinearity(m_SR_field)

Anova(m_SR_field)
summary(m_SR_field)

# check overdispersion
sum(residuals(m_SR_field, type = "pearson")^2) / df.residual(m_SR_field)
# negative binomial GLMM corrected for the overdispersion

# check model
plot(m2)
qqnorm(resid(m2))
qqline(resid(m2))

# check plots
plot_model(m2,type = "pred", terms = c("Grazing_intensity_A"), show.data=T, dot.alpha=0.3, title="") +
plot_model(m2,type = "pred", terms = c("SR_D_E_exper"), show.data=T, dot.alpha=0.3, title="") +
plot_model(m2,type = "pred", terms = c("Abund_D_E_exper"), show.data=T, dot.alpha=0.3, title="")

### check R2 ----
# R2 for the entire model 
# R2m and R2c are marginal (for fixed predictors) and 
## conditional (for fixed and random predictors) coefficients of determination
MuMIn::r.squaredGLMM(m2)

# partial R2
r2glmm::r2beta(m2,  partial = T, method = 'sgv')

# check in piecewiseSEM
summary(psem(m2))


## Replace NAs with "0" in experiment data ----
m3 <- glmer.nb (Plant_SR_vascular ~ Grazing_intensity_A + SR_D_E_exper +
                  Abund_D_E_exper + (1|Farm), 
                data = Dat %>% filter(!Parcel_name=="Brade 1") %>% 
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
               filter(!Parcel_name=="Brade 1")) # strong outlier in SR_D_E_exper and Abund_D_E_exper

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
               filter(!Parcel_name=="Brade 1") %>% 
                        mutate(SR_D_E_exper=case_when(is.na(SR_D_E_exper) ~ 0, .default=SR_D_E_exper),
                               Abund_D_E_exper=case_when(is.na(Abund_D_E_exper) ~ 0, .default=Abund_D_E_exper))
               )

# check plots
plot_model(m4b,type = "pred", terms = c("Abund_D_E_exper"), show.data=T, dot.alpha=0.3, title="")
              
Anova(m4b)





