
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

# set thepe for plots
set_theme(base = theme_bw(),
          axis.textsize.x = 1, axis.textsize.y = 1, axis.textcolor = "black",
          axis.title.color = "black", axis.title.size = 1.4,
          legend.pos = "None", geom.linetype = 1)


# Read data  ----
Dat <- read_csv("data/Panoara_Dat.csv")

str(Dat)
names(Dat)


# Test random effects ----
m_a <- glmer (Plant_SR_vascular ~ Grazing_intensity_A + SR_D_E_exper +
              Abund_D_E_exper +
              (1|Farm), family = "poisson", data = Dat %>% 
              filter(!Parcel_name=="Brade 1")) # strong outlier in SR_D_E_exper and Abund_D_E_exper

check_convergence(m_a)

ranef(m_a) # random effects are not zeros
hist(ranef(m_a)$`Farm`[,1])

# remove random effects
m_b <- glm (Plant_SR_vascular ~ Grazing_intensity_A + SR_D_E_exper +
              Abund_D_E_exper ,
            # (1|Farm),  
            family = "poisson", data = Dat%>% 
              filter(!Parcel_name=="Brade 1"))

# likelihood-ratio test
anova(m_a, m_b) 
# including random effect significantly improves the model

# check R2 for random effects
# R2m and R2c are marginal (for fixed predictors) and 
## conditional (for fixed and random predictors) coefficients of determination
MuMIn::r.squaredGLMM(m_a)
MuMIn::r.squaredGLMM(m_a)[1,"R2c"] - MuMIn::r.squaredGLMM(m_a)[1,"R2m"]
# [1] 0.2055548
# R2c shows that random effect explain 20% of data  


# Check model  ----

m1 <- glmer (Plant_SR_vascular ~ Grazing_intensity_A + SR_D_E_exper +
                Abund_D_E_exper +
                (1|Farm), family = "poisson", data = Dat %>% 
                filter(!Parcel_name=="Brade 1")) # strong outlier in SR_D_E_exper and Abund_D_E_exper

check_convergence(m1)

# check model
plot(m1)
qqnorm(resid(m1))
qqline(resid(m1))

# check multicolinearity
check_collinearity(m1)

Anova(m1)
summary(m1)

# check plots
plot_model(m1,type = "pred", terms = c("SR_D_E_exper"), show.data=T, dot.alpha=0.3, title="")
plot_model(m1,type = "pred", terms = c("Grazing_intensity_A"), show.data=T, dot.alpha=0.3, title="")
plot_model(m1,type = "pred", terms = c("Abund_D_E_exper"), show.data=T, dot.alpha=0.3, title="")

# check overdispersion
sum(residuals(m1, type = "pearson")^2) / df.residual(m1)
# slightly overdispersed

# Try negative binomial GLMM -----
m2 <- glmer.nb (Plant_SR_vascular ~ Grazing_intensity_A + SR_D_E_exper +
                  Abund_D_E_exper + (1|Farm), 
                data = Dat %>% filter(!Parcel_name=="Brade 1")) 

check_convergence(m2)

check_collinearity(m2)

Anova(m2)
summary(m2)
drop1(m2, test = "Chisq")

# check overdispersion
sum(residuals(m2, type = "pearson")^2) / df.residual(m2)
# negative binomial GLMM corrected for the overdispersion

# check model
plot(m2)
qqnorm(resid(m2))
qqline(resid(m2))

# check plots
plot_model(m2,type = "pred", terms = c("SR_D_E_exper"), show.data=T, dot.alpha=0.3, title="")
plot_model(m2,type = "pred", terms = c("Grazing_intensity_A"), show.data=T, dot.alpha=0.3, title="")
plot_model(m2,type = "pred", terms = c("Abund_D_E_exper"), show.data=T, dot.alpha=0.3, title="")

# check R2 ----
# R2 for the entire model 
# R2m and R2c are marginal (for fixed predictors) and 
## conditional (for fixed and random predictors) coefficients of determination
MuMIn::r.squaredGLMM(m2)

# partial R2
r2glmm::r2beta(m2,  partial = T, method = 'sgv')

# check in piecewiseSEM
summary(psem(m2))
