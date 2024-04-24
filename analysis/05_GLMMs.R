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
Dat <- read_csv("data/Panoara_Dat.csv") %>%
  mutate(Grazer_type=str_replace_all(Grazer_type, "_", ", ")) %>% 
  mutate(Mowing_delay=fct_relevel(Mowing_delay,c("no mowing","June","July-August"))) %>%
  arrange(Mowing_delay) %>% 
  mutate(Dung_cover_log=log(Dung_cover+1)) %>% 
  mutate(SR_D_E_exper=case_when(is.na(SR_D_E_exper) ~ 0, .default=SR_D_E_exper),
 Abund_D_E_exper=case_when(is.na(Abund_D_E_exper) ~ 0, .default=Abund_D_E_exper)) %>% 
   filter(!Parcel_name=="Brade_1") # extreme outlier

str(Dat)
names(Dat)


 
# (1) Dung Experiment Data-----

## Species Richness ----

m1_SR_exp <- glmer (SR_D_E_exper ~   Abund_D_E_exper + (1|Farm), 
                    family = "poisson", 
                    data = Dat) 
ranef(m1_SR_exp) # random effects are zeros
summary(m1_SR_exp)
# Remove random effects

m2_SR_exp <- glm (SR_D_E_exper ~   Abund_D_E_exper, family = "poisson", 
                    data = Dat)
anova(m1_SR_exp, m2_SR_exp)


Anova(m2_SR_exp)
summary(m2_SR_exp)

plot_model(m2_SR_exp,type = "pred", terms = c("Abund_D_E_exper"), show.data=T, dot.alpha=0.3, title="")

SR_exp_pred <- get_model_data(m2_SR_exp,type = "pred", terms="Abund_D_E_exper[0:50, by=0.1]")

ggplot(SR_exp_pred, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1)+
  geom_point(data=Dat,
             aes(Abund_D_E_exper, SR_D_E_exper), fill = "#64ABCE",
             size=3, alpha=0.7, pch=21)+
  labs(y="Species richness", x='Abundance')+
  geom_line(linetype=1, linewidth=1) 


#----#           
m3_SR_exp <- glm (SR_D_E_exper ~  # Abund_D_E_exper + 
                    log(Grazing_intensity_A+1) + Grazer_diversity, # Grazer_type + 
                      family = "poisson", 
                    data = Dat) # strong outlier in SR_D_E_exper and Abund_D_E_exper

check_collinearity(m3_SR_exp)

Anova(m3_SR_exp)
summary(m3_SR_exp)

# Grazing_intensity
plot_model(m3_SR_exp,type = "pred", terms = c("Grazing_intensity_A"), show.data=T, dot.alpha=0.3, title="")
SR_exp_pred2 <- get_model_data(m3_SR_exp,type = "pred", terms="Grazing_intensity_A[0:1000, by=0.1]")

ggplot(SR_exp_pred2, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1)+
  geom_point(data=Dat%>% 
               mutate(Grazer_type=str_replace_all(Grazer_type, "_", ", ")),
             aes(Grazing_intensity_A, SR_D_E_exper, fill = Grazer_type),
             size=3, alpha=0.7, pch=21)+
  scale_fill_viridis(discrete=TRUE, option = "D")+
  labs(y="Species richness", x='Grazing intensity')+
  geom_line(linetype=1, linewidth=1)



# Grazing_intensity
plot_model(m3_SR_exp,type = "pred", terms = c("Grazer_diversity"), show.data=T, dot.alpha=0.3, title="")
SR_exp_pred3 <- get_model_data(m3_SR_exp,type = "pred", terms="Grazer_diversity[1:4, by=0.1]")


#

ggplot(SR_exp_pred3, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1)+
  geom_point(data=Dat %>% mutate(Grazer_type=str_replace_all(Grazer_type, "_", ", ")),
             aes(Grazer_diversity, SR_D_E_exper, fill = Grazer_type),
             size=3, alpha=0.7, pch=21)+
  scale_fill_viridis(discrete=TRUE, option = "D")+
  labs(y="Species richness", x='Grazer diversity')+
  geom_line(linetype=1, linewidth=1)


# Grazer type

m4_SR_exp <- glm (SR_D_E_exper ~  Grazer_type,
                  family = "poisson", 
                  data = Dat%>%
                    mutate(Grazer_type=str_replace_all(Grazer_type, "_", ", ")))
Anova(m4_SR_exp)
summary(m4_SR_exp)


emmeans_m4_SR_exp <- cld(emmeans(m4_SR_exp, list(pairwise ~ Grazer_type)), 
                            Letters = letters) %>% arrange(Grazer_type)
emmeans_m4_SR_exp

ggplot(Dat %>%  filter(!Parcel_name=="Brade_1", !is.na(Abund_D_E_exper)) %>%
         mutate(Grazer_type=str_replace_all(Grazer_type, "_", ", ")), 
       aes(Grazer_type, SR_D_E_exper)) + 
  geom_boxplot(outlier.shape = NA, notch = F)+
  geom_point(aes(fill = Grazing_intensity_A,size=Grazing_intensity_A), alpha=0.7, pch=21,
             position=position_jitter(width = 0.05, height = 0)) +
  scale_fill_viridis(option = "C")+
  theme_minimal()+
  labs(x ="Grazer type used for grazing in the parcel", y="Species richess", fill="Grazing intensity", size="Grazing intensity")+
  theme(axis.text = element_text(size=10), axis.title=element_text(size=12)) +
  
  geom_text(data=emmeans_m4_SR_exp,aes(x=Grazer_type, y=c(9, 10, 9, 5, 8),
                                          label=emmeans_m4_SR_exp$.group),vjust=0.5, hjust=0, 
            size=4, col="black" , position=position_dodge(0))


## Abundance ----

m1_Abund_exp <- glmer (Abund_D_E_exper ~  log(Grazing_intensity_A+1) + Grazer_diversity +  # Grazer_type + 
                       (1|Farm), family = "poisson", 
                  data = Dat %>% filter(!Parcel_name=="Brade_1")) # strong outlier in SR_D_E_exper and Abund_D_E_exper

check_convergence(m1_Abund_exp)
ranef(m1_Abund_exp) # random effects are not zeros
hist(ranef(m1_Abund_exp)$`Farm`[,1])

check_collinearity(m1_Abund_exp)


m2_Abund_exp <- glm (Abund_D_E_exper ~  log(Grazing_intensity_A+1) + Grazer_diversity,  # Grazer_type + 
                        family = "poisson", 
                       data = Dat %>% filter(!Parcel_name=="Brade_1")) # strong outlier in SR_D_E_exper and Abund_D_E_exper
anova(m1_Abund_exp, m2_Abund_exp)



Anova(m1_Abund_exp)
summary(m1_Abund_exp)

# Grazing_intensity
plot_model(m1_Abund_exp,type = "pred", terms = c("Grazing_intensity_A"), show.data=T, dot.alpha=0.3, title="")
Abund_exp_pred1 <- get_model_data(m1_Abund_exp,type = "pred", terms="Grazing_intensity_A[0:1000, by=0.1]")

ggplot(Abund_exp_pred1, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1)+
  geom_point(data=Dat %>% filter(!Parcel_name=="Brade_1")%>% 
               mutate(Grazer_type=str_replace_all(Grazer_type, "_", ", ")),
             aes(Grazing_intensity_A, Abund_D_E_exper, fill = Grazer_type),
             size=3, alpha=0.7, pch=21)+
  scale_fill_viridis(discrete=TRUE, option = "D")+
  labs(y="Abundance", x='Grazing intensity', fill="Grazer type")+
  geom_line(linetype=1, linewidth=1)

# Grazer_diversity
plot_model(m1_Abund_exp,type = "pred", terms = c("Grazer_diversity"), show.data=T, dot.alpha=0.3, title="")
Abund_exp_pred2 <- get_model_data(m1_Abund_exp,type = "pred", terms="Grazer_diversity[1:4, by=0.1]")

ggplot(Abund_exp_pred2, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1)+
  geom_point(data=Dat %>% filter(!Parcel_name=="Brade_1")%>% 
               mutate(Grazer_type=str_replace_all(Grazer_type, "_", ", ")),
             aes(Grazer_diversity, Abund_D_E_exper, fill = Grazer_type),
             size=3, alpha=0.7, pch=21)+
  scale_fill_viridis(discrete=TRUE, option = "D")+
  labs(y="Abundance", x='Grazer diversity', fill="Grazer type")+
  geom_line(linetype=1, linewidth=1)

# Grazer_type
m3_Abund_exp <- glmer (Abund_D_E_exper ~  Grazer_type + 
                         (1|Farm), family = "poisson", 
                       data = Dat %>% filter(!Parcel_name=="Brade_1")) # strong outlier in SR_D_E_exper and Abund_D_E_exper
ranef(m3_Abund_exp) # random effects are zeros

m4_Abund_exp <- glm (Abund_D_E_exper ~  Grazer_type,
                       family = "poisson", 
                       data = Dat %>% filter(!Parcel_name=="Brade_1") %>% # strong outlier in SR_D_E_exper and Abund_D_E_exper
                       mutate(Grazer_type=str_replace_all(Grazer_type, "_", ", "))) # correct names in Grazer_type

Anova(m4_Abund_exp)
summary(m4_Abund_exp)


emmeans_m4_Abund_exp <- cld(emmeans(m4_Abund_exp, list(pairwise ~ Grazer_type)), 
                      Letters = letters) %>% arrange(Grazer_type)
emmeans_m4_Abund_exp

ggplot(Dat %>%  filter(!Parcel_name=="Brade_1", !is.na(Abund_D_E_exper)) %>%
         mutate(Grazer_type=str_replace_all(Grazer_type, "_", ", ")), 
       aes(Grazer_type, Abund_D_E_exper)) + 
  geom_boxplot(outlier.shape = NA, notch = F)+
  geom_point(aes(fill = Grazing_intensity_A,size=Grazing_intensity_A), alpha=0.7, pch=21,
             position=position_jitter(width = 0.05, height = 0)) +
  scale_fill_viridis(option = "C")+
  theme_minimal()+
  labs(x ="Grazer type used for grazing in the parcel", y="Abundance", fill="Grazing intensity", size="Grazing intensity")+
  theme(axis.text = element_text(size=10), axis.title=element_text(size=12)) +
  
  geom_text(data=emmeans_m4_Abund_exp,aes(x=Grazer_type, y=c(23, 44, 23, 15, 26),
                                    label=emmeans_m4_Abund_exp$.group),vjust=0.5, hjust=0, 
            size=4, col="black" , position=position_dodge(0))

# (2) Field Data -----


## 1) Dung cover ----

m1_dung <- lmer (log(Dung_cover+1) ~  Manuring_freq +  # Cow_dung_applied +
                   (1|Farm), 
                 data = Dat)

plot(m1_dung) 
qqnorm(resid(m1_dung))
qqline(resid(m1_dung))

Anova(m1_dung)
summary(m1_dung)

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

m2_Excrem <- lm(log(Excrements_cover+1) ~                      
                    Grazing_intensity_A +  # Grazer_diversity +  
                    Grazer_type, 
                  data = Dat )

par(mfrow=c(2,2))
plot(m2_Excrem) 


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

m2_Excrem <- lm(log(Excrements_cover+1) ~   Grazer_type, 
                  data = Dat )
Anova(m2_Excrem)


## 3) Species Richness ----

m1_SR_field <- glmer(Plant_SR_vascular ~   
                       Abund_D_E_exper +   SR_D_E_exper +
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

# mowing

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





