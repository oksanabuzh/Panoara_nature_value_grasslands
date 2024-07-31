# Purpose:  Structural equation modelling for Species richness

Soil_PC <- read_csv("data/soil_PC.csv") %>% 
  mutate(PC1_soil_log= log(PC1_soil +10),
         PC2_soil_log= log(PC2_soil +10))


Soil_PC2 <- read_csv("data/soil_PC2.csv")


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


# filter(!Parcel_name=="Brade_1") # extreme outlier


Dat <- Data %>% 
  left_join(Soil_PC, by="Parcel_name") %>% 
  left_join(Soil_PC2, by="Parcel_name")

Dat %>% 
  pull(Parcel_name, abundance_Exper)

SEM.dat <- Dat %>% filter(!Parcel_name=="Brade_1") %>%  # extrime outlyer
  mutate(SR_Exper=case_when(is.na(SR_Exper) ~ 0, .default=SR_Exper),
         Abund_D_E_exper=case_when(is.na(Abund_D_E_exper) ~ 0, .default=Abund_D_E_exper),
         abundance_Exper=case_when(is.na(abundance_Exper) ~ 0, .default=abundance_Exper),
         Evenness_Exper=case_when(is.na(Evenness_Exper) ~ 0, .default=Evenness_Exper),
         Shannon_Exper=case_when(is.na(Shannon_Exper) ~ 0, .default=Shannon_Exper)) %>% 
#  filter(!is.na(SR_D_E_exper)) %>% 
  mutate(Grazing_int_log = log1p(Grazing_intensity_A)) %>% 
 # mutate(Mowing_delay=case_when(Mowing_delay=="no mowing" ~ 0,
 #                               Mowing_delay=="July-August" ~ 1,
 #                               Mowing_delay=="June" ~ 2)) %>% 
filter(!Mowing_frequency==3) %>% 
  mutate(Mowing_frq_scal=scale(Mowing_frequency, center=T, scale=T)) %>% 
  mutate(Mowing_frq_sqrd=Mowing_frq_scal^2,
         SR_Exper_log = log1p(SR_Exper+1))




SEM.dat
summary(SEM.dat)
names(SEM.dat)

# Species Richness (field data) ----

m1_SR_field <- glmer(Plant_SR_vascular ~   
                       # Abund_D_E_exper + #
                       abundance_Exper +  
                       # Shannon_Exper + #  
                       # Evenness_Exper + # 
                       SR_Exper_log +
                       Grazing_int_log  +  #  Grazer_diversity +
                       Mowing_frequency +
                    #   Mowing_frq_sqrd +
                       Manuring_freq + 
                       humus_log + #  PC1_soil + 
                       PC1_soil_2 +
                       # humus_log  + PC1_soil_2 + #PC2_soil + #  +
                       (1|Farm),  family = "poisson", 
                     data = SEM.dat) 


check_convergence(m1_SR_field)
check_collinearity(m1_SR_field)

sum(residuals(m1_SR_field, type = "pearson")^2) / df.residual(m1_SR_field)
check_overdispersion(m1_SR_field)



Anova(m1_SR_field)
summary(m1_SR_field)

ranef(m1_SR_field) # random effects are not zeros
hist(ranef(m1_SR_field)$`Farm`[,1])

m2_SR_field <- glmer.nb(Plant_SR_vascular ~   
                          # Abund_D_E_exper + #
                          abundance_Exper +  
                          # Shannon_Exper + #  
                          # Evenness_Exper + # 
                          SR_Exper_log +
                          Grazing_int_log  +  #  Grazer_diversity +
                          Mowing_frequency +
                          Mowing_frq_sqrd +
                          Manuring_freq + 
                          humus_log   + PC1_soil_2 + #humus_log + PC2_soil + #  +
                          (1|Farm), data = SEM.dat) 


Anova(m2_SR_field)


ranef(m2_SR_field) # random effects are not zeros
hist(ranef(m2_SR_field)$`Farm`[,1])

check_convergence(m2_SR_field)
check_collinearity(m2_SR_field)


summary(m2_SR_field)
plot_model(m2_SR_field, type = "pred", terms = c("abundance_Exper"), show.data=T, dot.alpha=0.3, title="") 
plot_model(m2_SR_field, type = "pred", terms = c("SR_Exper"), show.data=T, dot.alpha=0.3, title="") 
plot_model(m2_SR_field, type = "pred", terms = c("Evenness_Exper"), show.data=T, dot.alpha=0.3, title="") 
plot_model(m2_SR_field, type = "pred", terms = c("Shannon_Exper"), show.data=T, dot.alpha=0.3, title="") 

plot_model(m2_SR_field, type = "pred", terms = c("Mowing_frequency"), show.data=T, dot.alpha=0.3, title="") 
plot_model(m2_SR_field, type = "pred", terms = c("Grazing_int_log"), show.data=T, dot.alpha=0.3, title="") 
plot_model(m2_SR_field, type = "pred", terms = c("Manuring_freq"), show.data=T, dot.alpha=0.3, title="") 
plot_model(m2_SR_field, type = "pred", terms = c("PC1_soil_log"), show.data=T, dot.alpha=0.3, title="") 

plot_model(m2_SR_field, type = "pred", terms = c("humus_log"), show.data=T, dot.alpha=0.3, title="") 
plot_model(m2_SR_field, type = "pred", terms = c("soil_P_acces"), show.data=T, dot.alpha=0.3, title="") 
plot_model(m2_SR_field, type = "pred", terms = c("soil_CN"), show.data=T, dot.alpha=0.3, title="") 


# Species Richness (experiment data) ----


m1_SR_exp <- lmer (SR_Exper_log ~ # Shannon_Exper ~   
                   #  abundance_Exper +
                     Grazing_int_log  +  # Grazer_diversity +
                     Mowing_frequency +
                     #  Mowing_frq_sqrd +
                     Manuring_freq + 
                     PC1_soil_2 +
                     (1|Farm), data = SEM.dat) 



Anova(m1_SR_exp)


plot(m1_SR_exp)
qqnorm(resid(m1_SR_exp))
qqline(resid(m1_SR_exp))

check_convergence(m1_SR_exp)

check_collinearity(m1_SR_exp)

ranef(m1_SR_exp) # random effects are not zeros
hist(ranef(m1_SR_exp)$`Farm`[,1])

m_SR_2_mowing_pred <- get_model_data(m1_SR_exp,type = "pred", terms="Mowing_frequency[0:2, by=0.1]")

ggplot(m_SR_2_mowing_pred, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1)+
  geom_point(data=SEM.dat, aes(Mowing_frequency, SR_Exper_log, fill = Mowing_delay), 
             size=3, alpha=0.5, pch=21,
             position=position_jitter(width=0.05, height=0))+
  scale_fill_manual(values=c("blue", "red", "green"))+
  theme(axis.title.x=element_text(vjust=-0.1), axis.title.y=element_text(vjust=2)) +
  labs(y="Species richness", x='Mowing frequency', fill="Mowing delay")+
  geom_line(linetype=1, linewidth=1) 




# Abundance (Experiment data) ----
m1_Abund_exp <- glmer (abundance_Exper ~ 
                         Grazing_int_log  + #  Grazer_diversity +
                         Mowing_frequency +
                         # Mowing_frq_sqrd +
                         Manuring_freq + 
                         (1|Farm), family = "poisson", 
                       data = SEM.dat) 

check_convergence(m1_Abund_exp)
check_collinearity(m1_Abund_exp)

Anova(m1_Abund_exp)
summary(m1_Abund_exp)

ranef(m1_Abund_exp) # random effects are not zeros
hist(ranef(m1_Abund_exp)$`Farm`[,1])

check_overdispersion(m1_Abund_exp)

m2_Abund_exp <- glmer.nb(abundance_Exper ~
                           Grazing_int_log  +  # Grazer_diversity +
                           Mowing_frequency +
                         # Mowing_frq_sqrd +
                           Manuring_freq + 
                           (1|Farm),
                         data=SEM.dat) 

check_convergence(m2_Abund_exp)
check_collinearity(m2_Abund_exp)
check_overdispersion(m2_Abund_exp)

Anova(m2_Abund_exp)
summary(m2_Abund_exp)


# Soil PC ----

m1_Soil_PC <- lmer(humus_log ~   
                     Grazing_int_log  +  # Grazer_diversity +
                     Mowing_frequency +
                 #    Mowing_frq_sqrd +
                     Manuring_freq +  
                     (1|Farm), data = SEM.dat) 

m1_Soil_PC <- lm(humus_log ~   
                   Grazing_int_log  +  # Grazer_diversity +
                   Mowing_frequency +
                 #  Mowing_frq_sqrd +
                   Manuring_freq, data = SEM.dat)

Anova(m1_Soil_PC)

plot(m1_Soil_PC)
qqnorm(resid(m1_Soil_PC))
qqline(resid(m1_Soil_PC))

check_convergence(m1_Soil_PC)

check_collinearity(m1_Soil_PC)

ranef(m1_Soil_PC) # random effects are not zeros
hist(ranef(m1_Soil_PC)$`Farm`[,1])

summary(m1_Soil_PC)

# 

m1_Soil_PC2 <- lm(PC1_soil_2 ~   
                    Grazing_int_log  +  #  Grazer_diversity +
                    Mowing_frequency +
                    #  Mowing_frq_sqrd +
                    Manuring_freq, data = SEM.dat) 

Anova(m1_Soil_PC2)



plot(m1_Soil_PC2)
qqnorm(resid(m1_Soil_PC2))
qqline(resid(m1_Soil_PC2))

check_convergence(m1_Soil_PC2)

check_collinearity(m1_Soil_PC2)

ranef(m1_Soil_PC2) # random effects are not zeros
hist(ranef(m1_Soil_PC2)$`Farm`[,1])




## PSEM ----
psem_model <- psem (m1_SR_field, 
                    m1_SR_exp, 
                    m2_Abund_exp,
                    m1_Soil_PC,
                    m1_Soil_PC2,
                    Mowing_frequency %~~% Mowing_frq_sqrd,
                    abundance_Exper %~~% SR_Exper_log,
                    humus_log %~~%  PC1_soil_2)



summary(psem_model, .progressBar =F, conserve = TRUE)


tibble(coefs(psem_model)) %>% 
  select(Std.Estimate)

plot(psem_model)
plot(psem_model, digits=2, layout = "dot", 
     node_attrs = list(shape = "rectangle", 
                       color = "black",fillcolor = "white", width=1.4, 
                       distortion=5)) 

