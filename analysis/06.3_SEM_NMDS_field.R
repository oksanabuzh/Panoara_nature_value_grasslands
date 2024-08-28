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
  mutate(Mowing_delay=case_when(Mowing_delay=="no mowing" ~ 0,
                                Mowing_delay=="July-August" ~ 1,
                                Mowing_delay=="June" ~ 2)) %>% 
  #  filter(!Mowing_frequency==3) %>% 
  mutate(Mowing_frq_scal=scale(Mowing_frequency, center=T, scale=T)) %>% 
  mutate(Mowing_frq_sqrd=Mowing_frq_scal^2,
         SR_Exper_log = log1p(SR_Exper+1))




SEM.dat
summary(SEM.dat)
names(SEM.dat)

# NMDS (field data) ----

m1_NMDS1_field <- lmer(NMDS1_VP_field ~   
                       NMDS1_exper +                                                                   
                       NMDS2_exper +
                       Grazing_int_log  +  #  Grazer_diversity +
                       Mowing_frequency +
                       Mowing_frq_sqrd +
                       Manuring_freq + 
                       humus_log + #  PC1_soil + 
                       PC1_soil_2 +
                       # humus_log  + PC1_soil_2 + #PC2_soil + #  +
                       (1|Farm), data = SEM.dat) 


ranef(m1_NMDS1_field) # random effects are not zeros
hist(ranef(m1_SR_field)$`Farm`[,1])

plot(m1_NMDS1_field)
qqnorm(resid(m1_NMDS1_field))
qqline(resid(m1_NMDS1_field))

check_convergence(m1_NMDS1_field)
check_collinearity(m1_NMDS1_field)

Anova(m1_NMDS1_field)
summary(m1_NMDS1_field)

plot_model(m1_NMDS1_field, type = "pred", terms = c("Mowing_frequency"), show.data=T, dot.alpha=0.3, title="") 



## NMDS 2 (field) ----

m1_NMDS2_field <- lmer(NMDS2_VP_field ~   
                      NMDS1_exper +                                                                   
                      NMDS2_exper +
                      Grazing_int_log  +  #  Grazer_diversity +
                      Mowing_frequency +
                      Mowing_frq_sqrd +
                      Manuring_freq + 
                      humus_log + #  PC1_soil + 
                      PC1_soil_2 +
                      # humus_log  + PC1_soil_2 + #PC2_soil + #  +
                      (1|Farm), data = SEM.dat) 


ranef(m1_NMDS2_field) # random effects are not zeros
hist(ranef(m1_NMDS2_field)$`Farm`[,1])

plot(m1_NMDS2_field)
qqnorm(resid(m1_NMDS2_field))
qqline(resid(m1_NMDS2_field))

check_convergence(m1_NMDS2_field)
check_collinearity(m1_NMDS2_field)

Anova(m1_NMDS2_field)
summary(m1_NMDS2_field)


# NMDS (experiment data) ----


m1_NMDS1_exp <- lmer (NMDS1_exper ~ # Shannon_Exper ~   
                     Grazing_int_log  +  # Grazer_diversity +
                   #  Mowing_frequency +
                   #  Mowing_frq_sqrd +
                     Manuring_freq + 
                     humus_log +
                     PC1_soil_2 + Grazing_season +
                     (1|Farm), data = SEM.dat) 

ranef(m1_NMDS1_exp) # random effects are not zeros
hist(ranef(m1_NMDS1_exp)$`Farm`[,1])

m2_NMDS1_exp <- lm (NMDS1_exper ~ # Shannon_Exper ~   
                     Grazing_int_log  +  # Grazer_diversity +
                      # Mowing_frequency +
               #      Mowing_frq_sqrd +
                     Manuring_freq, # + 
                    # humus_log +
                    # PC1_soil_2, 
               data = SEM.dat)

check_collinearity(m2_NMDS1_exp)

Anova(m2_NMDS1_exp)

par(mfrow = c(2, 2))
plot(m2_NMDS1_exp)
par(mfrow = c(1, 1))


## NMDS 2 Experiment
m1_NMDS2_exp <- lmer (NMDS2_exper ~ # Shannon_Exper ~   
                       Grazing_int_log  +  # Grazer_diversity +
                       Mowing_frequency +
                       # Mowing_frq_sqrd +
                       Manuring_freq + 
                        humus_log +
                       PC1_soil_2 +
                       (1|Farm), data = SEM.dat) 

ranef(m1_NMDS2_exp) # random effects are not zeros
hist(ranef(m1_NMDS2_exp)$`Farm`[,1])


m2_NMDS2_exp <- lm(NMDS2_exper ~   
                     Grazing_int_log  +
                     #  Mowing_frequency +
                   #  Mowing_frq_sqrd +
                       Manuring_freq , 
                   #  humus_log + PC1_soil_2, 
                     data = SEM.dat)

check_collinearity(m2_NMDS2_exp)

Anova(m2_NMDS2_exp)

par(mfrow = c(2, 2))
plot(m2_NMDS2_exp)
par(mfrow = c(1, 1))


# Soil PC ----

m1_Soil_PC <- lmer(humus_log ~   
                     Grazing_int_log  +  # Grazer_diversity +
                     Mowing_frequency +
                     #    Mowing_frq_sqrd +
                     Manuring_freq +  
                     (1|Farm), data = SEM.dat) 

ranef(m1_Soil_PC) # random effects are not zeros
hist(ranef(m1_Soil_PC)$`Farm`[,1])


m2_Soil_PC <- lm(humus_log ~   
                   Grazing_int_log  +  # Grazer_diversity +
                   Mowing_frequency +
                   #  Mowing_frq_sqrd +
                   Manuring_freq, data = SEM.dat)

check_collinearity(m2_Soil_PC)


par(mfrow = c(2, 2))
plot(m2_Soil_PC)
par(mfrow = c(1, 1))

Anova(m2_Soil_PC)
summary(m1_Soil_PC)

# 

m1_Soil_PC2 <- lm(PC1_soil_2 ~   
                    Grazing_int_log  +  #  Grazer_diversity +
                    Mowing_frequency +
                   # Mowing_frq_sqrd +
                    Manuring_freq, data = SEM.dat) 


Anova(m1_Soil_PC2)



par(mfrow = c(2, 2))
plot(m1_Soil_PC2)
par(mfrow = c(1, 1))



## PSEM ----
psem_model <- psem (m1_NMDS1_field, m1_NMDS2_field,
                    m2_NMDS1_exp, m2_NMDS2_exp,
                    m2_Soil_PC,
                    m1_Soil_PC2,
                    Mowing_frequency %~~% Mowing_frq_sqrd,
                    humus_log %~~%  PC1_soil_2)



summary(psem_model, .progressBar =F, conserve = TRUE)

tibble(coefs(psem_model)) %>% 
  select(Std.Estimate)

plot(psem_model)
plot(psem_model, digits=2, layout = "dot", 
     node_attrs = list(shape = "rectangle", 
                       color = "black",fillcolor = "white", width=1.4, 
                       distortion=5)) 

