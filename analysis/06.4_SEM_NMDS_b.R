#

SEM.dat <- Dat %>% filter(!Parcel_name=="Brade_1") %>% 
  mutate(SR_D_E_exper=case_when(is.na(SR_D_E_exper) ~ 0, .default=SR_D_E_exper),
       Abund_D_E_exper=case_when(is.na(Abund_D_E_exper) ~ 0, .default=Abund_D_E_exper)) %>% 
#          filter(!is.na(SR_D_E_exper)) %>% 
  mutate(Grazing_intensity=sqrt(Grazing_intensity_A+1))



# NMDS1_VP_field  (field data) ----

m1_NMDS1_field <- lmer(NMDS1_VP_field ~   
                         NMDS1_exper +   NMDS2_exper +
                       Grazing_intensity  +  Grazer_diversity +
                       Mowing_frequency + # Mowing_delay + #Mowing_method + 
                       Excrements_cover +  Manuring_freq  +
                       (1|Farm), data = SEM.dat) 

check_convergence(m1_NMDS1_field)
check_collinearity(m1_NMDS1_field)

ranef(m1_NMDS1_field) # random effects are not zeros
hist(ranef(m1_NMDS1_field)$`Farm`[,1])
ranova(m1_NMDS1_field)

summary(m1_NMDS1_field)

plot(m1_NMDS1_field)
qqnorm(resid(m1_NMDS1_field))
qqline(resid(m1_NMDS1_field))

Anova(m1_NMDS1_field)

# NMDS2_VP_field  (field data) ----

m1_NMDS2_field <- lmer(NMDS2_VP_fiels ~   
                         NMDS1_exper +   NMDS2_exper +
                         Grazing_intensity  +  Grazer_diversity +
                         Mowing_frequency + # Mowing_delay + #Mowing_method + 
                         Excrements_cover +  Manuring_freq  +
                         (1|Farm), data = SEM.dat) 

check_convergence(m1_NMDS2_field)
check_collinearity(m1_NMDS2_field)

summary(m1_NMDS2_field)

plot(m1_NMDS2_field)
qqnorm(resid(m1_NMDS2_field))
qqline(resid(m1_NMDS2_field))

Anova(m1_NMDS2_field)

ranef(m1_NMDS2_field) # random effects are not zeros
hist(ranef(m1_NMDS2_field)$`Farm`[,1])
ranova(m1_NMDS2_field)

m2_NMDS2_field <- lm(NMDS2_VP_fiels ~   
                         NMDS1_exper +   NMDS2_exper +
                         Grazing_intensity  +  Grazer_diversity +
                         Mowing_frequency + # Mowing_delay + #Mowing_method + 
                         Excrements_cover +  Manuring_freq, data = SEM.dat)


plot(m2_NMDS2_field)
Anova(m2_NMDS2_field)

# NMDS1_exper  (experiment data) ----

m1_NMDS1_exper <- lmer(NMDS1_exper ~   
                         Grazing_intensity  +  Grazer_diversity +
                         Mowing_frequency + # Mowing_delay + #Mowing_method + 
                         Excrements_cover +  Manuring_freq  +
                         (1|Farm), data = SEM.dat) 

check_convergence(m1_NMDS1_exper)
check_collinearity(m1_NMDS1_exper)

plot(m1_NMDS1_exper)
qqnorm(resid(m1_NMDS1_exper))
qqline(resid(m1_NMDS1_exper))

ranef(m1_NMDS1_exper) # random effects are  zeros
hist(ranef(m1_NMDS1_exper)$`Farm`[,1])
ranova(m1_NMDS1_exper)

m2_NMDS1_exper <- lm(NMDS1_exper ~   
                         Grazing_intensity  +  Grazer_diversity +
                         Mowing_frequency + # Mowing_delay + #Mowing_method + 
                         Excrements_cover +  Manuring_freq, data = SEM.dat) 

summary(m2_NMDS1_exper)
plot(m2_NMDS1_exper)
Anova(m2_NMDS1_exper)



# NMDS2_exper  (experiment data) ----

m1_NMDS2_exper <- lmer(NMDS2_exper ~   
                         Grazing_intensity  +  Grazer_diversity +
                         Mowing_frequency + # Mowing_delay + #Mowing_method + 
                         Excrements_cover +  Manuring_freq  +
                         (1|Farm), data = SEM.dat) 

check_convergence(m1_NMDS2_exper)
check_collinearity(m1_NMDS2_exper)

plot(m1_NMDS2_exper)
qqnorm(resid(m1_NMDS2_exper))
qqline(resid(m1_NMDS2_exper))

ranef(m1_NMDS2_exper) # random effects are  zeros
hist(ranef(m1_NMDS2_exper)$`Farm`[,1])
ranova(m1_NMDS2_exper)

m2_NMDS2_exper <- lm(NMDS2_exper ~   
                       Grazing_intensity  +  Grazer_diversity +
                       Mowing_frequency + # Mowing_delay + #Mowing_method + 
                       Excrements_cover +  Manuring_freq, data = SEM.dat) 

summary(m2_NMDS2_exper)
plot(m2_NMDS2_exper)
Anova(m2_NMDS2_exper)



# Excrement cover ----

m1_Excrem_cover <- lmer(Excrements_cover  ~                      
                          Grazing_intensity  +  # Mowing_frequency + # Grazer_type + 
                    (1|Farm), 
                  data = SEM.dat )

plot(m1_Excrem_cover)
qqnorm(resid(m1_Excrem_cover))
qqline(resid(m1_Excrem_cover))

Anova(m1_Excrem_cover)


ranova(m1_Excrem_cover)

m2_Excrem_cover <- lm(Excrements_cover ~                      
                          Grazing_intensity, 
                        data = SEM.dat )
# plot(m2_Excrem_cover)
Anova(m2_Excrem_cover)


psem_model <- psem (m1_NMDS1_field, m2_NMDS2_field,
                    m2_NMDS1_exper, m2_NMDS2_exper)



summary(psem_model, .progressBar =F)


plot(psem_model)
plot(psem_model, digits=2, layout = "dot", 
     node_attrs = list(shape = "rectangle", 
                       color = "black",fillcolor = "white", width=1.4, 
                       distortion=5)) 

