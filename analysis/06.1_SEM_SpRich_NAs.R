#

SEM.dat <- Dat %>% filter(!Parcel_name=="Brade_1") %>% 
 # mutate(SR_D_E_exper=case_when(is.na(SR_D_E_exper) ~ 0, .default=SR_D_E_exper),
#     Abund_D_E_exper=case_when(is.na(Abund_D_E_exper) ~ 0, .default=Abund_D_E_exper)) %>% 
          filter(!is.na(SR_D_E_exper)) %>% 
  mutate(Grazing_intensity=sqrt(Grazing_intensity_A+1)) %>% 
 # mutate(Grazer_diversity=Grazer_type) %>% 
  mutate(Mowing_delay=case_when(Mowing_delay=="no mowing" ~ 0,
                                Mowing_delay=="July-August" ~ 1,
                                Mowing_delay=="June" ~ 2)) 
  



# Species Richness (field data) ----

m1_SR_field <- glmer(Plant_SR_vascular ~   
                       Abund_D_E_exper +   SR_D_E_exper +
                       Grazing_intensity  +  Grazer_diversity +
                       Mowing_frequency + # Mowing_delay + # Mowing_method + 
                       Excrements_cover +  Manuring_freq  +
                       (1|Farm),  family = "poisson", 
                     data = SEM.dat) 

check_convergence(m1_SR_field)
check_collinearity(m1_SR_field)

sum(residuals(m1_SR_field, type = "pearson")^2) / df.residual(m1_SR_field)


Anova(m1_SR_field)
summary(m1_SR_field)

ranef(m1_SR_field) # random effects are not zeros
hist(ranef(m1_SR_field)$`Farm`[,1])

m2_SR_field <- glm(Plant_SR_vascular ~   
                       Abund_D_E_exper +   SR_D_E_exper +
                       Grazing_intensity  +  Grazer_diversity +
                       Mowing_frequency + # Mowing_delay + #Mowing_method + 
                       Excrements_cover +  Manuring_freq,  family = "poisson", 
                     data = SEM.dat)

anova(m1_SR_field, m2_SR_field)

Anova(m2_SR_field)
check_overdispersion(m2_SR_field)

m3_SR_field <- glm.nb(Plant_SR_vascular ~   
                     Abund_D_E_exper +   SR_D_E_exper +
                     Grazing_intensity  + Grazer_diversity +
                     Mowing_frequency + # Mowing_delay + #Mowing_method + 
                     Excrements_cover +  Manuring_freq, 
                   data = SEM.dat) 

Anova(m3_SR_field)
summary(m3_SR_field)

# Species Richness (experiment data) ----

m1_SR_exp <- glmer (SR_D_E_exper ~   Abund_D_E_exper +
                    Grazing_intensity  +   Grazer_diversity +  
                    Excrements_cover + Manuring_freq +
                      (1|Farm), family = "poisson", 
                  data = SEM.dat) 

check_convergence(m1_SR_exp)

check_collinearity(m1_SR_exp)
sum(residuals(m1_SR_exp, type = "pearson")^2) / df.residual(m1_SR_exp)

ranef(m1_SR_exp) # random effects are not zeros
hist(ranef(m1_SR_exp)$`Farm`[,1])

m2_SR_exp <- glm (SR_D_E_exper ~   Abund_D_E_exper +
                      Grazing_intensity  +   Grazer_diversity +  
                      Excrements_cover + Manuring_freq,
                  family = "poisson", 
                    data = SEM.dat)

summary(m2_SR_exp)
check_overdispersion(m2_SR_exp)
Anova(m2_SR_exp)


# Abundance (Experiment data) ----
m1_Abund_exp <- glmer (Abund_D_E_exper ~  Grazing_intensity  +  Grazer_diversity +  
                         Excrements_cover + Manuring_freq + 
                         (1|Farm), family = "poisson", 
                       data = SEM.dat) 

check_convergence(m1_Abund_exp)

Anova(m1_Abund_exp)
summary(m1_Abund_exp)

ranef(m1_Abund_exp) # random effects are not zeros
hist(ranef(m1_Abund_exp)$`Farm`[,1])

m2_Abund_exp <- glm (Abund_D_E_exper ~  Grazing_intensity  +  Grazer_diversity +  
                       Excrements_cover + Manuring_freq, family = "poisson", 
                       data = SEM.dat) 

anova(m1_Abund_exp, m2_Abund_exp)

Anova(m1_Abund_exp)


sum(residuals(m1_Abund_exp, type = "pearson")^2) / df.residual(m1_Abund_exp)
check_overdispersion(m1_Abund_exp)


m3_Abund_exp <- glmer.nb(Abund_D_E_exper ~  Grazing_intensity  +  Grazer_diversity +  
                          Excrements_cover + Manuring_freq  + (1|Farm),
                        data=SEM.dat) 

check_convergence(m3_Abund_exp)
check_overdispersion(m3_Abund_exp)
check_collinearity(m3_Abund_exp)

m4_Abund_exp <- glm.nb(Abund_D_E_exper ~  Grazing_intensity  +  Grazer_diversity +  
                           Excrements_cover + Manuring_freq ,
                         data=SEM.dat) 

anova(m3_Abund_exp, m4_Abund_exp)

Anova(m4_Abund_exp)
# or
m5_Abund_exp <- glm(Abund_D_E_exper ~  Grazing_intensity  +  Grazer_diversity +  
                         Excrements_cover + Manuring_freq , family="quasipoisson",
                       data=SEM.dat)
Anova(m5_Abund_exp)

# Excrement cover ----

m1_Excrem_cover <- lmer(Excrements_cover  ~                      
                          Grazing_intensity   + # Mowing_frequency + # 
                    (1|Farm), 
                  data = SEM.dat )

plot(m1_Excrem_cover)
qqnorm(resid(m1_Excrem_cover))
qqline(resid(m1_Excrem_cover))

Anova(m1_Excrem_cover)


ranova(m1_Excrem_cover)

m2_Excrem_cover <- lm(Excrements_cover ~                      
                          Grazing_intensity + Mowing_frequency, #+ Grazer_type , 
                        data = SEM.dat )
# plot(m2_Excrem_cover)
Anova(m2_Excrem_cover)

## PSEM ----
psem_model <- psem (m3_SR_field, 
                    m2_SR_exp, m4_Abund_exp,  
                    m2_Excrem_cover
                   # SR_D_E_exper %~~% Abund_D_E_exper,
                    #Grazing_intensity  %~~%  Grazer_diversity  
                   )



summary(psem_model, .progressBar =F)
coefs(psem_model)

plot(psem_model)
plot(psem_model, digits=2, layout = "dot", 
     node_attrs = list(shape = "rectangle", 
                       color = "black",fillcolor = "white", width=1.4, 
                       distortion=5)) 

