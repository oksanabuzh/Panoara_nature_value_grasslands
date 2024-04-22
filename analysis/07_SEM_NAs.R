#

SEM.dat <- Dat %>% filter(!Parcel_name=="Brade_1") %>% 
  # mutate(SR_D_E_exper=case_when(is.na(SR_D_E_exper) ~ 0, .default=SR_D_E_exper),
  #      Abund_D_E_exper=case_when(is.na(Abund_D_E_exper) ~ 0, .default=Abund_D_E_exper)) %>% 
            filter(!is.na(SR_D_E_exper)) %>% 
  mutate(Grazing_intensity=sqrt(Grazing_intensity_A+1)) 



#m1_SR_field

m1_SR_field <- glmer(Plant_SR_vascular ~   
                       Abund_D_E_exper +   SR_D_E_exper +
                       Grazing_intensity  +  Grazer_diversity +
                       Mowing_frequency + # Mowing_delay + #Mowing_method + 
                       Excrements_cover +  Manuring_freq  +
                       (1|Farm),  family = "poisson", 
                     data = SEM.dat) 


check_convergence(m1_SR_field)
check_collinearity(m1_SR_field)
Anova(m1_SR_field)
summary(m1_SR_field)
sum(residuals(m1_SR_field, type = "pearson")^2) / df.residual(m1_SR_field)


m2_SR_field <- glmer.nb(Plant_SR_vascular ~   
                       Abund_D_E_exper +   SR_D_E_exper +
                       Grazing_intensity  +  Grazer_diversity +
                       Mowing_frequency + # Mowing_delay + #Mowing_method + 
                       Excrements_cover +  Manuring_freq  +
                       (1|Farm),  family = "poisson", 
                     data = SEM.dat) 


check_convergence(m2_SR_field)
check_collinearity(m2_SR_field)

ranef(m2_SR_field) # random effects are not zeros
hist(ranef(m2_SR_field)$`Farm`[,1])

Anova(m2_SR_field)
summary(m2_SR_field)
sum(residuals(m2_SR_field, type = "pearson")^2) / df.residual(m2_SR_field)


m1_SR_exp <- glmer (SR_D_E_exper ~   Abund_D_E_exper +
                      Grazing_intensity  +  Grazer_diversity +  
                      Excrements_cover + Manuring_freq +
                      (1|Farm), family = "poisson", data = SEM.dat) 


ranef(m1_SR_exp) # random effects are not zeros
hist(ranef(m1_SR_exp)$`Farm`[,1])

Anova(m1_SR_exp)
summary(m1_SR_exp)


m2_SR_exp <- glm (SR_D_E_exper ~   Abund_D_E_exper +
                      Grazing_intensity  +  Grazer_diversity +  
                      Excrements_cover + Manuring_freq, family = "poisson", data = SEM.dat) 

summary(m2_SR_exp)
38.994/27

Anova(m2_SR_exp)


m1_Abund_exp <- glmer (Abund_D_E_exper ~  Grazing_intensity  +  Grazer_diversity +  
                         Excrements_cover + Manuring_freq + 
                        Mowing_frequency +
                         (1|Farm), family = "poisson", 
                       data = SEM.dat) 

ranef(m1_Abund_exp) # random effects are not zeros
hist(ranef(m1_Abund_exp)$`Farm`[,1])

check_convergence(m1_Abund_exp)
Anova(m1_Abund_exp)
summary(m1_Abund_exp)


check_overdispersion(m1_Abund_exp)
sum(residuals(m1_Abund_exp, type = "pearson")^2) / df.residual(m1_Abund_exp)


m2_Abund_exp <- glm (Abund_D_E_exper ~  Grazing_intensity  +  Grazer_diversity +  
                       Excrements_cover + Manuring_freq + 
                       Mowing_frequency, family = "poisson", 
                     data = SEM.dat) 

anova(m1_Abund_exp, m2_Abund_exp)

library(glmmTMB)

m3_Abund_exp <- glmer.nb (Abund_D_E_exper ~  Grazing_intensity  +  Grazer_diversity +  
                            Excrements_cover + Manuring_freq + 
                            Mowing_frequency + (1|Farm), 
                          data = SEM.dat) 
check_convergence(m3_Abund_exp)
Anova(m3_Abund_exp)
summary(m3_Abund_exp)

m3_Abund_exp <- glmer.nb(Abund_D_E_exper ~  Grazing_intensity  +  Grazer_diversity +  
                          Excrements_cover + Manuring_freq + 
                          Mowing_frequency + (1|Farm),
                         data=SEM.dat,
                         ziformula=~0 ,
                         family=nbinom1) #quasipoisson

check_overdispersion(m3_Abund_exp)
Anova(m3_Abund_exp)


#

m4_Abund_exp <- lmer(log(Abund_D_E_exper+1) ~  Grazing_intensity  +  Grazer_diversity +  
                            Excrements_cover + Manuring_freq + 
                            Mowing_frequency + (1|Farm), 
                          data = SEM.dat) 
check_convergence(m4_Abund_exp)
Anova(m4_Abund_exp)
summary(m4_Abund_exp)

plot(m4_Abund_exp)
qqnorm(resid(m4_Abund_exp))
qqline(resid(m4_Abund_exp))

ranova(m4_Abund_exp)

m5_Abund_exp <- lm(log(Abund_D_E_exper+1) ~  Grazing_intensity  +  Grazer_diversity +  
                       Excrements_cover + Manuring_freq + 
                       Mowing_frequency , 
                     data = SEM.dat)
plot(m5_Abund_exp)
anova(m4_Abund_exp, m5_Abund_exp)

Anova(m5_Abund_exp)

 #-------------------------------------------------------------------------------#

m1_Excrem_cover <- lmer(Excrements_cover~                      
                          Grazing_intensity  +  # Grazer_type + 
                          (1|Farm), 
                        data = SEM.dat )


psem_model <- psem (m1_SR_field, 
                    m3_Abund_exp,  m1_SR_exp,
                    m1_Excrem_cover
                    # SR_D_E_exper %~~% Abund_D_E_exper,
                    #Grazing_intensity  %~~%  Grazer_diversity  
)



summary(psem_model, .progressBar =F)



coefs(psem_model)

plot(psem_model)

# nicely, circle, tree, kk, and fr.

x11(height=19,width=18.5)

plot(psem_model, digits=2, layout = "dot", 
     node_attrs = list(shape = "rectangle", 
                       color = "black",fillcolor = "white", width=1.4, 
                       distortion=5)) 

