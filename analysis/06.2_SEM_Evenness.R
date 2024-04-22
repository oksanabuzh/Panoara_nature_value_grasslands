#

SEM.dat <- Dat %>% filter(!Parcel_name=="Brade_1") %>% 
  mutate(SR_D_E_exper=case_when(is.na(SR_D_E_exper) ~ 0, .default=SR_D_E_exper),
       Abund_D_E_exper=case_when(is.na(Abund_D_E_exper) ~ 0, .default=Abund_D_E_exper)) %>% 
#          filter(!is.na(SR_D_E_exper)) %>% 
  mutate(Grazing_intensity=sqrt(Grazing_intensity_A+1))


# Correlation between SR and evenness

m1 <- lmer(EvennessVP_field ~   Plant_SR_vascular + (1|Farm), data = SEM.dat) 
ranova(m1)
Anova(m1)

m1_pred <- get_model_data(m1,type = "pred", terms="Plant_SR_vascular")

ggplot(m1_pred, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1)+
  geom_point(data=SEM.dat,
             aes(Plant_SR_vascular, EvennessVP_field), 
             size=3, alpha=0.7, pch=21, fill="#64ABCE") +
#  scale_fill_viridis(discrete=TRUE, option = "D")+
  labs(y="Evenness", x='Species richness')+
  geom_line(linetype=1, linewidth=1) 




# EvennessVP_field  (field data) ----

m1_Evenness_field <- lmer(EvennessVP_field ~   
                       Abund_D_E_exper +   SR_D_E_exper +
                       Grazing_intensity  +  Grazer_diversity +
                       Mowing_frequency + # Mowing_delay + #Mowing_method + 
                       Excrements_cover +  Manuring_freq  +
                       (1|Farm), data = SEM.dat) 

check_convergence(m1_Evenness_field)
check_collinearity(m1_Evenness_field)

ranef(m1_Evenness_field) # random effects are not zeros
hist(ranef(m1_Evenness_field)$`Farm`[,1])
ranova(m1_Evenness_field)

summary(m1_Evenness_field)

plot(m1_Evenness_field)
qqnorm(resid(m1_Evenness_field))
qqline(resid(m1_Evenness_field))

Anova(m1_Evenness_field)



# Species Richness (experiment data) ----

m1_SR_exp <- glmer (SR_D_E_exper ~   Abund_D_E_exper +
                    Grazing_intensity  +  Grazer_diversity +  
                    Excrements_cover + Manuring_freq +
                      (1|Farm), family = "poisson", 
                  data = SEM.dat) 

check_convergence(m1_SR_exp)

check_collinearity(m1_SR_exp)

ranef(m1_SR_exp) # random effects are not zeros
hist(ranef(m1_SR_exp)$`Farm`[,1])

sum(residuals(m1_SR_exp, type = "pearson")^2) / df.residual(m1_SR_exp)

summary(m1_SR_exp)
Anova(m1_SR_exp)


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
                        data=SEM.dat , # to include zero inflation use ziformula=~1
                        family=nbinom2) #quasipoisson

check_convergence(m3_Abund_exp)

check_overdispersion(m3_Abund_exp)

Anova(m3_Abund_exp)
summary(m3_Abund_exp)

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


psem_model <- psem (m1_Evenness_field, 
                    m1_SR_exp, m3_Abund_exp,  
                    m1_Excrem_cover
                   # SR_D_E_exper %~~% Abund_D_E_exper,
                    #Grazing_intensity  %~~%  Grazer_diversity  
                   )



summary(psem_model, .progressBar =F)


plot(psem_model)
plot(psem_model, digits=2, layout = "dot", 
     node_attrs = list(shape = "rectangle", 
                       color = "black",fillcolor = "white", width=1.4, 
                       distortion=5)) 

