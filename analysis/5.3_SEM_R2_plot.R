# Purpose: Create R2 plot for the SEM (Fig. 4E)

library(tidyverse)

# set theme for plots
set_theme(base = theme_bw(),
          axis.textsize.x = 0.9, axis.textsize.y = 0.9, axis.textcolor = "black",
          axis.title.color = "black", axis.title.size = 1.2,
          geom.linetype = 1) #legend.pos = "None", 



# read data
R2_SR_partial <- read_csv("results/SEM_SR_R2_partial.csv")
R2_NMDS1_partial <- read_csv("results/SEM_NMDS_R2_partial.csv")

R2_SR_partial %>% 
  distinct(Effect)

R2_partial <- R2_SR_partial %>% 
  bind_rows(R2_NMDS1_partial) %>% 
  dplyr::select(Effect, responce, Rsq) %>% 
  filter(!Effect=="Model") %>% 
  mutate(predictor=recode_factor(Effect,
                                 Grazing_int_log = "Grazing intensity",
                                 Manuring_freq = "Manuring frequency",
                                 Mowing_frequency = "Mowing frequency",
                                 Mowing_frq_sqrd = "Mowing frequency",
                                 
                                 abundance_Exper = "Seedling abundance",
                                 SR_Exper_log = "Seedling species richness",
                                 NMDS1_exper = "Seedling composition, NMDS1",
                                 NMDS2_exper = "Seedling composition, NMDS2",
                                 humus_log = "Humus",
                                 soil_NPK = "Soil NPK"))%>% 
  mutate(predictor_group=recode_factor(Effect,
                                 Grazing_int_log = "Grazing",
                                 Manuring_freq = "Manuring",
                                 Mowing_frequency = "Mowing",
                                 Mowing_frq_sqrd = "Mowing",
                                 
                                 abundance_Exper = "Seed dispersal",
                                 SR_Exper_log = "Seed dispersal",
                                 NMDS1_exper = "Seed dispersal",
                                 NMDS2_exper = "Seed dispersal",
                                 humus_log = "Soil properties",
                                 soil_NPK = "Soil properties"))%>% 
  arrange(predictor) %>% 
  group_by(predictor, predictor_group, responce) %>% 
  summarise(Rsq=sum(Rsq)) %>% 
  ungroup() %>% 
  mutate(responce=recode_factor(responce,
                                SR = "Species richness",
                                NMDS1 = "NMDS 1",
                                NMDS2 = "NMDS 2"))
 # mutate(responce=fct_relevel(responce, c("SR", "NMDS1", "NMDS2")))
  

R2_partial



## > Plot R2 ----

dodge_width <- 0.5

ggplot(R2_partial,
               aes(y =reorder(predictor, desc(predictor)), x = Rsq, 
                   color = predictor_group)) +
  geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
  geom_point(position = position_dodge(width = dodge_width), size = 4) +
  geom_errorbarh(aes(xmin = 0, xmax = Rsq),
                 position = position_dodge(width = dodge_width), height = 0.1) +
  theme(legend.key=element_blank()) +
  # theme_bw()+
  xlim(-0.03,0.54)+
  facet_wrap(~responce) +
  scale_color_manual(values=c("#F8766D", "#9BB018","orange", "#C77CFF", "#BA7938"))+ 
  # MetBrewer::scale_color_met_d("Kandinsky") +
  theme(axis.text.x = element_text(size=8), axis.text.y = element_text(size=10),
        axis.title=element_text(size=13, face="bold"), 
        legend.text = element_text(size=10)) +
  #  theme(axis.title.x=element_text(vjust=-0.1), axis.title.y=element_text(vjust=2))+
  labs(y="", color="Driver type",
       x= expression(paste("Variance explained, ", R^{2})))

