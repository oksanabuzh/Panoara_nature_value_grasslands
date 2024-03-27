# Data wrangling and summary statistics

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

# Read data, tidying data  ----

Dat <- read_csv("data/Variables_selected.csv") %>% 
  rename(Parcel=`Parcel name`,
         Plant_SR_total = "number of all plants",
         Plant_SR_vascular = "number of vascular plants_10 m2",
         Plant_SR_cryptogams = "number of cryptogams_10 m2",
         Plant_SR_bryophytes ="number of bryophytes_10 m2",
         Plant_SR_lichens ="number of lichens_10 m2",
         Mowing_method=`Mowing (N no, H hand, M mowing machine, T tractor)`,
         Mowing_frequency="Mowing frequency per year",
         Mowing_delay="Date of the first cut",
         Grazer_type_specific="Grazing animal numbers (C cow, S sheep, G goar, H horse)",
         Grazing_season = "Grazing season (No, A, SA, W=whole season)") %>% 
  mutate(Mowing_method=case_when(Mowing_method=="M, T" ~ "M_T" ,
                                 Mowing_method=="T, M"~ "M_T" ,
                                 Mowing_method=="No"~ "no" ,
                                 .default =Mowing_method),
         ownership=factor(ownership),
         Farm=factor(Farm),
         Parcel_name=Parcel) %>% 
  separate_wider_delim(Parcel, " ", names=c("Farm_name", NA)) # Does this makes sense?

str(Dat)
names(Dat)

# Farms & parcels

Dat %>% pull(Farm) %>% unique()
Dat %>% pull(Farm_name) %>% unique()

Dat %>% pull(Parcel_name) %>% unique()
Dat %>% pull(ownership) %>% unique()

Dat %>% group_by(Farm) %>% count()

Dat %>% group_by(ownership, Farm) %>% count()
Dat %>% group_by(ownership, Farm_name) %>% count()


# Management variables ----

## Habitat ----
# habitat: the meadows are mowed, grazed and sometimes also manured, 
#          the pastures are only grazed.
Dat %>% pull(habitat)%>% unique()

## Mowing ----
# Mowing_method: no, "H" - hand; "M"- mowing machine; "T" - tractor; "M_T" - mowing machine and tractor
Dat %>% pull(Mowing_method) %>% unique()

Dat %>% group_by(habitat, Mowing_method) %>% count()

# Mowing_frequency
Dat %>% pull(Mowing_frequency)

Dat %>% group_by(habitat, Mowing_method) %>% 
  summarise(n = n(),
            min_Mowing_frequency=min(Mowing_frequency),
            max=max(Mowing_frequency),
            mean=mean(Mowing_frequency),
            sd=sd(Mowing_frequency)) 
  

Dat %>% mutate(Mowing_frequency=factor(Mowing_frequency)) %>% 
  group_by(habitat, Mowing_method, Mowing_frequency) %>% 
  count()

# Mowing_delay
Dat %>% pull(Mowing_delay) %>% unique()

Dat %>% mutate(Mowing_frequency=factor(Mowing_frequency)) %>% 
  group_by(Mowing_delay, Mowing_frequency, Mowing_method, habitat) %>% 
  count()


# Check of random effects & rank deficiency of fixed effcets
# we need to use Poisson for Plant_SR_total and random effects of "Farm" 
# For now I use lm / lmer  only to see the rank deficiency of the fixed effects 

m1 <- lm(
  Plant_SR_total ~  habitat + 
    Mowing_delay + Mowing_frequency + Mowing_method,
    data = Dat)

check_collinearity(m1)

summary(m1)

m1 <- lm(
  Plant_SR_total ~  habitat + 
    Mowing_delay + Mowing_frequency + Mowing_method,
  data = Dat %>% 
    filter(!Mowing_method=="no"))

check_collinearity(m1)

car::Anova(m1)
anova(m1)


# check the random effects

m2 <- lmer(
  Plant_SR_total ~ habitat + 
    Mowing_delay +  Mowing_frequency + # Mowing_method + 
      (1 | Farm), data = Dat%>% 
    filter(!Mowing_method=="no"))

# random effects
lmerTest::ranova(m2)


check_collinearity(m2)

summary(m2)


car::Anova(m2)

# Model assumptions
plot(m2) 
qqnorm(resid(m2))
qqline(resid(m2))

# Grazing ----

Dat %>% pull(Grazer_type)%>% unique()

Dat %>% pull(Grazer_type)%>% unique()

