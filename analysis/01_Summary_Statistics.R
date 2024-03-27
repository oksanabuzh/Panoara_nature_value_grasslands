# Data wrangling and summary statistics

library(tidyverse)
library(ggplot2)
library(sjPlot)
library(patchwork)

# Read data, tidying data  ----

Dat <- read_csv("data/Variables_selected.csv") %>% 
  rename(Parcel=`Parcel name`,
         Mowing_method=`Mowing (N no, H hand, M mowing machine, T tractor)`) %>% 
  mutate(Mowing_method=case_when(Mowing_method=="M, T" ~ "M_T" ,
                                 Mowing_method=="T, M"~ "M_T" ,
                                 Mowing_method=="No"~ "no" ,
                                 .default =Mowing_method))
str(Dat)

# Farms & parcels

Dat %>% pull(Farm) %>% unique()
Dat %>% pull(Parcel) %>% unique()

Dat %>% group_by(Farm) %>% count()

# Management variables ----

## Habitat ----
# habitat: the meadows are mowed, grazed and sometimes also manured, 
#          the pastures are only grazed.
Dat %>% pull(habitat)%>% unique()

## Mowing ----
# Mowing_method: no, "H" - hand; "M"- mowing machine; "T" - tractor; "M_T" - mowing machine and tractor
Dat %>% pull(Mowing_method) %>% unique()

Dat %>% select() %>% unique()
