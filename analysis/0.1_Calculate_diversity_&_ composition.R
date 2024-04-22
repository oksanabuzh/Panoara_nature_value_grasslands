# Purpose: Calculate evenness and community composition (NMDS) for each plot

# (1) Calculate evenness----

## 1.1 Field data----

Community_VP_field <- read_csv("data/Community_composition_VegetationPlots.csv") %>% 
  rename(layer="layer (VP - vascular, B - bryophyte, L - lichen)") %>% 
  filter(layer == "VP") %>% # only for vascular plants 
  dplyr::select(-layer) %>% 
  pivot_longer(!Plot_Parcel, names_to = "Parcel_name", values_to = "cover") %>% 
  rename(Species="Plot_Parcel") %>% 
  mutate(Parcel_name=str_replace_all(Parcel_name, " ", "_")) %>% 
  filter(!is.na(cover))
  
Community_VP_field %>% pull(Parcel_name)%>% unique()

str(Community_VP_field)

VP_field <- Community_VP_field %>% 
  group_by(Parcel_name) %>% 
  summarise(
    CoverVP_field=sum(cover),
    SR_VP_field = n_distinct(Species),
    EvennessVP_field = vegan::diversity(cover, index = "invsimpson")
    ) %>%
  ungroup()

VP_field


## 1.2 Experiment data----

Community_exper <- read_csv("data/Community_composition_DungExperiment.csv") %>% 
  rename(layer="layer (VP - vascular, B - bryophyte, L - lichen)") %>% 
  filter(layer == "VP") %>% # only for vascular plants 
  dplyr::select(-layer) %>% 
  pivot_longer(!Plot_Parcel, names_to = "Parcel_name", values_to = "abundance") %>% 
  rename(Species="Plot_Parcel") %>% 
  mutate(Parcel_name=str_replace_all(Parcel_name, " ", "_")) %>% 
  filter(!is.na(abundance))

Community_exper %>% pull(Parcel_name)%>% unique()

str(Community_exper)

VP_Exper <- Community_exper %>% 
  group_by(Parcel_name) %>% 
  summarise(
    abundance_Exper=sum(abundance),
    SR_Exper = n_distinct(Species),
    Evenness_Exper = vegan::diversity(abundance, index = "invsimpson")
  ) %>%
  ungroup()

VP_Exper


Diver_data <-  VP_field %>% 
  left_join(VP_Exper, by="Parcel_name")

#(2) NMDS ----

## 2.1 Field data----

compos_VP_field <- Community_VP_field %>% 
  pivot_wider(names_from = "Species", values_from = "cover") %>%
  replace(is.na(.), 0) 

# check if NA
which(is.na(compos_VP_field %>% 
              dplyr::select(-Parcel_name)))

#Check if any species with 0 cover 
which(colSums(compos_VP_field%>% 
                dplyr::select(-Parcel_name),
              na.rm=TRUE) %in% 0)

# NMDS analysis
nmds1<-vegan::metaMDS(compos_VP_field%>% 
                        dplyr::select(-Parcel_name), 
                      distance="bray",k=2,trymax=100)
nmds1 # the stress value shows how easy it was to condense multidimensional data into two dimensional space, below 0.2 is generally good
vegan::scores(nmds1)

as.data.frame(vegan::scores(nmds1)$sites)
as.data.frame(vegan::scores(nmds1)$species)

which(is.na(as.data.frame(vegan::scores(nmds1)$species)))

NMDS_VP_field <- compos_VP_field %>% 
  dplyr::select(Parcel_name) %>% 
  mutate(NMDS1_VP_field = as.data.frame(vegan::scores(nmds1)$sites)$NMDS1,
         NMDS2_VP_fiels = as.data.frame(vegan::scores(nmds1)$sites)$NMDS2)


Diver_NMDS_data <- Diver_data %>% 
  left_join(NMDS_VP_field, by="Parcel_name")

## 2.2 Experiment data----

compos_exper <- Community_exper %>% 
  pivot_wider(names_from = "Species", values_from = "abundance") %>%
  replace(is.na(.), 0) 

# check if NA
which(is.na(compos_exper %>% 
              dplyr::select(-Parcel_name)))

#Check if any species with 0 cover 
which(colSums(compos_exper%>% 
                dplyr::select(-Parcel_name),
              na.rm=TRUE) %in% 0)

# NMDS analysis
nmds2<-vegan::metaMDS(compos_exper%>% 
                        dplyr::select(-Parcel_name), 
                      distance="bray",k=2,trymax=100)
nmds2 # the stress value shows how easy it was to condense multidimensional data into two dimensional space, below 0.2 is generally good
vegan::scores(nmds2)

as.data.frame(vegan::scores(nmds2)$sites)
as.data.frame(vegan::scores(nmds2)$species)

which(is.na(as.data.frame(vegan::scores(nmds2)$species)))

NMDS_exper <- compos_exper %>% 
  dplyr::select(Parcel_name) %>% 
  mutate(NMDS1_exper = as.data.frame(vegan::scores(nmds2)$sites)$NMDS1,
         NMDS2_exper = as.data.frame(vegan::scores(nmds2)$sites)$NMDS2)


Diver_NMDS_data <- Diver_data %>% 
  left_join(NMDS_VP_field, by="Parcel_name") %>% 
  left_join(NMDS_exper, by="Parcel_name")


write_csv(Diver_NMDS_data, "data/Diversity_&_NMDS_data.csv")
