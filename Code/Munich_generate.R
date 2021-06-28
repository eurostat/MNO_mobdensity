library(data.table)
library(tidyverse)

# loading in census data
census.raw <- fread("Data/Zensus_Bevoelkerung_100m-Gitter.csv")

# loading in polygon data from germany to visually define polygon coordinates for the bounding box
germany.raw <- readRDS("Data/gadm36_DEU_1_sf.rds") 

germany <- germany.raw %>%  
  st_transform(crs = 3035)

# define the focus region (bounding box) and execute necessary population manipulation to yield about 1/3 of the population
munich <- census.raw %>% 
  dplyr::select(x = x_mp_100m, y = y_mp_100m, pop.raw = Einwohner) %>% 
  filter(between(y, 2760000, 2800000), 
         between(x, 4420000, 4460000)) %>% # bounding box of munich and near surroundings
  mutate(tile.id = row_number()) %>% 
  mutate(pop.pert = case_when(pop.raw == "-1" | is.na(pop.raw) ~ sample(0:1, n(), replace = T),
                              pop.raw %in% c(2:3) ~ sample(2:3, n(), replace = T),
                              TRUE ~ as.integer(pop.raw))) %>%  # perturbe necessary values according to census specific perturbation rules
  mutate(bin.helper = case_when(pop.pert <= 40 ~ sample(0:1, n(), prob = c(2/3, 1/3), replace = T), # binomial process for tiles with small pop numbers
                                pop.pert > 40 ~  NA_integer_)) %>% 
  mutate(pop = case_when(bin.helper == 1 ~ pop.pert,
                         bin.helper == 0 ~ as.integer(0),
                         is.na(bin.helper) ~ as.integer(round(pop.pert / 3, 0)))) %>%  # reduce population to 1/3 to mimic the mobile phone pop from one MNO
  mutate(elevation = 0)

munich %>% 
  filter(pop.raw > 50) %>% 
  st_as_sf(coords = c("x", "y")) %>% 
  ggplot() +
  geom_sf(aes(color = as.numeric(pop.raw)))

sum(munich$pop.pert)


saveRDS(munich, "Data/munich.rds")



