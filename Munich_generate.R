library(data.table)
library(tidyverse)


census.raw <- fread("Data/Zensus_Bevoelkerung_100m-Gitter.csv")


germany.raw <- readRDS("Data/gadm36_DEU_1_sf.rds") 

germany <- germany.raw %>%  
  st_transform(crs = 3035)

munich <- census.raw %>% 
  dplyr::select(x = x_mp_100m, y = y_mp_100m, pop.raw = Einwohner) %>% 
  filter(between(y, 2760000, 2800000), 
         between(x, 4420000, 4460000)) %>%
  mutate(tile.id = row_number()) %>% 
  mutate(pop = case_when(pop.raw == "-1" | is.na(pop.raw) ~ sample(0:1, n(), replace = T),
                         pop.raw %in% c(2:3) ~ sample(2:3, n(), replace = T),
                         TRUE ~ as.integer(pop.raw))) %>% 
  mutate(elevation = 0)
  
munich %>% 
  filter(pop.raw > 50) %>% 
  st_as_sf(coords = c("x", "y")) %>% 
  ggplot() +
  geom_sf(aes(color = as.numeric(pop.raw)))


saveRDS(munich, "Data/munich.rds")



