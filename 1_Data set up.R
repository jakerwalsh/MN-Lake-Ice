# Long-term trends in lake ice
# Hansen Lab Project (Revisit 2023)

# Set up repository folders ----
ifelse(!dir.exists("data"), 
       dir.create("data"), 
       "Folder exists already")
ifelse(!dir.exists("figures"), 
       dir.create("figures"), 
       "Folder exists already")
ifelse(!dir.exists("tables"), 
       dir.create("tables"), 
       "Folder exists already")
ifelse(!dir.exists("exploration"), 
       dir.create("exploration"), 
       "Folder exists already")
ifelse(!dir.exists("outputdata"), 
       dir.create("outputdata"), 
       "Folder exists already")

# Data set up ----

## Packages ----
# Apologies - some of these are left over from old scripts, and we may not use them.

# The "mnsentinellakes" package is downloaded through github using the package "remotes":
# install.packages("remotes")
# remotes::install_github("mnsentinellakes/mnsentinellakes")

# Data manipulation, plotting
library(tidyverse)
library(viridis)
library(gridExtra)
library(grid)
library(broom)
library(cowplot)
library(scales)

library("googledrive")
library(data.table)
library(RCurl)

# Nice for dealing with working directories
library(here)

# Time series breakpoint analysis
library(strucchange)
library(segmented)
library(trend)# Includes a function for Sen's Slope
library(pracma)

# Linear mixed effects modeling
library(lme4)

# GAM
library(mgcv)
library(gratia)
library(gamm4)
library(itsadug)

library(lubridate)
library(parallel)
library(mnsentinellakes)

library(sf)
library(maps)

# Wavelet analysis
library(vectorwavelet)
library(wavelets)
library(WaveletComp)


## Download data from Google Drive ----
# All of these end up in the "data" folder
# Wasn't sure how to download a folder, I think "drive_download" is just for files
drive_download("https://drive.google.com/file/d/12sfEMIu2Mq1TfPGDqz-ShC18-FjhHdKM/view?usp=drive_link",
               path="data/ice_duration_summarized.csv")
drive_download("https://drive.google.com/file/d/12pI3h3fAg8TpKswBFRwyIqksElCymjm6/view?usp=drive_link",
               path="data/Rainy&Kabetogama_ice_duration_summarized.csv")
drive_download("https://drive.google.com/file/d/10vZYWBhdwEH4msV6myq3EAnHgl2KwteD/view?usp=drive_link",
               path="data/ice_on_summarized.csv")
drive_download("https://drive.google.com/file/d/1NO7qQmKOJSZsKOOveukIdMA-n-0WPCsD/view?usp=drive_link",
               path="data/ice_off_summarized.csv")
drive_download("https://drive.google.com/file/d/1BSGrPcYpKYb9bJEOQAsdc8KL8gJV8FdO/view?usp=drive_link",
               path="data/Rainy&Kabetogama_ice_off_summarized.csv")
drive_download("https://drive.google.com/file/d/1y5l_Larn97P5rpHm_cK-DJ9DEkUjUMHG/view?usp=drive_link",
               path="data/mndow_lakes_sf_allDataUntransformed.rds")
drive_download("https://drive.google.com/file/d/1qzRnAhLthryc9Bj2x6bvTkOHNcs_Xu1E/view?usp=drive_link",
               path="data/QBO_NOAA.txt")
drive_download("https://drive.google.com/file/d/1p2QebafmkckZXC9-j2WUYDvcLDU5sr13/view?usp=drive_link",
               path="data/PDO_NOAA.txt")
drive_download("https://drive.google.com/file/d/1b9kvLD0H3Q2Rq_WXfWbNK0dN5D1FmKEN/view?usp=drive_link",
               path="data/SOLARCYCLE_NOAA.txt")
drive_download("https://drive.google.com/file/d/1MbYazwIzYclZJRCnkBq5wDCfepiufc40/view?usp=drive_link",
               path="data/NAO_DJFMwinter.txt")
drive_download("https://drive.google.com/file/d/1PeLKss2zbfKSg_ZlLi91IFV33C3srMdN/view?usp=drive_link",
               path="data/BESTlong_ENSO_NOAA.txt")
drive_download("https://drive.google.com/file/d/1sQ4zHuMTivxGsNWEwvXgKrIk7qMTa7O8/view?usp=drive_link",
               path="data/BEST_ENSO_NOAA.txt")

## Loading data into R project ----

### Ice cover data ----

ice <- read.csv("data/ice_duration_summarized.csv", 
                stringsAsFactors = F) %>%
  mutate(min_ice_on_date = as.Date(min_ice_on_date),
         max_ice_on_date = as.Date(max_ice_on_date),
         min_ice_off_date = as.Date(min_ice_off_date),
         max_ice_off_date = as.Date(max_ice_off_date),
         DOW = fixlakeid(DOW),
         fall.year = winter.year-1,
         fwinter.year=factor(winter.year),
         max_ice_on_julian2 = ifelse(max_ice_on_julian < 200, 365 + max_ice_on_julian, max_ice_on_julian), min_ice_on_julian2 = ifelse(min_ice_on_julian < 200, 365 + min_ice_on_julian, min_ice_on_julian)) %>%
  mutate(ID = factor(paste("dow", DOW, sep = "_"))) %>%
  dplyr::select(-X)

rk_ice <- read.csv("data/Rainy&Kabetogama_ice_duration_summarized.csv", 
                   stringsAsFactors = F) %>%
  mutate(min_ice_on_date = as.Date(min_ice_on_date),
         max_ice_on_date = as.Date(max_ice_on_date),
         min_ice_off_date = as.Date(min_ice_off_date),
         max_ice_off_date = as.Date(max_ice_off_date),
         DOW = fixlakeid(DOW),
         fall.year = winter.year-1,
         fwinter.year=factor(winter.year),
         max_ice_on_julian2 = ifelse(max_ice_on_julian < 200, 365 + max_ice_on_julian, max_ice_on_julian), min_ice_on_julian2 = ifelse(min_ice_on_julian < 200, 365 + min_ice_on_julian, min_ice_on_julian)) %>%
  mutate(ID = factor(paste("dow", DOW, sep = "_"))) %>%
  dplyr::select(-1)
head(rk_ice)

icein <- read.csv("data/ice_on_summarized.csv", 
                  stringsAsFactors = F) %>%
  mutate(min_ice_on_date = as.Date(min_ice_on_date),
         max_ice_on_date = as.Date(max_ice_on_date),
         DOW = fixlakeid(DOW),
         fall.year = winter.year-1,
         fwinter.year=factor(winter.year),
         max_ice_on_julian2 = ifelse(max_ice_on_julian < 200, 365 + max_ice_on_julian, max_ice_on_julian), min_ice_on_julian2 = ifelse(min_ice_on_julian < 200, 365 + min_ice_on_julian, min_ice_on_julian)) %>%
  mutate(ID = factor(paste("dow", DOW, sep = "_"))) %>%
  dplyr::select(-X)

iceout <- read.csv("data/ice_off_summarized.csv", 
                   stringsAsFactors = F) %>%
  mutate(min_ice_off_date = as.Date(min_ice_off_date),
         max_ice_off_date = as.Date(max_ice_off_date),
         DOW = fixlakeid(DOW),
         winter.year=year) %>%
  mutate(ID = factor(paste("dow", DOW, sep = "_")),
         fall.year = winter.year-1,
         fwinter.year=factor(winter.year)) %>%
  dplyr::select(-X)

rk_iceout <- read.csv("data/Rainy&Kabetogama_ice_off_summarized.csv", 
                      stringsAsFactors = F) %>%
  mutate(min_ice_off_date = as.Date(min_ice_off_date),
         max_ice_off_date = as.Date(max_ice_off_date),
         DOW = fixlakeid(DOW),
         winter.year=year) %>%
  mutate(ID = factor(paste("dow", DOW, sep = "_")),
         fall.year = winter.year-1,
         fwinter.year=factor(winter.year)) %>%
  dplyr::select(-1)

# Add rk_ice and rk_iceout to ice and iceout objects
ice <- rbind(ice, rk_ice)
iceout <- rbind(iceout, rk_iceout)

# Test for and remove duplicated Rainy & Kabetogama data
ice %>% filter(duplicated(ice)) # no duplicates
iceout %>% filter(duplicated(iceout)) # duplicates

iceout <- unique(iceout)

# There are still multiple observations for some iceout lake years in Rainy & Kabetogama
ice %>% group_by(winter.year, ID) %>% count() %>% filter(n> 1) # good!
iceout %>% group_by(winter.year, ID) %>% count() %>% filter(n> 1)

# select minimum duration and ice off values and maximum ice on values

iceout <- iceout %>%
  group_by(DOW, year, winter.year, ID, fall.year, fwinter.year) %>%
  summarize(min_ice_off_julian=min(min_ice_off_julian),
            max_ice_off_julian=max(max_ice_off_julian),
            min_ice_off_date=min(min_ice_off_date),
            max_ice_off_date=max(max_ice_off_date),
            N_ice_off=sum(N_ice_off),
            range_ice_off=max_ice_off_julian-min_ice_off_julian) %>%
  ungroup() %>%
  as.data.frame()

# check for multiple observations
iceout %>% group_by(winter.year, ID) %>% count() %>% filter(n> 1) # good!

# Number lakes for each ice metric
length(unique(ice$DOW)) # 1244 lakes
length(unique(icein$DOW)) # 1478 lakes
length(unique(iceout$DOW)) # 1871 lakes

### Lake size ----
MN.shape <-  readRDS("data/mndow_lakes_sf_allDataUntransformed.rds")
st_crs(MN.shape) <- 26915


# Aggregate polygons by DOW
MN.shape_agg <- MN.shape %>% group_by(dowlknum) %>% summarise()
#rm(MN.shape)

# Measure lake area
MN.shape_agg$AggDOW_Area_m2 <- st_area(MN.shape_agg)
MN.shape_agg$AggDOW_Area_acres <- units::set_units(MN.shape_agg$AggDOW_Area_m2, acre)

MN.shape_agg <- MN.shape_agg %>% mutate(DOW = fixlakeid(dowlknum)) %>%
  mutate(ID = factor(paste("dow", DOW, sep = "_")))

### Attach lake area to ice cover data----
# Ice duration
ice_spatial <- inner_join(ice, MN.shape_agg) %>% rename(Area_acres=AggDOW_Area_acres) %>% dplyr::select(-dowlknum)
ice_spatial$ID <- as.factor(ice_spatial$ID)
ice_spatial$Area_acres <- as.numeric(ice_spatial$Area_acres)
ice_spatial$AggDOW_Area_m2 <- as.numeric(ice_spatial$AggDOW_Area_m2)
ice_spatial$lnArea_acres <- log(ice_spatial$Area_acres)

# Iceout data
iceout_spatial <- inner_join(iceout, MN.shape_agg) %>% rename(Area_acres=AggDOW_Area_acres) %>%  dplyr::select(-dowlknum)
str(iceout_spatial)
iceout_spatial$ID <- as.factor(iceout_spatial$ID)
iceout_spatial$Area_acres <- as.numeric(iceout_spatial$Area_acres)
iceout_spatial$AggDOW_Area_m2 <- as.numeric(iceout_spatial$AggDOW_Area_m2)
iceout_spatial$lnArea_acres <- log(iceout_spatial$Area_acres)


### Political boundaries----
states <- st_as_sf(maps::map("state", plot = FALSE, fill = TRUE))
states <- st_transform(states, crs = st_crs(MN.shape_agg))

MN_outline <- st_as_sf(states %>% filter(ID == "minnesota"))

# MN Counties
MN_counties <- st_as_sf(maps::map("county", regions = "minnesota", plot = FALSE, fill = TRUE))
MN_counties <- st_transform(MN_counties, crs = st_crs(MN.shape_agg))

MN_arrowheadcounties <- MN_counties %>% filter(ID%in%c("minnesota,carlton", "minnesota,cook", "minnesota,lake", "minnesota,st louis"))
MN_arrowheadcounties_outline <- st_union(MN_arrowheadcounties)



### Climate indices ----

# Quasi-biennial oscillations
qbo <- read.table("data/QBO_NOAA.txt")
colnames(qbo) <- c("year", paste("month", 1:12, sep = ""))

# Pacific decadal oscillations
pdo <- read.table("data/PDO_NOAA.txt")
colnames(pdo) <- c("year", paste("month", 1:12, sep = ""))

# Solar cycles (important in Sharma/Magnuson models)
sun <- read.table("data/SOLARCYCLE_NOAA.txt")
colnames(sun) <- c("year", paste("month", 1:12, sep = ""))

# North Atlantic oscillations (important in Sharma/Magnuson models) - now using winter version (Dec, Jan, Feb, March)
nao <- read.table("data/NAO_DJFMwinter.txt") 
colnames(nao) <- c("year", "index")

# El Nino-Southern oscillations (important in Sharma/Magnuson models)
best_enso_old <- read.table("data/BESTlong_ENSO_NOAA.txt")
colnames(best_enso_old) <- c("year", paste("month", 1:12, sep = ""))

best_enso_new <- read.table("data/BEST_ENSO_NOAA.txt")
colnames(best_enso_new) <- c("year", paste("month", 1:12, sep = ""))

enso <- rbind(best_enso_old %>% filter(year < 1948), best_enso_new)

#### Summarize and detrend monthly indices ----
# annual (Jan-Dec) averages (e.g., QBO) and 
# July (fall year) -June (winter year) averages (e.g., QBOjj)
# ie, July is prior to formation and June is after breakup

# detrend indices (e.g., QBO_dt or QBOjj_dt)

# QBO
QBO <- qbo %>% 
  filter(year <= 2022) %>%
  pivot_longer(month1:month12, names_to="month_lab", values_to="QBO") %>%
  mutate(winter.year=year) %>%
  group_by(winter.year) %>%
  summarize(QBO=mean(QBO, na.rm=T)) %>%
  mutate(QBO_dt=detrend(QBO, tt='linear'))

QBOjj <- qbo %>% 
  filter(year <= 2022) %>%
  pivot_longer(month1:month12, names_to="month_lab", values_to="QBO") %>%
  separate_wider_delim(month_lab, "month", names=c("text", "month")) %>%
  dplyr::select(-text) %>%
  mutate(month=as.numeric(month),
         fall.year=year-1,
         winter.year=ifelse(month<7, fall.year, year)) %>%
  group_by(winter.year) %>%
  summarize(QBOjj=mean(QBO, na.rm=T)) %>%
  mutate(QBOjj_dt=detrend(QBOjj, tt='linear'))

# view QBO
ggplot() + 
  geom_line(aes(x=QBO$winter.year, y=QBO$QBO)) + 
  geom_line(aes(x=QBO$winter.year, y=QBO$QBO_dt), lty=2) + 
  geom_line(aes(x=QBOjj$winter.year, y=QBOjj$QBOjj), col='red') + 
  geom_line(aes(x=QBOjj$winter.year, y=QBOjj$QBOjj_dt), lty=2, col='red') + 
  theme_classic(12)

# PDO
PDO <- pdo %>% 
  filter(year <= 2022) %>%
  pivot_longer(month1:month12, names_to="month_lab", values_to="PDO") %>%
  mutate(winter.year=year) %>%
  group_by(winter.year) %>%
  summarize(PDO=mean(PDO, na.rm=T)) %>%
  mutate(PDO_dt=detrend(PDO, tt='linear'))

PDOjj <- pdo %>% 
  filter(year <= 2022) %>%
  pivot_longer(month1:month12, names_to="month_lab", values_to="PDO") %>%
  separate_wider_delim(month_lab, "month", names=c("text", "month")) %>%
  dplyr::select(-text) %>%
  mutate(month=as.numeric(month),
         fall.year=year-1,
         winter.year=ifelse(month<7, fall.year, year)) %>%
  group_by(winter.year) %>%
  summarize(PDOjj=mean(PDO, na.rm=T)) %>%
  mutate(PDOjj_dt=detrend(PDOjj, tt='linear'))

# view PDO
ggplot() + 
  geom_line(aes(x=PDO$winter.year, y=PDO$PDO)) + 
  geom_line(aes(x=PDO$winter.year, y=PDO$PDO_dt), lty=2) + 
  geom_line(aes(x=PDOjj$winter.year, y=PDOjj$PDOjj), col='red') + 
  geom_line(aes(x=PDOjj$winter.year, y=PDOjj$PDOjj_dt), lty=2, col='red') + 
  theme_classic(12)

# SUN
SUN <- sun %>% 
  filter(year <= 2022) %>%
  pivot_longer(month1:month12, names_to="month_lab", values_to="SUN") %>%
  mutate(winter.year=year) %>%
  group_by(winter.year) %>%
  summarize(SUN=mean(SUN, na.rm=T)) %>%
  mutate(SUN_dt=detrend(SUN, tt='linear'))

SUNjj <- sun %>% 
  filter(year <= 2022) %>%
  pivot_longer(month1:month12, names_to="month_lab", values_to="SUN") %>%
  separate_wider_delim(month_lab, "month", names=c("text", "month")) %>%
  dplyr::select(-text) %>%
  mutate(month=as.numeric(month),
         fall.year=year-1,
         winter.year=ifelse(month<7, fall.year, year)) %>%
  group_by(winter.year) %>%
  summarize(SUNjj=mean(SUN, na.rm=T)) %>%
  mutate(SUNjj_dt=detrend(SUNjj, tt='linear'))

# view SUN
ggplot() + 
  geom_line(aes(x=SUN$winter.year, y=SUN$SUN)) + 
  geom_line(aes(x=SUN$winter.year, y=SUN$SUN_dt), lty=2) + 
  geom_line(aes(x=SUNjj$winter.year, y=SUNjj$SUNjj), col='red') + 
  geom_line(aes(x=SUNjj$winter.year, y=SUNjj$SUNjj_dt), lty=2, col='red') + 
  theme_classic(12)


# NAO
NAO <- data.frame(winter.year=nao$year,
                  NAO=nao$index,
                  NAO_dt=detrend(nao$index, tt='linear'))

# view NAO
ggplot() + 
  geom_line(aes(x=NAO$winter.year, y=NAO$NAO)) + 
  geom_line(aes(x=NAO$winter.year, y=NAO$NAO_dt), lty=2) + 
  theme_classic(12)

# ENSO
ENSO <- enso %>% 
  filter(year <= 2022) %>%
  pivot_longer(month1:month12, names_to="month_lab", values_to="ENSO") %>%
  mutate(winter.year=year) %>%
  group_by(winter.year) %>%
  summarize(ENSO=mean(ENSO, na.rm=T)) %>%
  mutate(ENSO_dt=detrend(ENSO, tt='linear'))

ENSOjj <- enso %>% 
  filter(year <= 2022) %>%
  pivot_longer(month1:month12, names_to="month_lab", values_to="ENSO") %>%
  separate_wider_delim(month_lab, "month", names=c("text", "month")) %>%
  dplyr::select(-text) %>%
  mutate(month=as.numeric(month),
         fall.year=year-1,
         winter.year=ifelse(month<7, fall.year, year)) %>%
  group_by(winter.year) %>%
  summarize(ENSOjj=mean(ENSO, na.rm=T)) %>%
  mutate(ENSOjj_dt=detrend(ENSOjj, tt='linear'))

# view ENSO
ggplot() + 
  geom_line(aes(x=ENSO$winter.year, y=ENSO$ENSO)) + 
  geom_line(aes(x=ENSO$winter.year, y=ENSO$ENSO_dt), lty=2) + 
  geom_line(aes(x=ENSOjj$winter.year, y=ENSOjj$ENSOjj), col='red') + 
  geom_line(aes(x=ENSOjj$winter.year, y=ENSOjj$ENSOjj_dt), lty=2, col='red') + 
  theme_classic(12)


### Join climate indices to ice cover data----

ice_spatial <- left_join(ice_spatial, QBO, by='winter.year') %>%
  left_join(., QBOjj, by='winter.year') %>%
  left_join(., PDO, by='winter.year') %>%
  left_join(., PDOjj, by='winter.year') %>%
  left_join(., SUN, by='winter.year') %>%
  left_join(., SUNjj, by='winter.year') %>%
  left_join(., NAO, by='winter.year') %>%
  left_join(., ENSO, by='winter.year') %>%
  left_join(., ENSOjj, by='winter.year')

iceout_spatial <- left_join(iceout_spatial, QBO, by='winter.year') %>%
  left_join(., QBOjj, by='winter.year') %>%
  left_join(., PDO, by='winter.year') %>%
  left_join(., PDOjj, by='winter.year') %>%
  left_join(., SUN, by='winter.year') %>%
  left_join(., SUNjj, by='winter.year') %>%
  left_join(., NAO, by='winter.year') %>%
  left_join(., ENSO, by='winter.year') %>%
  left_join(., ENSOjj, by='winter.year')


### Add x and y coordinates to ice cover data----

ice_spatial$x <- st_coordinates(st_centroid(st_as_sf(ice_spatial)))[,"X"]
ice_spatial$y <- st_coordinates(st_centroid(st_as_sf(ice_spatial)))[,"Y"]

iceout_spatial$x <- st_coordinates(st_centroid(st_as_sf(iceout_spatial)))[,"X"]
iceout_spatial$y <- st_coordinates(st_centroid(st_as_sf(iceout_spatial)))[,"Y"]

## Summarize number of ice cover observations for each lake in different time periods----

### All three metrics----
# All duration data
ice_lakesum <- ice_spatial  %>%  
  group_by(ID) %>%  
  summarise(n_years=n_distinct(winter.year),  n = n(), 
            Number50s = n_distinct(winter.year[winter.year > 1949 & winter.year <1960]), # 1950-1959
            Number60s = n_distinct(winter.year[winter.year > 1959 & winter.year <1970]), # 1960-1969
            Number70s = n_distinct(winter.year[winter.year > 1969 & winter.year <1980]), # 1970-1979
            Number80s = n_distinct(winter.year[winter.year > 1979 & winter.year <1990]), # 1980-1989
            Number90s = n_distinct(winter.year[winter.year > 1989 & winter.year <2000]), # 1990-1999
            Number00s = n_distinct(winter.year[winter.year > 1999 & winter.year <2010]), # 2000-2009
            Number10s = n_distinct(winter.year[winter.year > 2009 & winter.year <2020]), # 2010-2019
            Number20s = n_distinct(winter.year[winter.year > 2019 & winter.year <2030])) 

# Duration data since 1949
ice_lakesum1949 <- ice_spatial %>% filter(winter.year > 1948) %>%
  group_by(ID) %>%  
  summarise(n_years=n_distinct(winter.year), n = n(), 
            Number50s = n_distinct(winter.year[winter.year > 1949 & winter.year <1960]), # 1950-1959
            Number60s = n_distinct(winter.year[winter.year > 1959 & winter.year <1970]), # 1960-1969
            Number70s = n_distinct(winter.year[winter.year > 1969 & winter.year <1980]), # 1970-1979
            Number80s = n_distinct(winter.year[winter.year > 1979 & winter.year <1990]), # 1980-1989
            Number90s = n_distinct(winter.year[winter.year > 1989 & winter.year <2000]), # 1990-1999
            Number00s = n_distinct(winter.year[winter.year > 1999 & winter.year <2010]), # 2000-2009
            Number10s = n_distinct(winter.year[winter.year > 2009 & winter.year <2020]), # 2010-2019
            Number20s = n_distinct(winter.year[winter.year > 2019 & winter.year <2030])) 

# Duration data since 1970
ice_lakesum1970  <- ice_spatial %>% filter(winter.year > 1969) %>%
  group_by(ID) %>%  
  summarise(n_years=n_distinct(winter.year), n = n(), 
            Number70s = n_distinct(winter.year[winter.year > 1969 & winter.year <1980]), # 1970-1979
            Number80s = n_distinct(winter.year[winter.year > 1979 & winter.year <1990]), # 1980-1989
            Number90s = n_distinct(winter.year[winter.year > 1989 & winter.year <2000]), # 1990-1999
            Number00s = n_distinct(winter.year[winter.year > 1999 & winter.year <2010]), # 2000-2009
            Number10s = n_distinct(winter.year[winter.year > 2009 & winter.year <2020]), # 2010-2019
            Number20s = n_distinct(winter.year[winter.year > 2019 & winter.year <2030])) 

### Longer ice out data----
iceout_lakesum <- iceout_spatial  %>%  
  group_by(ID) %>%  
  summarise(n_years=n_distinct(winter.year), n = n(), 
            Number50s = n_distinct(winter.year[winter.year > 1949 & winter.year <1960]), # 1950-1959
            Number60s = n_distinct(winter.year[winter.year > 1959 & winter.year <1970]), # 1960-1969
            Number70s = n_distinct(winter.year[winter.year > 1969 & winter.year <1980]), # 1970-1979
            Number80s = n_distinct(winter.year[winter.year > 1979 & winter.year <1990]), # 1980-1989
            Number90s = n_distinct(winter.year[winter.year > 1989 & winter.year <2000]), # 1990-1999
            Number00s = n_distinct(winter.year[winter.year > 1999 & winter.year <2010]), # 2000-2009
            Number10s = n_distinct(winter.year[winter.year > 2009 & winter.year <2020]), # 2010-2019
            Number20s = n_distinct(winter.year[winter.year > 2019 & winter.year <2030])) 

iceout_lakesum1949 <- iceout_spatial %>% filter(winter.year > 1948) %>%
  group_by(ID) %>%  
  summarise(n_years=n_distinct(winter.year), n = n(), 
            Number50s = n_distinct(winter.year[winter.year > 1949 & winter.year <1960]), # 1950-1959
            Number60s = n_distinct(winter.year[winter.year > 1959 & winter.year <1970]), # 1960-1969
            Number70s = n_distinct(winter.year[winter.year > 1969 & winter.year <1980]), # 1970-1979
            Number80s = n_distinct(winter.year[winter.year > 1979 & winter.year <1990]), # 1980-1989
            Number90s = n_distinct(winter.year[winter.year > 1989 & winter.year <2000]), # 1990-1999
            Number00s = n_distinct(winter.year[winter.year > 1999 & winter.year <2010]), # 2000-2009
            Number10s = n_distinct(winter.year[winter.year > 2009 & winter.year <2020]), # 2010-2019
            Number20s = n_distinct(winter.year[winter.year > 2019 & winter.year <2030])) 

iceout_lakesum1970  <- iceout_spatial %>% filter(winter.year > 1969) %>%
  group_by(ID) %>%  
  summarise(n_years=n_distinct(winter.year), n = n(), 
            Number70s = n_distinct(winter.year[winter.year > 1969 & winter.year <1980]), # 1970-1979
            Number80s = n_distinct(winter.year[winter.year > 1979 & winter.year <1990]), # 1980-1989
            Number90s = n_distinct(winter.year[winter.year > 1989 & winter.year <2000]), # 1990-1999
            Number00s = n_distinct(winter.year[winter.year > 1999 & winter.year <2010]), # 2000-2009
            Number10s = n_distinct(winter.year[winter.year > 2009 & winter.year <2020]), # 2010-2019
            Number20s = n_distinct(winter.year[winter.year > 2019 & winter.year <2030])) 



## Subset lakes with 10 and 30 years of ice cover data since 1949----
### All three metrics----
lakes_w10y <- ice_lakesum1949$ID[ice_lakesum1949$n_years >= 10]
ice_w10y <- ice_spatial %>% filter(ID %in% lakes_w10y, winter.year > 1948) 

lakes_w30y <- ice_lakesum1949$ID[ice_lakesum1949$n_years >= 30]
ice_w30y <- ice_spatial %>% filter(ID %in% lakes_w30y, winter.year > 1948) 

# Pull out data for lakes with at least 10 yrs since 1949 **** ONLY COMPLETE CASES used for modeling****
ice_w10y_condense <- ice_w10y %>% dplyr::select(min_duration, max_duration, winter.year, fwinter.year, ID, x, y, ENSO, ENSOjj, NAO, QBO, QBOjj, SUN, SUNjj, PDO, PDOjj, lnArea_acres, Area_acres, min_ice_on_julian, max_ice_on_julian, max_ice_on_julian2, min_ice_off_julian, max_ice_off_julian, min_duration, max_duration)
ice_w10y_condense <- ice_w10y_condense[complete.cases(ice_w10y_condense),]

ice_w30y_condense <- ice_w30y %>% dplyr::select(min_duration, max_duration, winter.year, fwinter.year, ID, x, y, ENSO, ENSOjj, NAO, QBO, QBOjj, SUN, SUNjj, PDO, PDOjj, lnArea_acres, Area_acres, min_ice_on_julian, max_ice_on_julian, max_ice_on_julian2, min_ice_off_julian, max_ice_off_julian, min_duration, max_duration)
ice_w30y_condense <- ice_w30y_condense[complete.cases(ice_w30y_condense),]

length(unique(ice_w10y_condense$ID)) # 262 lakes with 10 yrs ice duration data since 1949
dim(ice_w10y_condense) # 5,601 observations
range(ice_w10y_condense$winter.year) # 1949-2022

length(unique(ice_w30y_condense$ID)) # 54 lakes with 50 yrs ice duration data since 1949
dim(ice_w30y_condense) # 2247 observations
range(ice_w30y_condense$winter.year) # 1949-2022


ice_all_condense <- ice_spatial %>% dplyr::select(min_duration, max_duration, winter.year, fwinter.year, ID, x, y, ENSO, ENSOjj, NAO, QBO, QBOjj, SUN, SUNjj, PDO, PDOjj, lnArea_acres, Area_acres, min_ice_on_julian, max_ice_on_julian, max_ice_on_julian2, min_ice_off_julian, max_ice_off_julian, min_duration, max_duration)
ice_all_condense <- ice_all_condense[complete.cases(ice_all_condense),]


### Longer ice out data----
lakesout_w10y <- iceout_lakesum1949$ID[iceout_lakesum1949$n_years >= 10]
iceout_w10y <- iceout_spatial %>% filter(ID %in% lakesout_w10y, winter.year > 1948) 

lakesout_w30y <- iceout_lakesum1949$ID[iceout_lakesum1949$n_years >= 30]
iceout_w30y <- iceout_spatial %>% filter(ID %in% lakesout_w30y, winter.year > 1948) 

# Pull out data for lakes with at least 10 yrs since 1949 - **** ONLY COMPLETE CASES used for modeling****
iceout_w10y_condense <- iceout_w10y %>% dplyr::select(winter.year, fwinter.year, ID, x, y, ENSO, ENSOjj, NAO, QBO, QBOjj, SUN, SUNjj, PDO, PDOjj, lnArea_acres, Area_acres, min_ice_off_julian, max_ice_off_julian)
iceout_w10y_condense <- iceout_w10y_condense[complete.cases(iceout_w10y_condense),]

iceout_w30y_condense <- iceout_w30y %>% dplyr::select(winter.year, fwinter.year, ID, x, y, ENSO, ENSOjj, NAO, QBO, QBOjj, SUN, SUNjj, PDO, PDOjj, lnArea_acres, Area_acres, min_ice_off_julian, max_ice_off_julian)
iceout_w30y_condense <- iceout_w30y_condense[complete.cases(iceout_w30y_condense),]

length(unique(iceout_w10y_condense$ID))  # 586 lakes with 10 yrs ice off since 1949
dim(iceout_w10y_condense) # 14428
range(iceout_w10y_condense$winter.year) # 1949-2022

length(unique(iceout_w30y_condense$ID))  # 163 lakes with 50 yrs ice off since 1949
dim(iceout_w30y_condense) # 7800 observations
range(iceout_w30y_condense$winter.year) # 1949-2022

iceout_all_condense <- iceout_spatial %>% dplyr::select(winter.year, fwinter.year, ID, x, y, ENSO, ENSOjj, NAO, QBO, QBOjj, SUN, SUNjj, PDO, PDOjj, lnArea_acres, Area_acres, min_ice_off_julian, max_ice_off_julian)
iceout_all_condense <- iceout_all_condense[complete.cases(iceout_all_condense),]


## Subset lakes with 10 years of data, at least one observation in 70s, 80s, 90s, 00s, 10s----

IDseachDecade <- ice_lakesum1970 %>% filter(Number70s > 0, Number80s > 0, Number90s > 0, Number00s > 0, Number10s > 0, n_years >= 10) 
ice_w10y1970 <- ice_spatial %>% filter(ID %in% IDseachDecade$ID, winter.year > 1969) 

ice_w10y1970_condense <- ice_w10y1970 %>% dplyr::select(min_duration, max_duration, winter.year, fwinter.year, ID, x, y, ENSO, ENSOjj, NAO, QBO, QBOjj, SUN, SUNjj, PDO, PDOjj, lnArea_acres, Area_acres, min_ice_on_julian, max_ice_on_julian, max_ice_on_julian2, min_ice_off_julian, max_ice_off_julian, min_duration, max_duration)
ice_w10y1970_condense <- ice_w10y1970_condense[complete.cases(ice_w10y1970_condense),]
ice_w10y1970_condense$ID <- droplevels(ice_w10y1970_condense$ID )


dim(ice_w10y1970_condense) # 1531 observations
length(unique(ice_w10y1970_condense$ID)) # 41 lakes

## Plot median lake area and sample size over time -- lakes w/ 10 yrs since 1949----
size_sum <- ice_w10y_condense %>% group_by(winter.year) %>% summarise(avg_size=mean(Area_acres), avg_logsize=log(mean(Area_acres)), med_size=median(Area_acres), med_logsize=log(median(Area_acres)), min_size=min(Area_acres), max_size=max(Area_acres), min_logsize=log(min(Area_acres)), max_logsize=log(max(Area_acres)), n=n())

# Plot size over time
plot(med_size~winter.year, data=size_sum, type='l', xlim=c(1949, 2019), ylim=c(0, 3000), las=1, ylab="Area (acres)", xlab="Year", cex.lab=1.2, cex.axis=1.2)

# Plot number lakes over time: NumLakes_vsYear_1950.png
plot(n~winter.year, data=size_sum, type='l', xlim=c(1949, 2019), las=1, ylab="Number of lakes", xlab="Year", ylim=c(0, 200))

