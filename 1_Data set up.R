# Variable phenology but consistent loss of ice cover of 1,213 Minnesota lakes
# Walsh et al. submitted to L&O Letters
# Data set up

# Set up repository folders ----

# The data folder has some of the climate indices and the lake shapefile
# ifelse(!dir.exists("data"), 
#        dir.create("data"), 
#        "Folder exists already")
ifelse(!dir.exists("figures"), 
       dir.create("figures"), 
       "Folder exists already")
ifelse(!dir.exists("tables"), 
       dir.create("tables"), 
       "Folder exists already")
ifelse(!dir.exists("outputdata"), 
       dir.create("outputdata"), 
       "Folder exists already")

# Data set up ----

## Packages ----
# Data manipulation
library(tidyverse)
library(viridis)
library(lubridate)
library(data.table)
library(RCurl)

# Plotting
library(gridExtra)
library(grid)
library(broom)
library(cowplot)
library(scales)

# Function for detrending climate indices
library(pracma)

# GAM
library(mgcv)
library(gratia)
library(itsadug)
library(parallel)

# Spatial
library(sf)
library(maps)

# The "mnsentinellakes" package is downloaded through github using the package "remotes":
install.packages("remotes")
remotes::install_github("mnsentinellakes/mnsentinellakes")

library(mnsentinellakes)


## Download data from DRUM ----
# DRUM Repository:
# https://conservancy.umn.edu/items/940c6197-486b-437c-96ca-cca9b4534dfa

# Download the file and save it to the data folder
download.file("https://conservancy.umn.edu/bitstreams/f66ca8ac-dedc-46c0-bf29-e136b0723a4a/download",
               destfile="data/ice_duration.csv")
download.file("https://conservancy.umn.edu/bitstreams/0a71aa4b-071a-4685-a84c-ea5e0fe0691c/download",
               destfile="data/ice_in.csv")
download.file("https://conservancy.umn.edu/bitstreams/32150e3a-62be-4f95-911f-9afbd9011fa4/download",
               destfile="data/ice_out.csv")


## Loading data into R project ----

### Ice cover data ----

# Read ice files into R, classify relevant columns
ice <- read.csv("data/ice_duration.csv", 
                stringsAsFactors = F) %>%
  mutate(min_ice_on_date = as.Date(min_ice_on_date),
         max_ice_on_date = as.Date(max_ice_on_date),
         min_ice_off_date = as.Date(min_ice_off_date),
         max_ice_off_date = as.Date(max_ice_off_date),
         DOW = fixlakeid(DOW),
         fwinter.year=factor(winter.year),
         ID=factor(ID))

icein <- read.csv("data/ice_in.csv", 
                  stringsAsFactors = F) %>%
  mutate(min_ice_on_date = as.Date(min_ice_on_date),
         max_ice_on_date = as.Date(max_ice_on_date),
         DOW = fixlakeid(DOW),
         fwinter.year=factor(winter.year)) %>%
  mutate(ID = factor(paste("dow", DOW, sep = "_")))

iceout <- read.csv("data/ice_out.csv", 
                   stringsAsFactors = F) %>%
  mutate(min_ice_off_date = as.Date(min_ice_off_date),
         max_ice_off_date = as.Date(max_ice_off_date),
         DOW = fixlakeid(DOW),
         fwinter.year=factor(winter.year),
         ID=factor(ID))

# Number lakes for each ice metric
length(unique(ice$DOW)) # 1244 lakes
length(unique(icein$DOW)) # 1479 lakes
length(unique(iceout$DOW)) # 1871 lakes

### Lake size ----

# Shapefile contains lake polygons
MN.shape <-  readRDS("data/mndow_lakes_sf_allDataUntransformed.rds")
st_crs(MN.shape) <- 26915


# Aggregate polygons by DOW
MN.shape_agg <- MN.shape %>% group_by(dowlknum) %>% summarise()

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
iceout_spatial$ID <- as.factor(iceout_spatial$ID)
iceout_spatial$Area_acres <- as.numeric(iceout_spatial$Area_acres)
iceout_spatial$AggDOW_Area_m2 <- as.numeric(iceout_spatial$AggDOW_Area_m2)
iceout_spatial$lnArea_acres <- log(iceout_spatial$Area_acres)


### Political boundaries----
states <- st_as_sf(maps::map("state", plot = FALSE, fill = TRUE))
states <- st_transform(states, crs = st_crs(MN.shape_agg))

# MN
MN_outline <- st_as_sf(states %>% filter(ID == "minnesota"))

# MN Counties
MN_counties <- st_as_sf(maps::map("county", regions = "minnesota", plot = FALSE, fill = TRUE))
MN_counties <- st_transform(MN_counties, crs = st_crs(MN.shape_agg))

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


### Add climate indices to ice cover data----
ice_spatial$ENSO <- NA
ice_spatial$NAO <- NA
ice_spatial$QBO <- NA
ice_spatial$SUN <- NA
ice_spatial$PDO <- NA
ice_spatial$ENSOjj <- NA
ice_spatial$QBOjj <- NA
ice_spatial$SUNjj <- NA
ice_spatial$PDOjj <- NA

# Match winter.year with mean annual climate index values (climate index values are constant across lake)
for(i in 1:nrow(ice_spatial)){
  # Year out averages
  if(ice_spatial$winter.year[i] %in% enso$year) ice_spatial$ENSO[i] <- rowMeans(enso[enso$year == ice_spatial$winter.year[i], 2:13])
  if(ice_spatial$winter.year[i] %in% qbo$year) ice_spatial$QBO[i] <- rowMeans(qbo[qbo$year == ice_spatial$winter.year[i], 2:13])
  if(ice_spatial$winter.year[i] %in% sun$year) ice_spatial$SUN[i] <- rowMeans(sun[sun$year == ice_spatial$winter.year[i], 2:13])
  if(ice_spatial$winter.year[i] %in% pdo$year) ice_spatial$PDO[i] <- rowMeans(pdo[pdo$year == ice_spatial$winter.year[i], 2:13])
  if(ice_spatial$winter.year[i] %in% nao$year) ice_spatial$NAO[i] <- nao$index[nao$year == ice_spatial$winter.year[i]]
  
  # June-June averages: winter.year=months 1-6 (cols 2-7), fall.year=months 7-12 (cols 8-13)
  if(ice_spatial$winter.year[i] %in% enso$year & ice_spatial$fall.year[i] %in% enso$year) ice_spatial$ENSOjj[i] <- rowMeans(cbind(enso[enso$year == ice_spatial$fall.year[i], 8:13], enso[enso$year == ice_spatial$winter.year[i], 2:7]))
  if(ice_spatial$winter.year[i] %in% qbo$year & ice_spatial$fall.year[i] %in% qbo$year) ice_spatial$QBOjj[i] <- rowMeans(cbind(qbo[qbo$year == ice_spatial$fall.year[i], 8:13], qbo[qbo$year == ice_spatial$winter.year[i], 2:7]))
  if(ice_spatial$winter.year[i] %in% sun$year & ice_spatial$fall.year[i] %in% sun$year) ice_spatial$SUNjj[i] <- rowMeans(cbind(sun[sun$year == ice_spatial$fall.year[i], 8:13], sun[sun$year == ice_spatial$winter.year[i], 2:7]))
  if(ice_spatial$winter.year[i] %in% pdo$year & ice_spatial$fall.year[i] %in% pdo$year) ice_spatial$PDOjj[i] <- rowMeans(cbind(pdo[pdo$year == ice_spatial$fall.year[i], 8:13], pdo[pdo$year == ice_spatial$winter.year[i], 2:7]))
}

# Add climate indices to ice out
iceout_spatial$ENSO <- NA
iceout_spatial$NAO <- NA
iceout_spatial$QBO <- NA
iceout_spatial$SUN <- NA
iceout_spatial$PDO <- NA
iceout_spatial$ENSOjj <- NA
iceout_spatial$QBOjj <- NA
iceout_spatial$SUNjj <- NA
iceout_spatial$PDOjj <- NA

# Match winter.year with mean annual climate index values (climate index values are constant across lake)
for(i in 1:nrow(iceout_spatial)){
  # Year out averages
  if(iceout_spatial$winter.year[i] %in% enso$year) iceout_spatial$ENSO[i] <- rowMeans(enso[enso$year == iceout_spatial$winter.year[i], 2:13])
  if(iceout_spatial$winter.year[i] %in% qbo$year) iceout_spatial$QBO[i] <- rowMeans(qbo[qbo$year == iceout_spatial$winter.year[i], 2:13])
  if(iceout_spatial$winter.year[i] %in% sun$year) iceout_spatial$SUN[i] <- rowMeans(sun[sun$year == iceout_spatial$winter.year[i], 2:13])
  if(iceout_spatial$winter.year[i] %in% pdo$year) iceout_spatial$PDO[i] <- rowMeans(pdo[pdo$year == iceout_spatial$winter.year[i], 2:13])
  if(iceout_spatial$winter.year[i] %in% nao$year) iceout_spatial$NAO[i] <- nao$index[nao$year == iceout_spatial$winter.year[i]]
  
  # June-July averages: winter.year=months 1-6 (cols 2-7), fall.year=months 7-12 (cols 8-13)
  if(iceout_spatial$winter.year[i] %in% enso$year & iceout_spatial$fall.year[i] %in% enso$year) iceout_spatial$ENSOjj[i] <- rowMeans(cbind(enso[enso$year == iceout_spatial$fall.year[i], 8:13], enso[enso$year == iceout_spatial$winter.year[i], 2:7]))
  if(iceout_spatial$winter.year[i] %in% qbo$year & iceout_spatial$fall.year[i] %in% qbo$year) iceout_spatial$QBOjj[i] <- rowMeans(cbind(qbo[qbo$year == iceout_spatial$fall.year[i], 8:13], qbo[qbo$year == iceout_spatial$winter.year[i], 2:7]))
  if(iceout_spatial$winter.year[i] %in% sun$year & iceout_spatial$fall.year[i] %in% sun$year) iceout_spatial$SUNjj[i] <- rowMeans(cbind(sun[sun$year == iceout_spatial$fall.year[i], 8:13], sun[sun$year == iceout_spatial$winter.year[i], 2:7]))
  if(iceout_spatial$winter.year[i] %in% pdo$year & iceout_spatial$fall.year[i] %in% pdo$year) iceout_spatial$PDOjj[i] <- rowMeans(cbind(pdo[pdo$year == iceout_spatial$fall.year[i], 8:13], pdo[pdo$year == iceout_spatial$winter.year[i], 2:7]))
}


### Detrend indices----
# "detrend" is in the package "pracma"

ice_spatial <- ice_spatial %>%
  mutate(ENSO_dt=detrend(ENSO, tt="linear"),
         NAO_dt=detrend(NAO, tt="linear"),
         QBO_dt=detrend(QBO, tt="linear"),
         SUN_dt=detrend(SUN, tt="linear"),
         PDO_dt=detrend(PDO, tt="linear"),
         ENSOjj_dt=detrend(ENSOjj, tt="linear"),
         QBOjj_dt=detrend(QBOjj, tt="linear"),
         SUNjj_dt=detrend(SUNjj, tt="linear"),
         PDOjj_dt=detrend(PDOjj, tt="linear"))

iceout_spatial <- iceout_spatial %>%
  mutate(ENSO_dt=detrend(ENSO, tt="linear"),
         NAO_dt=detrend(NAO, tt="linear"),
         QBO_dt=detrend(QBO, tt="linear"),
         SUN_dt=detrend(SUN, tt="linear"),
         PDO_dt=detrend(PDO, tt="linear"),
         ENSOjj_dt=detrend(ENSOjj, tt="linear"),
         QBOjj_dt=detrend(QBOjj, tt="linear"),
         SUNjj_dt=detrend(SUNjj, tt="linear"),
         PDOjj_dt=detrend(PDOjj, tt="linear"))

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

# Condense ice_all data
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
dim(iceout_w10y_condense) # 14291
range(iceout_w10y_condense$winter.year) # 1949-2022

length(unique(iceout_w30y_condense$ID))  # 163 lakes with 50 yrs ice off since 1949
dim(iceout_w30y_condense) # 7663 observations
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