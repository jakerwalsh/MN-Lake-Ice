# Variable phenology but consistent loss of ice cover of 1,213 Minnesota lakes
# Walsh et al. submitted to L&O Letters
# Model fitting

# Dependent on "1_Data set up.R"

# NOTE: These are big models that take a long time to fit. Models were fit using
# a machine with 32 GB of RAM. Even with more RAM, model crashing was not 
# uncommon without managing the RAM being used. 

# Model fitting ----
# Packages ----
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


# Fitting models with "method='ML'" for comparison----

# Gaussian Normal models ----

## Duration----
gam_all_duration_smoothci <- gam(min_duration ~ s(winter.year, k=5) + s(fwinter.year, bs='re') + s(lnArea_acres, k=5) + s(ID, bs='re') + s(x, y, k=5) + 
                                    s(ENSOjj, k=5) + s(NAO, k=5) + s(QBOjj, k=5) + s(SUNjj, k=5) + s(PDOjj, k=5),
                                  data = ice_all_condense,
                                  method='ML', select=TRUE)

## Ice on ----
gam_all_iceon_smoothci <- gam(max_ice_on_julian2 ~ s(winter.year, k=5) + s(fwinter.year, bs='re') + s(lnArea_acres, k=5) + s(ID, bs='re') + s(x, y, k=5) + 
                                 s(ENSOjj, k=5) + s(NAO, k=5) + s(QBOjj, k=5) + s(SUNjj, k=5) + s(PDOjj, k=5),
                               data = ice_all_condense,
                               method='ML', select=TRUE)

## Ice out ----
gam_all_iceout_smoothci <- gam(min_ice_off_julian ~ s(winter.year, k=5) + s(fwinter.year, bs='re') + s(lnArea_acres, k=5) + s(ID, bs='re') + s(x, y, k=5) + 
                                  s(ENSOjj, k=5) + s(NAO, k=5) + s(QBOjj, k=5) + s(SUNjj, k=5) + s(PDOjj, k=5),
                                data = ice_all_condense,
                                method='ML', select=TRUE)

## Ice out long ----
gam_all_iceoutlong_smoothci <- gam(min_ice_off_julian ~ s(winter.year, k=5) + s(fwinter.year, bs='re') + s(lnArea_acres, k=5) + s(ID, bs='re') + s(x, y, k=5) + 
                                      s(ENSOjj, k=5) + s(NAO, k=5) + s(QBOjj, k=5) + s(SUNjj, k=5) + s(PDOjj, k=5),
                                    data = iceout_all_condense,
                                    method='ML', select=TRUE)


# Heteroskedastic models ----

## Duration----
gamhs_all_duration_smoothci <- gam(list(min_duration ~ s(winter.year, k=5) + s(fwinter.year, bs='re') + s(lnArea_acres, k=5) + s(ID, bs='re') + s(x, y, k=5) + 
                                           s(ENSOjj, k=5) + s(NAO, k=5) + s(QBOjj, k=5) + s(SUNjj, k=5) + s(PDOjj, k=5),
                                         ~s(winter.year)), 
                                    data = ice_all_condense,
                                    family=gaulss(), 
                                    method='ML', select=TRUE)

## Ice on ----
gamhs_all_iceon_smoothci <- gam(list(max_ice_on_julian2 ~ s(winter.year, k=5) + s(fwinter.year, bs='re') + s(lnArea_acres, k=5) + s(ID, bs='re') + s(x, y, k=5) + 
                                        s(ENSOjj, k=5) + s(NAO, k=5) + s(QBOjj, k=5) + s(SUNjj, k=5) + s(PDOjj, k=5),
                                      ~s(winter.year)), 
                                 data = ice_all_condense,
                                 family=gaulss(), 
                                 method='ML', select=TRUE)

## Ice out ----
gamhs_all_iceout_smoothci <- gam(list(min_ice_off_julian ~ s(winter.year, k=5) + s(fwinter.year, bs='re') + s(lnArea_acres, k=5) + s(ID, bs='re') + s(x, y, k=5) + 
                                         s(ENSOjj, k=5) + s(NAO, k=5) + s(QBOjj, k=5) + s(SUNjj, k=5) + s(PDOjj, k=5),
                                       ~s(winter.year)), 
                                  data = ice_all_condense,
                                  family=gaulss(), 
                                  method='ML', select=TRUE)

## Ice out long ----
gamhs_all_iceoutlong_smoothci <- gam(list(min_ice_off_julian ~ s(winter.year, k=5) + s(fwinter.year, bs='re') + s(lnArea_acres, k=6) + s(ID, bs='re') + s(x, y, k=6) + 
                                             s(ENSOjj, k=5) + s(NAO, k=5) + s(QBOjj, k=5) + s(SUNjj, k=5) + s(PDOjj, k=5),
                                           ~s(winter.year)), 
                                      data = iceout_all_condense,
                                      family=gaulss(), 
                                      method='ML', select=TRUE)

# Creating a table to reference ML models based on components ----
ls()[grepl("^gam_", ls())|grepl("^gamhs_", ls())]
length(ls()[grepl("^gam_", ls())|grepl("^gamhs_", ls())])

ml_gam_models_meta <- data.frame(models=ls()[grepl("^gam_", ls())|grepl("^gamhs_", ls())]) %>%
  separate(col=models, into=c("approach", "data", "metric", "climate_index"), sep="_", remove=FALSE) %>%
  mutate(ci_smooth=ifelse(grepl("lin", climate_index), "Linear", 
                          ifelse(grepl("smooth", climate_index), "Smooth",
                                 ifelse(grepl("no", climate_index), "Not included", NA))),
         ci_window=ifelse(grepl("jd", climate_index), "January - December", "July - June"),
         metric2=ifelse(grepl("iceout_", models), "Day of ice breakup", 
                              ifelse(grepl("iceoutlong", metric), "Day of ice breakup (all data)",
                                     ifelse(grepl("iceon", metric), "Day of ice formation",
                                            ifelse(grepl("duration", metric), "Ice cover duration", NA)))),
         approach2=ifelse(grepl("^gam_", models), "Gaussian Normal", "Heteroskedastic"))

## Pulling model AICs ----
ml_gam_models_meta$aic <- NA

for(i in ml_gam_models_meta$models){
  
  ml_gam_models_meta$aic[ml_gam_models_meta$models==i] <- get(i)$aic
  
}

### TABLE S1: Maximum likelihood AICs ----
ml_gam_output <- ml_gam_models_meta %>%
  group_by(metric2, data) %>%
  arrange(aic, .by_group=TRUE) %>%
  dplyr::select(models, metric2, data, approach2, ci_smooth, ci_window, aic) %>%
  mutate(delta_aic=aic-min(aic)) %>%
  ungroup()

write.csv(x=ml_gam_output, file="tables/ml_gam_output.csv")



# Saving and unloading ML GAMs ----
save(list=ls()[grepl("^gam_", ls())|grepl("^gamhs_", ls())],
     file="outputdata/ml_gams_alldata.Rdata")

rm(list=ls()[grepl("^gam_", ls())|grepl("^gamhs_", ls())])

# Final models fit with REML ----
## Duration ----
gamreml_all_duration <- gam(list(min_duration ~ s(winter.year, k=5) + s(fwinter.year, bs='re') + s(lnArea_acres, k=5) + s(ID, bs='re') + s(x, y, k=5),
                                 ~s(winter.year)), 
                            data = ice_all_condense,
                            family=gaulss(),
                            method='REML')
gamreml_w10y_duration <- gam(list(min_duration ~ s(winter.year, k=5) + s(fwinter.year, bs='re') + s(lnArea_acres, k=5) + s(ID, bs='re') + s(x, y, k=5),
                                  ~s(winter.year)),
                             data = ice_w10y_condense,
                             family=gaulss(),
                             method='REML')
gamreml_w30y_duration <- gam(list(min_duration ~ s(winter.year, k=5) + s(fwinter.year, bs='re') + s(lnArea_acres, k=5) + s(ID, bs='re') + s(x, y, k=5),
                                  ~s(winter.year)), 
                             data = ice_w30y_condense,
                             family=gaulss(),
                             method='REML')

## Ice on ----
gamreml_all_iceon <- gam(list(max_ice_on_julian2 ~ s(winter.year, k=5) + s(fwinter.year, bs='re') + s(lnArea_acres, k=5) + s(ID, bs='re') + s(x, y, k=5),
                              ~s(winter.year)), 
                         data = ice_all_condense,
                         family=gaulss(),
                         method='REML')
gamreml_w10y_iceon <- gam(list(max_ice_on_julian2 ~ s(winter.year, k=5) + s(fwinter.year, bs='re') + s(lnArea_acres, k=5) + s(ID, bs='re') + s(x, y, k=5),
                               ~s(winter.year)), 
                          data = ice_w10y_condense,
                          family=gaulss(),
                          method='REML')
gamreml_w30y_iceon <- gam(list(max_ice_on_julian2 ~ s(winter.year, k=5) + s(fwinter.year, bs='re') + s(lnArea_acres, k=5) + s(ID, bs='re') + s(x, y, k=5),
                               ~s(winter.year)), 
                          data = ice_w30y_condense,
                          family=gaulss(),
                          method='REML')

## Ice out ----
gamreml_all_iceout <- gam(list(min_ice_off_julian ~ s(winter.year, k=5) + s(fwinter.year, bs='re') + s(lnArea_acres, k=5) + s(ID, bs='re') + s(x, y, k=5),
                               ~s(winter.year)), 
                          data = ice_all_condense,
                          family=gaulss(),
                          method='REML')
gamreml_w10y_iceout <- gam(list(min_ice_off_julian ~ s(winter.year, k=5) + s(fwinter.year, bs='re') + s(lnArea_acres, k=5) + s(ID, bs='re') + s(x, y, k=5),
                                ~s(winter.year)), 
                           data = ice_w10y_condense,
                           family=gaulss(),
                           method='REML')
gamreml_w30y_iceout <- gam(list(min_ice_off_julian ~ s(winter.year, k=5) + s(fwinter.year, bs='re') + s(lnArea_acres, k=5) + s(ID, bs='re') + s(x, y, k=5),
                                ~s(winter.year)), 
                           data = ice_w30y_condense,
                           family=gaulss(),
                           method='REML')

## Ice out long ----
gamreml_all_iceoutlong <- gam(list(min_ice_off_julian ~ s(winter.year, k=5) + s(fwinter.year, bs='re') + s(lnArea_acres, k=5) + s(ID, bs='re') + s(x, y, k=5),
                                   ~s(winter.year)), 
                              data = iceout_all_condense,
                              family=gaulss(),
                              method='REML')
gamreml_w10y_iceoutlong <- gam(list(min_ice_off_julian ~ s(winter.year, k=5) + s(fwinter.year, bs='re') + s(lnArea_acres, k=5) + s(ID, bs='re') + s(x, y, k=5),
                                    ~s(winter.year)), 
                               data = iceout_w10y_condense,
                               family=gaulss(),
                               method='REML')
gamreml_w30y_iceoutlong <- gam(list(min_ice_off_julian ~ s(winter.year, k=5) + s(fwinter.year, bs='re') + s(lnArea_acres, k=5) + s(ID, bs='re') + s(x, y, k=5),
                                    ~s(winter.year)), 
                               data = iceout_w30y_condense,
                               family=gaulss(),
                               method='REML')


# Lake-specific models----

gamreml_lakespecific_duration <- gam(list(min_duration ~ s(winter.year, k=5) + s(fwinter.year, bs='re') + 
                                            s(winter.year, by=ID, k=5) + 
                                            s(lnArea_acres, k=5) + s(ID, bs='re') + s(x, y, k=5),
                                          ~s(winter.year)), 
                                     data = ice_w10y1970_condense,
                                     family=gaulss(),
                                     method='REML', select=TRUE)
gamreml_lakespecific_iceon <- gam(list(max_ice_on_julian2 ~ s(winter.year, k=6) + s(fwinter.year, bs='re') + 
                                         s(winter.year, by=ID, k=6) + 
                                         s(lnArea_acres, k=6) + s(ID, bs='re') + s(x, y, k=6),
                                       ~s(winter.year)), 
                                  data = ice_w10y1970_condense,
                                  family=gaulss(),
                                  method='REML', select=TRUE)
gamreml_lakespecific_iceout <- gam(list(min_ice_off_julian ~ s(winter.year, k=5) + s(fwinter.year, bs='re') + 
                                          s(winter.year, by=ID, k=5) + 
                                          s(lnArea_acres, k=5) + s(ID, bs='re') + s(x, y, k=5),
                                        ~s(winter.year)), 
                                   data = ice_w10y1970_condense,
                                   family=gaulss(),
                                   method='REML', select=TRUE)
gamreml_lakespecific_iceoutlong <- gam(list(min_ice_off_julian ~ s(winter.year, k=5) + s(fwinter.year, bs='re') + 
                                              s(winter.year, by=ID, k=5) + 
                                              s(lnArea_acres, k=5) + s(ID, bs='re') + s(x, y, k=5),
                                            ~s(winter.year)), 
                                       data = ice_w10y1970_condense,
                                       family=gaulss(),
                                       method='REML', select=TRUE)


# Creating a table to reference final REML models based on components ----
ls()[grepl("gamreml_", ls())]
length(ls()[grepl("gamreml_", ls())])

reml_gam_final_models_meta <- data.frame(models=ls()[grepl("gamreml_", ls())]) %>%
  separate(col=models, into=c("approach", "data", "metric"), sep="_", remove=FALSE) %>%
  filter(data!="lakespecific")

reml_gam_final_lakespecific_models_meta <- data.frame(models=ls()[grepl("gamreml_lakespecific", ls())]) %>%
  separate(col=models, into=c("approach", "data", "metric"), sep="_", remove=FALSE)

# Saving final REML models ----
# These are needed for 3_Data visualization.R and 4_Supplemental Information.R,
# so aren't unloaded
save(list=ls()[grepl("gamreml", ls()) & !grepl("lakespecific", ls())], 
     file="outputdata/reml_gams.Rdata")

save(list=ls()[grepl("gamreml", ls()) & grepl("lakespecific", ls())], 
     file="outputdata/reml_gams_lakespecific.Rdata")