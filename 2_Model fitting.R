# Long-term trends in lake ice
# Hansen Lab Project (Revisit 2023)

# Dependent on "1_Data set up.R"

# Model fitting ----

# Model naming convention:
#     gam_LAKESUBSET_METRIC_INDEX or gam_LAKESUBSET_METRIC_INDEX_#
#     where "#" is explained in the code for each specific model if more specificity is needed

# Packages ----
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

# Fitting models with "method='ML'" for comparison----

# Gaussian Normal models ----

## Duration----
gam_all_duration_smoothci <- gam(min_duration ~ s(winter.year, k=5) + s(fwinter.year, bs='re') + s(lnArea_acres, k=5) + s(ID, bs='re') + s(x, y, k=5) + 
                                    s(ENSOjj, k=5) + s(NAO, k=5) + s(QBOjj, k=5) + s(SUNjj, k=5) + s(PDOjj, k=5),
                                  data = ice_all_condense,
                                  method='ML', select=TRUE)
save.image()
gc()

## Ice on ----
gam_all_iceon_smoothci <- gam(max_ice_on_julian2 ~ s(winter.year, k=5) + s(fwinter.year, bs='re') + s(lnArea_acres, k=5) + s(ID, bs='re') + s(x, y, k=5) + 
                                 s(ENSOjj, k=5) + s(NAO, k=5) + s(QBOjj, k=5) + s(SUNjj, k=5) + s(PDOjj, k=5),
                               data = ice_all_condense,
                               method='ML', select=TRUE)
save.image()
gc()

## Ice out ----
gam_all_iceout_smoothci <- gam(min_ice_off_julian ~ s(winter.year, k=5) + s(fwinter.year, bs='re') + s(lnArea_acres, k=5) + s(ID, bs='re') + s(x, y, k=5) + 
                                  s(ENSOjj, k=5) + s(NAO, k=5) + s(QBOjj, k=5) + s(SUNjj, k=5) + s(PDOjj, k=5),
                                data = ice_all_condense,
                                method='ML', select=TRUE)
save.image()
gc()

## Ice out long ----
gam_all_iceoutlong_smoothci <- gam(min_ice_off_julian ~ s(winter.year, k=5) + s(fwinter.year, bs='re') + s(lnArea_acres, k=5) + s(ID, bs='re') + s(x, y, k=5) + 
                                      s(ENSOjj, k=5) + s(NAO, k=5) + s(QBOjj, k=5) + s(SUNjj, k=5) + s(PDOjj, k=5),
                                    data = iceout_all_condense,
                                    method='ML', select=TRUE)
save.image()
gc()

# Heteroskedastic models: Over time ----

## Duration----
gamhs_all_duration_smoothci <- gam(list(min_duration ~ s(winter.year, k=5) + s(fwinter.year, bs='re') + s(lnArea_acres, k=5) + s(ID, bs='re') + s(x, y, k=5) + 
                                           s(ENSOjj, k=5) + s(NAO, k=5) + s(QBOjj, k=5) + s(SUNjj, k=5) + s(PDOjj, k=5),
                                         ~s(winter.year, k=10)), 
                                    data = ice_all_condense,
                                    family=gaulss(), 
                                    method='ML', select=TRUE)
save.image()
gc()

## Ice on ----
gamhs_all_iceon_smoothci <- gam(list(max_ice_on_julian2 ~ s(winter.year, k=5) + s(fwinter.year, bs='re') + s(lnArea_acres, k=5) + s(ID, bs='re') + s(x, y, k=5) + 
                                        s(ENSOjj, k=5) + s(NAO, k=5) + s(QBOjj, k=5) + s(SUNjj, k=5) + s(PDOjj, k=5),
                                      ~s(winter.year, k=10)), 
                                 data = ice_all_condense,
                                 family=gaulss(), 
                                 method='ML', select=TRUE)
save.image()
gc()

## Ice out ----
gamhs_all_iceout_smoothci <- gam(list(min_ice_off_julian ~ s(winter.year, k=5) + s(fwinter.year, bs='re') + s(lnArea_acres, k=5) + s(ID, bs='re') + s(x, y, k=5) + 
                                         s(ENSOjj, k=5) + s(NAO, k=5) + s(QBOjj, k=5) + s(SUNjj, k=5) + s(PDOjj, k=5),
                                       ~s(winter.year, k=10)), 
                                  data = ice_all_condense,
                                  family=gaulss(), 
                                  method='ML', select=TRUE)
save.image()
gc()

## Ice out long ----
gamhs_all_iceoutlong_smoothci <- gam(list(min_ice_off_julian ~ s(winter.year, k=5) + s(fwinter.year, bs='re') + s(lnArea_acres, k=6) + s(ID, bs='re') + s(x, y, k=6) + 
                                             s(ENSOjj, k=5) + s(NAO, k=5) + s(QBOjj, k=5) + s(SUNjj, k=5) + s(PDOjj, k=5),
                                           ~s(winter.year, k=10)), 
                                      data = iceout_all_condense,
                                      family=gaulss(), 
                                      method='ML', select=TRUE)
save.image()
gc()

# Heteroskedastic models: Over time and space ----

## Duration----
gamhss_all_duration_smoothci <- gam(list(min_duration ~ s(winter.year, k=5) + s(fwinter.year, bs='re') + s(lnArea_acres, k=5) + s(ID, bs='re') + s(x, y, k=5) + 
                                          s(ENSOjj, k=5) + s(NAO, k=5) + s(QBOjj, k=5) + s(SUNjj, k=5) + s(PDOjj, k=5),
                                        ~s(winter.year, k=10) + s(x, y, k=5)), 
                                   data = ice_all_condense,
                                   family=gaulss(), 
                                   method='ML', select=TRUE)
save.image()
gc()

### What about lake area for the SD term? ----
# Doesn't seem important enought to include, really just confirming the 
# spatial pattern is due to space rather than covariance between
# space and lake area.
## Duration----
gamhssa_all_duration_smoothci <- gam(list(min_duration ~ s(winter.year, k=5) + s(fwinter.year, bs='re') + s(lnArea_acres, k=5) + s(ID, bs='re') + s(x, y, k=5) + 
                                           s(ENSOjj, k=5) + s(NAO, k=5) + s(QBOjj, k=5) + s(SUNjj, k=5) + s(PDOjj, k=5),
                                         ~s(winter.year, k=10) + s(x, y, k=5) + s(lnArea_acres, k=5)), 
                                    data = ice_all_condense,
                                    family=gaulss(), 
                                    method='ML', select=TRUE)
save.image()
gc()

## Ice on ----
gamhss_all_iceon_smoothci <- gam(list(max_ice_on_julian2 ~ s(winter.year, k=5) + s(fwinter.year, bs='re') + s(lnArea_acres, k=5) + s(ID, bs='re') + s(x, y, k=5) + 
                                       s(ENSOjj, k=5) + s(NAO, k=5) + s(QBOjj, k=5) + s(SUNjj, k=5) + s(PDOjj, k=5),
                                     ~s(winter.year, k=10) + s(x, y, k=5)), 
                                data = ice_all_condense,
                                family=gaulss(), 
                                method='ML', select=TRUE)
save.image()
gc()

## Ice out ----
gamhss_all_iceout_smoothci <- gam(list(min_ice_off_julian ~ s(winter.year, k=5) + s(fwinter.year, bs='re') + s(lnArea_acres, k=5) + s(ID, bs='re') + s(x, y, k=5) + 
                                        s(ENSOjj, k=5) + s(NAO, k=5) + s(QBOjj, k=5) + s(SUNjj, k=5) + s(PDOjj, k=5),
                                      ~s(winter.year, k=10) + s(x, y, k=5)), 
                                 data = ice_all_condense,
                                 family=gaulss(), 
                                 method='ML', select=TRUE)
save.image()
gc()

## Ice out long ----
gamhss_all_iceoutlong_smoothci <- gam(list(min_ice_off_julian ~ s(winter.year, k=5) + s(fwinter.year, bs='re') + s(lnArea_acres, k=6) + s(ID, bs='re') + s(x, y, k=6) + 
                                            s(ENSOjj, k=5) + s(NAO, k=5) + s(QBOjj, k=5) + s(SUNjj, k=5) + s(PDOjj, k=5),
                                          ~s(winter.year, k=10) + s(x, y, k=5)), 
                                     data = iceout_all_condense,
                                     family=gaulss(), 
                                     method='ML', select=TRUE)
save.image()
gc()


# Creating a table to reference ML models based on components ----
ls()[grepl("^gam_", ls())|grepl("^gamhs_", ls())|grepl("^gamhss_", ls())|grepl("^gamhssa_", ls())]
length(ls()[grepl("^gam_", ls())|grepl("^gamhs_", ls())|grepl("^gamhss_", ls())|grepl("^gamhssa_", ls())])

ml_gam_models_meta <- data.frame(models=ls()[grepl("^gam_", ls())|grepl("^gamhs_", ls())|grepl("^gamhss_", ls())|grepl("^gamhssa_", ls())]) %>%
  separate(col=models, into=c("approach", "data", "metric", "climate_index"), sep="_", remove=FALSE) %>%
  mutate(ci_smooth=ifelse(grepl("lin", climate_index), "Linear", 
                          ifelse(grepl("smooth", climate_index), "Smooth",
                                 ifelse(grepl("no", climate_index), "Not included", NA))),
         ci_window=ifelse(grepl("jd", climate_index), "January - December", "July - June"),
         metric2=ifelse(grepl("iceout_", models), "Day of ice breakup", 
                              ifelse(grepl("iceoutlong", metric), "Day of ice breakup (all data)",
                                     ifelse(grepl("iceon", metric), "Day of ice formation",
                                            ifelse(grepl("duration", metric), "Ice cover duration", NA)))),
         approach2=ifelse(grepl("^gam_", models), "Gaussian Normal", 
                          ifelse(grepl("^gamhs_", models), "Heteroskedastic: Time", 
                                 ifelse(grepl("gamhss_", models), "Heteroskedastic: Time & Space", "Heteroskedastic: Time, Space, Area"))))

## Pulling model AICs ----
ml_gam_models_meta$aic <- NA

for(i in ml_gam_models_meta$models){
  
  ml_gam_models_meta$aic[ml_gam_models_meta$models==i] <- get(i)$aic
  
}

## Table: Maximum likelihood AICs ----
ml_gam_output <- ml_gam_models_meta %>%
  group_by(metric2, data) %>%
  arrange(aic, .by_group=TRUE) %>%
  dplyr::select(models, metric2, data, approach2, ci_smooth, ci_window, aic) %>%
  mutate(delta_aic=aic-min(aic)) %>%
  ungroup()

ml_gam_output %>%
  dplyr::select(models, aic, delta_aic) %>%
  as.data.frame()

write.csv(x=ml_gam_output, file="tables/ml_gam_output.csv")

# Saving and unloading ML GAMs ----

save(list=ls()[grepl("^gam_", ls())|grepl("^gamhs_", ls())|grepl("^gamhss_", ls())|grepl("^gamhssa_", ls())],
     file="outputdata/ml_gams_alldata.Rdata")

rm(list=ls()[grepl("^gam_", ls())|grepl("^gamhs_", ls())|grepl("^gamhss_", ls())|grepl("^gamhssa_", ls())])

gc()

# Final models fit with method = 'REML' ----

### Test k for spatial smooth ----
gamtest_all_duration_higherKforSpatial_k5 <- gam(list(min_duration ~ s(winter.year, k=5) + s(fwinter.year, bs='re') + s(lnArea_acres, k=5) + s(ID, bs='re') + s(x, y, k=5),
                                                      ~s(winter.year, k=10) + s(x, y, k=5)), 
                                                 data = ice_all_condense,
                                                 family=gaulss(),
                                                 method='REML')
save.image()
gam.check(gamtest_all_duration_higherKforSpatial_k5)
draw(gamtest_all_duration_higherKforSpatial_k5, select='s(x,y)')
draw(gamtest_all_duration_higherKforSpatial_k5, select='s.1(x,y)')

gamtest_all_duration_higherKforSpatial_k10 <- gam(list(min_duration ~ s(winter.year, k=5) + s(fwinter.year, bs='re') + s(lnArea_acres, k=5) + s(ID, bs='re') + s(x, y, k=10),
                                                       ~s(winter.year, k=10) + s(x, y, k=5)), 
                                                  data = ice_all_condense,
                                                  family=gaulss(),
                                                  method='REML')
save.image()
gam.check(gamtest_all_duration_higherKforSpatial_k10)
draw(gamtest_all_duration_higherKforSpatial_k10, select='s(x,y)')
draw(gamtest_all_duration_higherKforSpatial_k10, select='s.1(x,y)')

gamtest_all_duration_higherKforSpatial_k20 <- gam(list(min_duration ~ s(winter.year, k=5) + s(fwinter.year, bs='re') + s(lnArea_acres, k=5) + s(ID, bs='re') + s(x, y, k=20),
                                                       ~s(winter.year, k=10) + s(x, y, k=5)), 
                                                  data = ice_all_condense,
                                                  family=gaulss(),
                                                  method='REML')
save.image()
gam.check(gamtest_all_duration_higherKforSpatial_k20)
draw(gamtest_all_duration_higherKforSpatial_k20, select='s(x,y)')
draw(gamtest_all_duration_higherKforSpatial_k20, select='s.1(x,y)')

gamtest_all_duration_higherKforSpatial_k30 <- gam(list(min_duration ~ s(winter.year, k=5) + s(fwinter.year, bs='re') + s(lnArea_acres, k=5) + s(ID, bs='re') + s(x, y, k=30),
                                                       ~s(winter.year, k=10) + s(x, y, k=5)), 
                                                  data = ice_all_condense,
                                                  family=gaulss(),
                                                  method='REML')
save.image()
gam.check(gamtest_all_duration_higherKforSpatial_k30)
draw(gamtest_all_duration_higherKforSpatial_k30, select='s(x,y)')
draw(gamtest_all_duration_higherKforSpatial_k30, select='s.1(x,y)')

# comparison plot

plot_grid(draw(gamtest_all_duration_higherKforSpatial_k5, select='s(x,y)'), 
          draw(gamtest_all_duration_higherKforSpatial_k10, select='s(x,y)'),
          draw(gamtest_all_duration_higherKforSpatial_k20, select='s(x,y)'),
          draw(gamtest_all_duration_higherKforSpatial_k30, select='s(x,y)'),
          labels=c("k=5", "k=10", "k=20", "k=30"))


### Duration ----
gamreml_w10y_duration <- gam(list(min_duration ~ s(winter.year, k=5) + s(fwinter.year, bs='re') + s(lnArea_acres, k=5) + s(ID, bs='re') + s(x, y, k=20),
                                  ~s(winter.year, k=10) + s(x, y, k=5)),
                             data = ice_w10y_condense,
                             family=gaulss(),
                             method='REML')
gamreml_w30y_duration <- gam(list(min_duration ~ s(winter.year, k=5) + s(fwinter.year, bs='re') + s(lnArea_acres, k=5) + s(ID, bs='re') + s(x, y, k=20),
                                  ~s(winter.year, k=10) + s(x, y, k=5)), 
                             data = ice_w30y_condense,
                             family=gaulss(),
                             method='REML')

### Ice on ----
gamreml_w10y_iceon <- gam(list(max_ice_on_julian2 ~ s(winter.year, k=5) + s(fwinter.year, bs='re') + s(lnArea_acres, k=5) + s(ID, bs='re') + s(x, y, k=20),
                               ~s(winter.year, k=10) + s(x, y, k=5)), 
                          data = ice_w10y_condense,
                          family=gaulss(),
                          method='REML')
gamreml_w30y_iceon <- gam(list(max_ice_on_julian2 ~ s(winter.year, k=5) + s(fwinter.year, bs='re') + s(lnArea_acres, k=5) + s(ID, bs='re') + s(x, y, k=20),
                               ~s(winter.year, k=10) + s(x, y, k=5)), 
                          data = ice_w30y_condense,
                          family=gaulss(),
                          method='REML')

### Ice out ----
gamreml_w10y_iceout <- gam(list(min_ice_off_julian ~ s(winter.year, k=5) + s(fwinter.year, bs='re') + s(lnArea_acres, k=5) + s(ID, bs='re') + s(x, y, k=20),
                                ~s(winter.year, k=10) + s(x, y, k=5)), 
                           data = ice_w10y_condense,
                           family=gaulss(),
                           method='REML')
gamreml_w30y_iceout <- gam(list(min_ice_off_julian ~ s(winter.year, k=5) + s(fwinter.year, bs='re') + s(lnArea_acres, k=5) + s(ID, bs='re') + s(x, y, k=20),
                                ~s(winter.year, k=10) + s(x, y, k=5)), 
                           data = ice_w30y_condense,
                           family=gaulss(),
                           method='REML')

### Ice out long ----
gamreml_w10y_iceoutlong <- gam(list(min_ice_off_julian ~ s(winter.year, k=5) + s(fwinter.year, bs='re') + s(lnArea_acres, k=5) + s(ID, bs='re') + s(x, y, k=20),
                                    ~s(winter.year, k=10) + s(x, y, k=5)), 
                               data = iceout_w10y_condense,
                               family=gaulss(),
                               method='REML')
gamreml_w30y_iceoutlong <- gam(list(min_ice_off_julian ~ s(winter.year, k=5) + s(fwinter.year, bs='re') + s(lnArea_acres, k=5) + s(ID, bs='re') + s(x, y, k=20),
                                    ~s(winter.year, k=10) + s(x, y, k=5)), 
                               data = iceout_w30y_condense,
                               family=gaulss(),
                               method='REML')
save.image()

## All lakes--SLOW-FITTING, HUGE GAMs ----

# duration
gamreml_all_duration <- gam(list(min_duration ~ s(winter.year, k=5) + s(fwinter.year, bs='re') + s(lnArea_acres, k=5) + s(ID, bs='re') + s(x, y, k=20),
                                 ~s(winter.year, k=10) + s(x, y, k=5)), 
                            data = ice_all_condense,
                            family=gaulss(),
                            method='REML')
save.image()

# iceon
gamreml_all_iceon <- gam(list(max_ice_on_julian2 ~ s(winter.year, k=5) + s(fwinter.year, bs='re') + s(lnArea_acres, k=5) + s(ID, bs='re') + s(x, y, k=20),
                              ~s(winter.year, k=10) + s(x, y, k=5)), 
                         data = ice_all_condense,
                         family=gaulss(),
                         method='REML')
save.image()

# iceout
gamreml_all_iceout <- gam(list(min_ice_off_julian ~ s(winter.year, k=5) + s(fwinter.year, bs='re') + s(lnArea_acres, k=5) + s(ID, bs='re') + s(x, y, k=20),
                               ~s(winter.year, k=10) + s(x, y, k=5)), 
                          data = ice_all_condense,
                          family=gaulss(),
                          method='REML')
save.image()

# iceoutlong
gamreml_all_iceoutlong <- gam(list(min_ice_off_julian ~ s(winter.year, k=5) + s(fwinter.year, bs='re') + s(lnArea_acres, k=5) + s(ID, bs='re') + s(x, y, k=20),
                                   ~s(winter.year, k=10) + s(x, y, k=5)), 
                              data = iceout_all_condense,
                              family=gaulss(),
                              method='REML')
save.image()


# Lake-specific models----

gamreml_lakespecific_duration <- gam(list(min_duration ~ s(winter.year, k=5) + s(fwinter.year, bs='re') + 
                                            s(winter.year, by=ID, k=5) + 
                                            s(lnArea_acres, k=5) + s(ID, bs='re') + s(x, y, k=20),
                                          ~s(winter.year, k=10) + s(x, y, k=5)), 
                                     data = ice_w10y1970_condense,
                                     family=gaulss(),
                                     method='REML', select=TRUE)
save.image()
gc()

gamreml_lakespecific_iceon <- gam(list(max_ice_on_julian2 ~ s(winter.year, k=6) + s(fwinter.year, bs='re') +
                                         s(winter.year, by=ID, k=6) +
                                         s(lnArea_acres, k=6) + s(ID, bs='re') + s(x, y, k=20),
                                       ~s(winter.year, k=10) + s(x, y, k=6)),
                                  data = ice_w10y1970_condense,
                                  family=gaulss(),
                                  method='REML', select=TRUE)
save.image()
gc()

gamreml_lakespecific_iceout <- gam(list(min_ice_off_julian ~ s(winter.year, k=5) + s(fwinter.year, bs='re') + 
                                          s(winter.year, by=ID, k=5) + 
                                          s(lnArea_acres, k=5) + s(ID, bs='re') + s(x, y, k=20),
                                        ~s(winter.year, k=10) + s(x, y, k=5)), 
                                   data = ice_w10y1970_condense,
                                   family=gaulss(),
                                   method='REML', select=TRUE)
save.image()
gc()

gamreml_lakespecific_iceoutlong <- gam(list(min_ice_off_julian ~ s(winter.year, k=5) + s(fwinter.year, bs='re') + 
                                              s(winter.year, by=ID, k=5) + 
                                              s(lnArea_acres, k=5) + s(ID, bs='re') + s(x, y, k=20),
                                            ~s(winter.year, k=10) + s(x, y, k=5)), 
                                       data = ice_w10y1970_condense,
                                       family=gaulss(),
                                       method='REML', select=TRUE)
save.image()
gc()

# Creating a table to reference final REML models based on components ----
ls()[grepl("gamreml_", ls())]
length(ls()[grepl("gamreml_", ls())])

reml_gam_final_models_meta <- data.frame(models=ls()[grepl("gamreml_", ls())]) %>%
  separate(col=models, into=c("approach", "data", "metric"), sep="_", remove=FALSE) %>%
  filter(data!="lakespecific")

reml_gam_final_lakespecific_models_meta <- data.frame(models=ls()[grepl("gamreml_lakespecific", ls())]) %>%
  separate(col=models, into=c("approach", "data", "metric"), sep="_", remove=FALSE)

# Saving final REML models ----
save(list=ls()[grepl("gamreml", ls()) & !grepl("lakespecific", ls())], 
     file="outputdata/reml_gams.Rdata")

save(list=ls()[grepl("gamreml", ls()) & grepl("lakespecific", ls())], 
     file="outputdata/reml_gams_lakespecific.Rdata")