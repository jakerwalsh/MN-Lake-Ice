# Variable phenology but consistent loss of ice cover of 1,213 Minnesota lakes
# Walsh et al. submitted to L&O Letters
# Data visualization

# Dependent on "1_Data set up.R", "2_Model fitting.R"

# Data visualization ----

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

# Model predictors ----
reml_gam_final_models_meta

fun_confint <- function(object, parm){
  
  gam_i <- get(object)
  
  confint_i <- confint(object=gam_i, parm=parm) %>%
    mutate(models=object)
  
  colnames(confint_i)[4] <- "parm"
  
  return(confint_i)
  
}

fun_confint_xy <- function(object, parm){
  
  gam_i <- get(object)
  
  confint_i <- confint(object=gam_i, parm=parm) %>%
    mutate(models=object)
  
  return(confint_i)
  
}

fun_confint_SD <- function(object, parm){
  
  gam_i <- get(object)
  
  intercept_SD_i <- coef(gam_i)["(Intercept).1"]
  
  confint_i <- confint(object=gam_i, parm=parm) %>%
    mutate(models=object,
           intercept_SD=exp(intercept_SD_i)+0.01,
           est_shift=exp(est+intercept_SD_i)+0.01,
           lower_shift=exp(lower+intercept_SD_i)+0.01,
           upper_shift=exp(upper+intercept_SD_i)+0.01)
  
  colnames(confint_i)[4] <- "parm"
  
  return(confint_i)
  
}

## winter.year ----
predictor_winter.year <- bind_rows(lapply(X=reml_gam_final_models_meta$models, parm="s(winter.year)", FUN=fun_confint))
predictor_winter.year <- left_join(predictor_winter.year, reml_gam_final_models_meta, by=c("models"))

p_predictor_winter.year <- ggplot(predictor_winter.year %>% filter(data=='all' & metric != "iceoutlong"), 
       aes(x=parm, y=est)) + 
  facet_grid(metric~.) + 
  geom_hline(yintercept=0, col='grey') + 
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2) + 
  geom_line() + 
  ggtitle("a) s(winter.year)") + 
  scale_y_continuous(name="Days from average") +
  scale_x_continuous(name="winter.year (continuous)") + 
  theme_classic(8) + theme()

## fwinter.year ----
predictor_fwinter.year <- bind_rows(lapply(X=reml_gam_final_models_meta$models, parm="s(fwinter.year)", FUN=fun_confint))
predictor_fwinter.year <- left_join(predictor_fwinter.year, reml_gam_final_models_meta, by=c("models"))

p_predictor_fwinter.year <- ggplot(predictor_fwinter.year %>% filter(data=='all' & metric != "iceoutlong"), 
       aes(x=as.numeric(as.character(parm)), y=est)) + 
  facet_grid(metric~.) + 
  geom_hline(yintercept=0, col='grey') + 
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0, col='grey') + 
  geom_point() + 
  geom_smooth(method='gam', formula=y~s(x,k=20), col='dodgerblue') + 
  ggtitle("b) s(fwinter.year, bs='re')") + 
  scale_y_continuous(name="Days from average") +
  scale_x_continuous(name="fwinter.year (random intercept)") + 
  theme_classic(8) + theme()

## lnArea_acres ----
predictor_lnArea_acres <- bind_rows(lapply(X=reml_gam_final_models_meta$models, parm="s(lnArea_acres)", FUN=fun_confint))
predictor_lnArea_acres <- left_join(predictor_lnArea_acres, reml_gam_final_models_meta, by=c("models"))

p_predictor_lnArea_acres <- ggplot(predictor_lnArea_acres %>% filter(data=='all' & metric != "iceoutlong"), 
       aes(x=parm, y=est)) + 
  facet_grid(metric~.) + 
  geom_hline(yintercept=0, col='grey') + 
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2) + 
  geom_line() + 
  ggtitle("d) s(lnArea_acres)") + 
  scale_y_continuous(name="Days from average") +
  scale_x_continuous(name="lnArea_acres") + 
  theme_classic(8) + theme()

## ID ----
predictor_ID <- bind_rows(lapply(X=reml_gam_final_models_meta$models, parm="s(ID)", FUN=fun_confint))
predictor_ID <- left_join(predictor_ID, reml_gam_final_models_meta, by=c("models")) %>% 
  group_by(models) %>%
  arrange(desc(est)) %>%
  mutate(rank=dense_rank(desc(est))) %>%
  ungroup()

p_predictor_ID <- ggplot(predictor_ID %>% filter(data=='all' & metric != "iceoutlong"), 
       aes(x=rank, y=est)) + 
  facet_grid(metric~.) + 
  geom_hline(yintercept=0, col='grey') + 
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0, col='grey') + 
  geom_point() + 
  ggtitle("e) s(ID, bs='re')") + 
  scale_y_continuous(name="Days from average") +
  scale_x_continuous(name="Ranked random intercept") + 
  theme_classic(8) + theme()

## xy ----
 
predictor_xy <- bind_rows(lapply(X=reml_gam_final_models_meta$models, parm="s(x,y)", FUN=fun_confint_xy))
predictor_xy <- left_join(predictor_xy, reml_gam_final_models_meta, by=c("models"))

p_predictor_xy <- ggplot(predictor_xy %>% filter(data=='all' & metric != "iceoutlong"), 
                         aes(x=x/10000, y=y/100000, fill=est)) + 
  facet_grid(metric~.) + 
  geom_raster(interpolate=T) + 
  ggtitle("f) s(x,y)") + 
  scale_y_continuous(name="y*100,000") +
  scale_x_continuous(name="x*10,000") + 
  scale_fill_gradient2(name="Days\nfrom avg.") + 
  theme_classic(8) + theme() 

## winter.year SD ----
predictor_winter.year_SD <- bind_rows(lapply(X=reml_gam_final_models_meta$models, parm="s.1(winter.year)", FUN=fun_confint_SD))
predictor_winter.year_SD <- left_join(predictor_winter.year_SD, reml_gam_final_models_meta, by=c("models"))

p_predictor_winter.year_SD <- ggplot(predictor_winter.year_SD %>% filter(data=='all' & metric != "iceoutlong"), 
                                  aes(x=parm, y=est_shift)) + 
  facet_grid(metric~.) + 
  #geom_hline(yintercept=0, col='grey') + 
  geom_ribbon(aes(ymin=lower_shift, ymax=upper_shift), alpha=0.2) + 
  geom_line() + 
  ggtitle("c) S.D. s(winter.year)") + 
  scale_y_continuous(name="S.D.") +
  scale_x_continuous(name="winter.year (continuous)") + 
  theme_classic(8) + theme()

### FIGURE S1--Multi-panel plot of predictors ----

p_predictor_all <- plot_grid(p_predictor_winter.year, p_predictor_fwinter.year, p_predictor_winter.year_SD,
                             p_predictor_lnArea_acres, p_predictor_ID, p_predictor_xy,
          nrow=2, rel_widths = c(1, 1, 1, 1, 1, 1.25))
  
png("figures/FigureS1_ModelPredictors.png", width=7, height=5, units='in', res=1200)
p_predictor_all
dev.off()


# Increasing variance over time plot ----

variance_over_time_df <- expand.grid(winter.year=1953:2022,
                                     metric=c("duration", "iceon", "iceout"),
                                      sd_resid=NA,
                                      sd_fwinter.year=NA,
                                      sd_data=NA)

for(i in 1953:2022){
  #duration
  variance_over_time_df$sd_fwinter.year[variance_over_time_df$winter.year==i & variance_over_time_df$metric=="duration"] <- sd(predictor_fwinter.year$est[as.numeric(as.character(predictor_fwinter.year$parm)) <= i & 
                                                                                                                         as.numeric(as.character(predictor_fwinter.year$parm)) > i - 5 & 
                                                                                                                           predictor_fwinter.year$metric=="duration"])
  
  variance_over_time_df$sd_resid[variance_over_time_df$winter.year==i & variance_over_time_df$metric=="duration"] <- sd(resid(gamreml_all_duration)[gamreml_all_duration$model$winter.year <= i &
                                                                                                                              gamreml_all_duration$model$winter.year > i - 5])
  
  variance_over_time_df$sd_data[variance_over_time_df$winter.year==i & variance_over_time_df$metric=="duration"] <- sd(ice_all_condense$min_duration[ice_all_condense$winter.year<= i 
                                                                                                            & ice_all_condense$winter.year > i-5])
  
  #iceon
  variance_over_time_df$sd_fwinter.year[variance_over_time_df$winter.year==i & variance_over_time_df$metric=="iceon"] <- sd(predictor_fwinter.year$est[as.numeric(as.character(predictor_fwinter.year$parm)) <= i & 
                                                                                                                                                            as.numeric(as.character(predictor_fwinter.year$parm)) > i - 5 & 
                                                                                                                                                         predictor_fwinter.year$metric=="iceon"])
  
  variance_over_time_df$sd_resid[variance_over_time_df$winter.year==i & variance_over_time_df$metric=="iceon"] <- sd(resid(gamreml_all_iceon)[gamreml_all_iceon$model$winter.year <= i &
                                                                                                                                                       gamreml_all_iceon$model$winter.year > i - 5])
  
  variance_over_time_df$sd_data[variance_over_time_df$winter.year==i & variance_over_time_df$metric=="iceon"] <- sd(ice_all_condense$max_ice_on_julian2[ice_all_condense$winter.year<= i 
                                                                                                                                                      & ice_all_condense$winter.year > i-5])
  
  
  #iceout
  variance_over_time_df$sd_fwinter.year[variance_over_time_df$winter.year==i & variance_over_time_df$metric=="iceout"] <- sd(predictor_fwinter.year$est[as.numeric(as.character(predictor_fwinter.year$parm)) <= i & 
                                                                                                                                                         as.numeric(as.character(predictor_fwinter.year$parm)) > i - 5 & 
                                                                                                                                                         predictor_fwinter.year$metric=="iceout"])
  
  variance_over_time_df$sd_resid[variance_over_time_df$winter.year==i & variance_over_time_df$metric=="iceout"] <- sd(resid(gamreml_all_iceout)[gamreml_all_iceout$model$winter.year <= i &
                                                                                                                                                    gamreml_all_iceout$model$winter.year > i - 5])
  
  variance_over_time_df$sd_data[variance_over_time_df$winter.year==i & variance_over_time_df$metric=="iceout"] <- sd(ice_all_condense$min_ice_off_julian[ice_all_condense$winter.year<= i 
                                                                                                                                                         & ice_all_condense$winter.year > i-5])
  
  
  
  print(i)
  
}

variance_over_time_df <- variance_over_time_df %>%
  pivot_longer(sd_resid:sd_data, names_to="variable", values_to="SD")

variance_over_time_df$variable2 <- factor(ifelse(variance_over_time_df$variable=="sd_resid", "Residuals",
                                                  ifelse(variance_over_time_df$variable=="sd_fwinter.year", "Winter Year R.E.",
                                                         ifelse(variance_over_time_df$variable=="sd_data", "Raw Data", NA))),
                                           levels=c("Raw Data", "Winter Year R.E.", "Residuals"))

variance_over_time_df$metric2 <- factor(ifelse(variance_over_time_df$metric=="duration", "Duration",
                                               ifelse(variance_over_time_df$metric=="iceon", "Formation",
                                                      ifelse(variance_over_time_df$metric=="iceout", "Breakup", NA))),
                                        levels=c("Duration", "Formation", "Breakup"))

p_varianceincrease_rawdata <- ggplot(data=variance_over_time_df %>% filter(variable2=="Raw Data"), 
                                      aes(x=winter.year, y=SD, col=metric2, lty=metric2)) +
  facet_grid(variable2~., scales='free_y') + 
  geom_line() + 
  scale_x_continuous("Winter Year") + 
  scale_y_continuous(expression("S.D."["5-year"])) + 
  scale_linetype_manual("Metric", values=c(1, 3, 2)) + 
  scale_color_manual("Metric", values=c("black", "orange", "lightblue")) + 
  theme_classic(9) + theme(axis.title.x=element_blank(),
                           legend.position='top')

p_varianceincrease_re <- ggplot(data=variance_over_time_df %>% filter(variable2=="Winter Year R.E."), 
                                 aes(x=winter.year, y=SD, col=metric2, lty=metric2)) +
  facet_grid(variable2~., scales='free_y') + 
  geom_line() + 
  scale_x_continuous("Winter Year") + 
  scale_y_continuous(expression("S.D."["5-year"])) + 
  scale_linetype_manual("Metric", values=c(1, 3, 2)) + 
  scale_color_manual("Metric", values=c("black", "orange", "lightblue")) + 
  theme_classic(9) + theme(axis.title.x=element_blank(),
                           strip.background = element_rect(fill='lightgreen'),
                           legend.position='none')

p_varianceincrease_sdterm <- ggplot(data=predictor_winter.year_SD %>% 
                      filter(data=="all" & parm>=1953 & metric!="iceoutlong") %>%
                      mutate(variable2="Smooth Function",
                             metric2=factor(ifelse(metric=="duration", "Duration",
                                                   ifelse(metric=="iceon", "Formation",
                                                          ifelse(metric=="iceout", "Breakup", NA))),
                                            levels=c("Duration", "Formation", "Breakup"))), 
                    aes(x=parm, y=est_shift, lty=metric2, col=metric2, fill=metric2)) + 
  facet_grid(variable2~.) + 
  geom_ribbon(aes(ymin=lower_shift, ymax=upper_shift), alpha=0.2, col=NA) + 
  geom_line() + 
  scale_x_continuous("Winter Year") + 
  scale_y_continuous(expression(sigma[italic("i,t")]), 
                     breaks=seq(3, 9, 2), 
                     labels=c("   3", "   5", "   7", "   9")) +
  scale_linetype_manual("Metric", values=c(1, 3, 2)) + 
  scale_color_manual("Metric", values=c("black", "orange", "lightblue")) + 
  scale_fill_manual("Metric", values=c("black", "orange", "lightblue")) + 
  theme_classic(9) + theme(axis.title.x=element_blank(),
                           strip.background=element_rect(fill='lightgreen'),
                           legend.position='none')

p_varianceincrease_resid <- ggplot(data=variance_over_time_df %>% filter(variable2=="Residuals"), 
                                    aes(x=winter.year, y=SD, lty=metric2, col=metric2)) +
  facet_grid(variable2~., scales='free_y') + 
  geom_line() + 
  scale_x_continuous("Winter Year") + 
  scale_y_continuous(expression("S.D."["5-year"])) + 
  scale_linetype_manual("Metric", values=c(1, 3, 2)) + 
  scale_color_manual("Metric", values=c("black", "orange", "lightblue")) + 
  theme_classic(9) + theme(strip.background = element_rect(fill='lightyellow'),
                           legend.position='none')

## FIGURE 1--Increasing variance over time ----
png("figures/Figure1_IncreasingVarianceOverTime.png", width=3.5, height=6, units='in', res=1200)
plot_grid(p_varianceincrease_rawdata,
          p_varianceincrease_re,
          p_varianceincrease_sdterm,
          p_varianceincrease_resid,
          nrow=4, labels='auto', label_size=9)
dev.off()


# Trends ----

## Trend estimation ----
trends_df <- expand.grid(models=reml_gam_final_models_meta$models,
                        winter.year=1949:2022) %>%
  mutate(est=NA, se=NA)

for(i in reml_gam_final_models_meta$models){
  
  trends_df$est[trends_df$models==i] <- predict(get(i), 
                                              newdata=data.frame(winter.year=1949:2022), 
                                              exclude=c("s(fwinter.year)", "s(lnArea_acres)", "s(ID)", "s(x,y)", "s.1(winter.year)"),
                                              newdata.guaranteed=T)[,1]
  
  trends_df$se[trends_df$models==i] <- predict(get(i), 
                                              newdata=data.frame(winter.year=1949:2022), 
                                              exclude=c("s(fwinter.year)", "s(lnArea_acres)", "s(ID)", "s(x,y)", "s.1(winter.year)"),
                                             newdata.guaranteed=T,
                                             se.fit=T)$se.fit[,1]
  
}

trends_df <- trends_df %>%
  mutate(lower = est - 2*se, upper = est + 2*se)

trends_df <- left_join(trends_df, reml_gam_final_models_meta, by="models")

fun_trends <- function(object){
  
  gam_i <- get(object)
  
  trends_i <- confint(object=gam_i, parm="s(winter.year)", n=74, type='simultaneous', shift=TRUE) %>%
    mutate(models=object)
  
  return(trends_i)
  
}

trends_df <- bind_rows(lapply(X=reml_gam_final_models_meta$models, FUN=fun_trends)) %>%
  dplyr::select(models, winter.year, est, se, lower, upper) %>%
  group_by(models, winter.year) %>%
  summarize_all(mean)

trends_df <- left_join(trends_df, reml_gam_final_models_meta, by='models')


## Derivative estimation ----

fun_derivatives <- function(object){
  
  gam_i <- get(object)
  
  derivatives_i <- derivatives(object=gam_i, term="s(winter.year)", n=73, eps=1, type='backward', interval='simultaneous') %>%
    mutate(models=object)
  
  return(derivatives_i)
  
}

derivatives_df <- bind_rows(lapply(X=reml_gam_final_models_meta$models, FUN=fun_derivatives)) %>%
  mutate(winter.year=data) %>%
  dplyr::select(models, winter.year, derivative, se, lower, upper) %>%
  group_by(models, winter.year) %>%
  summarize_all(mean)

derivatives_df <- left_join(derivatives_df, reml_gam_final_models_meta, by='models')

# setting up the data frame to calculate relative change (% change)

pred_df <- expand.grid(models=unique(derivatives_df$models),
                       winter.year=1949:2022)

pred_df$est <- NA
pred_df$est_se <- NA

for(i in unique(derivatives_df$models)){
  
  pred_df$est[pred_df$models==i] <- predict(get(i), 
                                             newdata=data.frame(winter.year=pred_df$winter.year[pred_df$models==i]),
                                             exclude=c("s(fwinter.year)", "s(lnArea_acres)", "s(ID)", "s(x,y)", "s.1(winter.year)"),
                                             newdata.guaranteed=TRUE)[,1]
  
  pred_df$est_se[pred_df$models==i] <- predict(get(i), 
                                            newdata=data.frame(winter.year=pred_df$winter.year[pred_df$models==i]),
                                            exclude=c("s(fwinter.year)", "s(lnArea_acres)", "s(ID)", "s(x,y)", "s.1(winter.year)"),
                                            newdata.guaranteed=TRUE, se.fit=TRUE)$se.fit[,1]
  
}

derivatives_df <- full_join(derivatives_df, 
                            pred_df,
                            by=c("models", "winter.year"))

derivatives_df$perc_change <- NA
derivatives_df$perc_change_lower <- NA
derivatives_df$perc_change_upper <- NA

for(i in 1:dim(derivatives_df)[1]) {
  
  if(derivatives_df$winter.year[i] > 1949){
    
    derivatives_df$perc_change[i] <- (derivatives_df$derivative[i]/
      derivatives_df$est[derivatives_df$winter.year==derivatives_df$winter.year[i]-1 &
                           derivatives_df$models==derivatives_df$models[i]])*100
    
    derivatives_df$perc_change_lower[i] <- (derivatives_df$lower[i]/
                                        derivatives_df$est[derivatives_df$winter.year==derivatives_df$winter.year[i]-1 &
                                                             derivatives_df$models==derivatives_df$models[i]])*100
    
    derivatives_df$perc_change_upper[i] <- (derivatives_df$upper[i]/
                                        derivatives_df$est[derivatives_df$winter.year==derivatives_df$winter.year[i]-1 &
                                                             derivatives_df$models==derivatives_df$models[i]])*100
    
  }
  
}

## Trends plot for main text ----

### duration

p_trend_duration_trend <- ggplot(trends_df%>%filter(metric=="duration"& data=="all"), 
                                 aes(x=winter.year, y=est, lty=data)) + 
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2, col=NA) + 
    geom_line() + 
    ggtitle("Ice cover duration") + 
    scale_y_continuous("Days", breaks=c(135, 145, 155), labels=c(135, 145,155)) + 
    scale_linetype_discrete("Lake data") + 
    theme_classic(9) + theme(legend.position="none",
                             axis.title.x=element_blank(),
                             plot.margin=unit(c(1,0.25, 0.25, 0.25), units='lines'))
  
p_trend_duration_derivative <- ggplot(derivatives_df%>%filter(metric=="duration"& data=="all"), 
                                      aes(x=winter.year, y=derivative, lty=data)) + 
    geom_hline(yintercept=0) + 
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2, col=NA) + 
    geom_line() + 
    scale_y_continuous("Days/year", limits=c(-0.43, 0.43)) + 
    scale_x_continuous("Winter Year") + 
    scale_linetype_discrete("Lake data") + 
    theme_classic(9) + theme(legend.position='none', axis.title.x=element_blank(), plot.margin=unit(c(2,0.25, 0.25, 0.25), units="lines"))

### iceon

p_trend_iceon_trend <- ggplot(trends_df%>%filter(metric=="iceon"& data=="all"), 
                              aes(x=winter.year, y=est, lty=data)) + 
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2, col=NA) + 
    geom_line() + 
    ggtitle("Day of ice formation") + 
    scale_y_continuous("Day of year") + 
    scale_linetype_discrete("Lake data") + 
    theme_classic(9) + theme(legend.position='none', axis.title.x=element_blank(), plot.margin=unit(c(1,0.25, 0.25, 0.25), units="lines"))
  
p_trend_iceon_derivative <- ggplot(derivatives_df%>%filter(metric=="iceon"& data=="all"), 
                                   aes(x=winter.year, y=derivative, lty=data)) + 
    geom_hline(yintercept=0) + 
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2, col=NA) + 
    geom_line() + 
    scale_y_continuous("Days/year", limits=c(-0.43, 0.43)) + 
    scale_x_continuous("Winter Year") + 
    scale_linetype_discrete("Lake data") + 
    theme_classic(9) + theme(legend.position='none', axis.title.x=element_blank(), plot.margin=unit(c(2,0.25, 0.25, 0.25), units="lines"))

### iceout

p_trend_iceout_trend <- ggplot(trends_df%>%filter(metric=="iceout"& data=="all"), 
                               aes(x=winter.year, y=est, lty=data)) + 
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2, col=NA) + 
    geom_line() + 
    ggtitle("Day of ice breakup") + 
    scale_y_continuous("Day of year") + 
    scale_linetype_discrete("Lake data") + 
    theme_classic(9) + theme(legend.position='none', axis.title.x=element_blank(), plot.margin=unit(c(1,0.25, 0.25, 0.25), units="lines"))
  
p_trend_iceout_derivative <- ggplot(derivatives_df%>%filter(metric=="iceout"& data=="all"), 
                                    aes(x=winter.year, y=derivative, lty=data)) + 
    geom_hline(yintercept=0) + 
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2, col=NA) + 
    geom_line() + 
    scale_y_continuous("Days/year", limits=c(-0.43, 0.43)) + 
    scale_x_continuous("Winter Year") + 
    scale_linetype_discrete("Lake data") + 
    theme_classic(9) + theme(legend.position='none', axis.title.x=element_blank(), plot.margin=unit(c(2,0.25, 0.25, 0.25), units="lines"))
  
### iceoutlong

p_trend_iceoutlong_trend <- ggplot(trends_df%>%filter(metric=="iceoutlong"& data=="all"), 
                                   aes(x=winter.year, y=est, lty=data)) + 
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2, col=NA) + 
    geom_line() + 
    ggtitle("Day of ice breakup\n(all data)") + 
    scale_y_continuous("Day of year") + 
    scale_linetype_discrete("Lake data") +  
    theme_classic(9) + theme(legend.position='none', axis.title.x=element_blank(), plot.margin=unit(c(1,0.25, 0.25, 0.25), units="lines"))
  
p_trend_iceoutlong_derivative <- ggplot(derivatives_df%>%filter(metric=="iceoutlong"& data=="all"), 
                                        aes(x=winter.year, y=derivative, lty=data)) + 
    geom_hline(yintercept=0) + 
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2, col=NA) + 
    geom_line() + 
    scale_y_continuous("Days/year", limits=c(-0.43, 0.43)) + 
    scale_x_continuous("Winter Year") + 
    scale_linetype_discrete("Lake data") + 
    theme_classic(9) + theme(legend.position='none', axis.title.x=element_blank(), plot.margin=unit(c(2,0.25, 0.25, 0.25), units="lines"))

### FIGURE 2--Multi-panel plot of trends ----

png("figures/Figure2_Trends_Vertical.png", width=3.5, height=5.5, units='in', res=1200)
grid.arrange(plot_grid(p_trend_duration_trend , p_trend_duration_derivative,
                       p_trend_iceon_trend, p_trend_iceon_derivative,
                       p_trend_iceout_trend, p_trend_iceout_derivative,
                       nrow=3, labels='auto'), bottom="Winter Year")
dev.off()

## Trends plot for supplement ----

### duration

p_supp_trend_duration_trend <- ggplot(trends_df%>%filter(metric=="duration"), 
                                      aes(x=winter.year, y=est, lty=data)) + 
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2, col=NA) + 
  geom_line() + 
  ggtitle("Ice cover duration") + 
  scale_y_continuous("Duration (days)") + 
  scale_linetype_discrete("Lake data") + 
  theme_classic(9) + theme(legend.position=c(0.5, 1.25), 
                           legend.direction='horizontal',
                           axis.title.x=element_blank(),
                           plot.margin=unit(c(2,0.25, 0.25, 0.25), units='lines'))

p_supp_trend_duration_derivative <- ggplot(derivatives_df%>%filter(metric=="duration"), 
                                           aes(x=winter.year, y=derivative, lty=data)) + 
  geom_hline(yintercept=0) + 
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2, col=NA) + 
  geom_line() + 
  scale_y_continuous("Change in duration\n(days/year)", limits=c(-0.43, 0.43)) + 
  scale_x_continuous("Winter Year") + 
  scale_linetype_discrete("Lake data") + 
  theme_classic(9) + theme(legend.position='none', axis.title.x=element_blank(), plot.margin=unit(c(3,0.25, 0.25, 0.25), units="lines"))

### iceon

p_supp_trend_iceon_trend <- ggplot(trends_df%>%filter(metric=="iceon"), 
                                   aes(x=winter.year, y=est, lty=data)) + 
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2, col=NA) + 
  geom_line() + 
  ggtitle("Day of ice formation") + 
  scale_y_continuous("Day of year") + 
  scale_linetype_discrete("Lake data") + 
  theme_classic(9) + theme(legend.position='none', axis.title.x=element_blank(), plot.margin=unit(c(1,0.25, 0.25, 0.25), units="lines"))

p_supp_trend_iceon_derivative <- ggplot(derivatives_df%>%filter(metric=="iceon"), 
                                        aes(x=winter.year, y=derivative, lty=data)) + 
  geom_hline(yintercept=0) + 
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2, col=NA) + 
  geom_line() + 
  scale_y_continuous("Change in day of year\n(days/year)", limits=c(-0.43, 0.43)) + 
  scale_x_continuous("Winter Year") + 
  scale_linetype_discrete("Lake data") + 
  theme_classic(9) + theme(legend.position='none', axis.title.x=element_blank(), plot.margin=unit(c(2,0.25, 0.25, 0.25), units="lines"))

### iceoutlong

p_supp_trend_iceoutlong_trend <- ggplot(trends_df%>%filter(metric=="iceoutlong"), 
                                        aes(x=winter.year, y=est, lty=data)) + 
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2, col=NA) + 
  geom_line() + 
  ggtitle("Day of ice breakup (all data)") + 
  scale_y_continuous("Day of year") + 
  scale_linetype_discrete("Lake data") +  
  theme_classic(9) + theme(legend.position='none', axis.title.x=element_blank(), plot.margin=unit(c(1,0.25, 0.25, 0.25), units="lines"))

p_supp_trend_iceoutlong_derivative <- ggplot(derivatives_df%>%filter(metric=="iceoutlong"), 
                                             aes(x=winter.year, y=derivative, lty=data)) + 
  geom_hline(yintercept=0) + 
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2, col=NA) + 
  geom_line() + 
  scale_y_continuous("Change in day of year\n(days/year)", limits=c(-0.43, 0.43)) + 
  scale_x_continuous("Winter Year") + 
  scale_linetype_discrete("Lake data") + 
  theme_classic(9) + theme(legend.position='none', axis.title.x=element_blank(), plot.margin=unit(c(2,0.25, 0.25, 0.25), units="lines"))

### FIGURE S2--Multi-panel plot of trends with 10y and 30y datasets and larger iceout dataset ----
png("figures/FigureS2_Trends_Vertical_AllBreakup.png", width=6, height=7, units='in', res=1200)
grid.arrange(plot_grid(p_supp_trend_duration_trend , p_supp_trend_duration_derivative,
                       p_supp_trend_iceon_trend, p_supp_trend_iceon_derivative,
                       p_supp_trend_iceoutlong_trend, p_supp_trend_iceoutlong_derivative,
                       nrow=3, labels='auto'), bottom="Winter Year")
dev.off()

# Estimates of change from beginning to end of dataset ----

change_df <- reml_gam_final_models_meta%>%
  mutate(early=NA, early_lower=NA, early_upper=NA,
         late=NA, late_lower=NA, late_upper=NA,
         change=NA, change_lower=NA, change_upper=NA)

set.seed(8675309)

for(i in 1:dim(reml_gam_final_models_meta)[1]){
  
  early_i <- colMeans(simulate(get(reml_gam_final_models_meta$models[i]),data=data.frame(winter.year=1949:1953), 
                             exclude=c("s(fwinter.year)", "s(lnArea_acres)", "s(ID)", "s(x,y)"),
                             newdata.guaranteed=T, nsim=10000)[, 1:10000])
  
  late_i <- colMeans(simulate(get(reml_gam_final_models_meta$models[i]),data=data.frame(winter.year=2018:2022), 
                            exclude=c("s(fwinter.year)", "s(lnArea_acres)", "s(ID)", "s(x,y)"),
                            newdata.guaranteed=T, nsim=10000)[, 1:10000])
  
  change_i <- late_i - early_i
  
  change_df$early[i] <- mean(early_i)
  change_df$early_lower[i] <- quantile(early_i, probs=0.025)
  change_df$early_upper[i] <- quantile(early_i, probs=0.975)
  
  change_df$late[i] <- mean(late_i)
  change_df$late_lower[i] <- quantile(late_i, probs=0.025)
  change_df$late_upper[i] <- quantile(late_i, probs=0.975)
  
  change_df$change[i] <- mean(change_i)
  change_df$change_lower[i] <- quantile(change_i, probs=0.025)
  change_df$change_upper[i] <- quantile(change_i, probs=0.975)
  
  print(i)
  
}

change_df$report <- paste(round(change_df$change, digits=1), " (",
                         round(change_df$change_lower, digits=1), " - ",
                         round(change_df$change_upper, digits=1), ")", sep="")
change_df$reportCentury <- paste(round(change_df$change*(100/69), digits=1), " (",
                          round(change_df$change_lower*(100/69), digits=1), " - ",
                          round(change_df$change_upper*(100/69), digits=1), ")", sep="")

# percent change

change_df$perc_change <- (change_df$change/change_df$early)*100
change_df$perc_change_lower <- (change_df$change_lower/change_df$early)*100
change_df$perc_change_upper <- (change_df$change_upper/change_df$early)*100

## TABLE 1, S2--change estimates----

write.csv(x=change_df, file="tables/change_estimates.csv")

# Lake specific trends ----

reml_gam_final_lakespecific_models_meta

## duration

stab_gamreml_lakespecific_duration <- summary(gamreml_lakespecific_duration)$s.table

lakespecific_duration_df <- expand.grid(winter.year=1970:2022,
                               ID=IDseachDecade$ID)

lakespecific_duration_df <- left_join(lakespecific_duration_df,
                                      ice_w10y1970_condense %>% 
                                        dplyr::select(ID, lnArea_acres, x, y) %>%
                                        group_by(ID) %>%
                                        summarize_all(mean),
                                      by="ID")

lakespecific_duration_df <- lakespecific_duration_df%>%
  mutate(est=NA, se=NA)


lakespecific_duration_df$est <- predict(gamreml_lakespecific_duration, newdata=lakespecific_duration_df[, 1:5], 
                               exclude=c("s(fwinter.year)"),
                               newdata.guaranteed=T)[, 1]

lakespecific_duration_df$se <- predict(gamreml_lakespecific_duration, newdata=lakespecific_duration_df[, 1:5], 
                               exclude=c("s(fwinter.year)"),
                               newdata.guaranteed=T, se.fit=T)$se.fit[, 1]

lakespecific_duration_df$p_value <- NA
for(i in 1:dim(lakespecific_duration_df)[1]){
  
  lakespecific_duration_df$p_value[i] <- stab_gamreml_lakespecific_duration[grepl(lakespecific_duration_df$ID[i], 
                                                                                       rownames(stab_gamreml_lakespecific_duration)), 4]
  
}

lakespecific_duration_df$sig <- ifelse(lakespecific_duration_df$p_value >= 0.1, "p ≥ 0.1", 
                                       ifelse(lakespecific_duration_df$p_value < 0.1 & lakespecific_duration_df$p_value >= 0.05, "0.1 > p ≥ 0.05", 
                                              ifelse(lakespecific_duration_df$p_value < 0.05 & lakespecific_duration_df$p_value >= 0.01, "0.05 > p ≥ 0.01",
                                                     ifelse(lakespecific_duration_df$p_value < 0.01 , "p < 0.01", NA))))

lakespecific_duration_df$fsig <- factor(as.factor(lakespecific_duration_df$sig), levels=c("p ≥ 0.1", "0.1 > p ≥ 0.05", "0.05 > p ≥ 0.01", "p < 0.01"))

p_lakespecific_duration_condensed <- ggplot() + 
  facet_grid(~fsig) + 
  geom_line(data=lakespecific_duration_df, aes(x=winter.year, col=fsig, y=est, lty=ID)) + 
  scale_color_manual(name="Significance of difference from global trend", values=c("black", plasma(4)[1:3])) + 
  scale_y_continuous(" \n") + 
  scale_x_continuous(name="Winter Year", breaks=c(1970, 1990, 2010)) + 
  scale_linetype_manual(values=rep(1, 41)) + 
  ggtitle("a) Ice cover duration") + 
  theme_classic(9) + theme(legend.position='none')

p_lakespecific_duration <- ggplot() + 
  facet_wrap(~ID) + 
  geom_point(data=ice_w10y1970_condense, aes(x=winter.year, y=min_duration), col='grey') + 
  geom_ribbon(data=lakespecific_duration_df, aes(x=winter.year, ymin=est-2*se, ymax=est+2*se, fill=fsig), alpha=0.5) + 
  geom_line(data=lakespecific_duration_df, aes(x=winter.year, col=fsig, y=est)) + 
  scale_color_manual(name="Significance of difference from global trend", values=c("black", plasma(4)[1:3])) + 
  scale_fill_manual(name="Significance of difference from global trend", values=c("black", plasma(4)[1:3])) + 
  scale_y_continuous(name="Predicted ice cover duration trend (annual random effect excluded)") + 
  scale_x_continuous(name="Winter Year", breaks=c(1970, 1990, 2010)) + 
  theme_classic(8) + theme(legend.position='top')

## iceon

stab_gamreml_lakespecific_iceon <- summary(gamreml_lakespecific_iceon)$s.table

lakespecific_iceon_df <- expand.grid(winter.year=1970:2022,
                                        ID=IDseachDecade$ID)

lakespecific_iceon_df <- left_join(lakespecific_iceon_df,
                                      ice_w10y1970_condense %>% 
                                        dplyr::select(ID, lnArea_acres, x, y) %>%
                                        group_by(ID) %>%
                                        summarize_all(mean),
                                      by="ID")

lakespecific_iceon_df <- lakespecific_iceon_df%>%
  mutate(est=NA, se=NA)


lakespecific_iceon_df$est <- predict(gamreml_lakespecific_iceon, newdata=lakespecific_iceon_df[, 1:5], 
                                        exclude=c("s(fwinter.year)"),
                                        newdata.guaranteed=T)[,1]

lakespecific_iceon_df$se <- predict(gamreml_lakespecific_iceon, newdata=lakespecific_iceon_df[, 1:5], 
                                       exclude=c("s(fwinter.year)"),
                                       newdata.guaranteed=T, se.fit=T)$se.fit[,1]

lakespecific_iceon_df$p_value <- NA
for(i in 1:dim(lakespecific_iceon_df)[1]){
  
  lakespecific_iceon_df$p_value[i] <- stab_gamreml_lakespecific_iceon[grepl(lakespecific_iceon_df$ID[i], 
                                                                                       rownames(stab_gamreml_lakespecific_iceon)), 4]
  
}

lakespecific_iceon_df$sig <- ifelse(lakespecific_iceon_df$p_value >= 0.1, "p ≥ 0.1", 
                                       ifelse(lakespecific_iceon_df$p_value < 0.1 & lakespecific_iceon_df$p_value >= 0.05, "0.1 > p ≥ 0.05", 
                                              ifelse(lakespecific_iceon_df$p_value < 0.05 & lakespecific_iceon_df$p_value >= 0.01, "0.05 > p ≥ 0.01",
                                                     ifelse(lakespecific_iceon_df$p_value < 0.01, "p < 0.01", NA))))

lakespecific_iceon_df$fsig <- factor(as.factor(lakespecific_iceon_df$sig), levels=c("p ≥ 0.1", "0.1 > p ≥ 0.05", "0.05 > p ≥ 0.01", "p < 0.01"))

p_lakespecific_iceon_condensed <- ggplot() + 
  facet_grid(~fsig) + 
  geom_line(data=lakespecific_iceon_df, aes(x=winter.year, col=fsig, y=est, lty=ID)) + 
  scale_color_manual(name="Significance of difference from global trend", values=c("black", plasma(4)[1:3])) + 
  scale_y_continuous("Predicted trend\n(annual random effect excluded)") + 
  scale_x_continuous(name="Winter Year", breaks=c(1970, 1990, 2010)) + 
  scale_linetype_manual(values=rep(1, 41)) + 
  ggtitle("b) Day of ice formation") + 
  theme_classic(9) + theme(legend.position='none')


## iceout

stab_gamreml_lakespecific_iceout <- summary(gamreml_lakespecific_iceout)$s.table

lakespecific_iceout_df <- expand.grid(winter.year=1970:2022,
                                     ID=IDseachDecade$ID)

lakespecific_iceout_df <- left_join(lakespecific_iceout_df,
                                   ice_w10y1970_condense %>% 
                                     dplyr::select(ID, lnArea_acres, x, y) %>%
                                     group_by(ID) %>%
                                     summarize_all(mean),
                                   by="ID")

lakespecific_iceout_df <- lakespecific_iceout_df%>%
  mutate(est=NA, se=NA)


lakespecific_iceout_df$est <- predict(gamreml_lakespecific_iceout, newdata=lakespecific_iceout_df[, 1:5], 
                                     exclude=c("s(fwinter.year)"),
                                     newdata.guaranteed=T)[,1]

lakespecific_iceout_df$se <- predict(gamreml_lakespecific_iceout, newdata=lakespecific_iceout_df[, 1:5], 
                                    exclude=c("s(fwinter.year)"),
                                    newdata.guaranteed=T, se.fit=T)$se.fit[,1]

lakespecific_iceout_df$p_value <- NA
for(i in 1:dim(lakespecific_iceout_df)[1]){
  
  lakespecific_iceout_df$p_value[i] <- stab_gamreml_lakespecific_iceout[grepl(lakespecific_iceout_df$ID[i], 
                                                                                 rownames(stab_gamreml_lakespecific_iceout)), 4]
  
}

lakespecific_iceout_df$sig <- ifelse(lakespecific_iceout_df$p_value >= 0.1, "p ≥ 0.1", 
                                    ifelse(lakespecific_iceout_df$p_value < 0.1 & lakespecific_iceout_df$p_value >= 0.05, "0.1 > p ≥ 0.05", 
                                           ifelse(lakespecific_iceout_df$p_value < 0.05 & lakespecific_iceout_df$p_value >= 0.01, "0.05 > p ≥ 0.01",
                                                  ifelse(lakespecific_iceout_df$p_value < 0.01, "p < 0.01", NA))))

lakespecific_iceout_df$fsig <- factor(as.factor(lakespecific_iceout_df$sig), levels=c("p ≥ 0.1", "0.1 > p ≥ 0.05", "0.05 > p ≥ 0.01", "p < 0.01"))

p_lakespecific_iceout_condensed <- ggplot() + 
  facet_grid(~fsig) + 
  geom_line(data=lakespecific_iceout_df, aes(x=winter.year, col=fsig, y=est, lty=ID)) + 
  scale_color_manual(name="Significance of difference from global trend", values=c("black", plasma(4)[1:3])) + 
  scale_y_continuous(name=" \n ") + 
  scale_x_continuous(name="Winter Year", breaks=c(1970, 1990, 2010)) + 
  scale_linetype_manual(values=rep(1, 41)) + 
  ggtitle("c) Day of ice breakup") + 
  theme_classic(9) + theme(legend.position='none')

## FIGURE S3--Lake specific trends for all three metrics ----

p_lakespecific_condensed <- plot_grid(
  
  p_lakespecific_duration_condensed, 
  p_lakespecific_iceon_condensed, 
  p_lakespecific_iceout_condensed, 
  
  nrow=3)

png("figures/FigureS3_LakeSpecificTrends_Condensed.png", width=7, height=7, units='in', res=1200)
p_lakespecific_condensed
dev.off()

# Maps and spatial analysis ----

# Predictions
# Here, we're looking at just the long-term trend, so fwinter.year is removed

## Predictions for average trend from 1949 to 2022 across all lakes ----

### duration
new.data_duration <- expand.grid(winter.year=c(1949, 2022), 
                                           ID=unique(gamreml_all_duration$model$ID))

new.data_duration <- left_join(new.data_duration, 
                                         gamreml_all_duration$model %>%
                                           dplyr::select(ID, x, y, lnArea_acres) %>%
                                           group_by(ID) %>%
                                           summarize_all(mean),
                                         by="ID")

pred_duration <- predict.gam(gamreml_all_duration, 
                                       newdata=new.data_duration, 
                                       exclude=c("s(fwinter.year)"),
                                       se.fit=TRUE, newdata.guaranteed=TRUE)
new.data_duration$pred <- pred_duration$fit[,1]
new.data_duration$se <- pred_duration$se.fit[,1]

### iceon
new.data_iceon <- expand.grid(winter.year=c(1949, 2022), 
                                           ID=unique(gamreml_all_iceon$model$ID))

new.data_iceon <- left_join(new.data_iceon, 
                                         gamreml_all_iceon$model %>%
                                           dplyr::select(ID, x, y, lnArea_acres) %>%
                                           group_by(ID) %>%
                                           summarize_all(mean),
                                         by="ID")

pred_iceon <- predict.gam(gamreml_all_iceon, 
                                       newdata=new.data_iceon, 
                                       exclude=c("s(fwinter.year)"),
                                       se.fit=TRUE, newdata.guaranteed=TRUE)
new.data_iceon$pred <- pred_iceon$fit[,1]
new.data_iceon$se <- pred_iceon$se.fit[,1]

### iceout
new.data_iceout <- expand.grid(winter.year=c(1949, 2022), 
                                           ID=unique(gamreml_all_iceout$model$ID))

new.data_iceout <- left_join(new.data_iceout, 
                                         gamreml_all_iceout$model %>%
                                           dplyr::select(ID, x, y, lnArea_acres) %>%
                                           group_by(ID) %>%
                                           summarize_all(mean),
                                         by="ID")

pred_iceout <- predict.gam(gamreml_all_iceout, 
                                       newdata=new.data_iceout, 
                                       exclude=c("s(fwinter.year)"),
                                       se.fit=TRUE, newdata.guaranteed=TRUE)
new.data_iceout$pred <- pred_iceout$fit[,1]
new.data_iceout$se <- pred_iceout$se.fit[,1]

### iceoutlong
new.data_iceoutlong <- expand.grid(winter.year=c(1949, 2022), 
                                           ID=unique(gamreml_all_iceoutlong$model$ID))

new.data_iceoutlong <- left_join(new.data_iceoutlong, 
                                         gamreml_all_iceoutlong$model %>%
                                           dplyr::select(ID, x, y, lnArea_acres) %>%
                                           group_by(ID) %>%
                                           summarize_all(mean),
                                         by="ID")

pred_iceoutlong <- predict.gam(gamreml_all_iceoutlong, 
                                       newdata=new.data_iceoutlong, 
                                       exclude=c("s(fwinter.year)"),
                                       se.fit=TRUE, newdata.guaranteed=TRUE)
new.data_iceoutlong$pred <- pred_iceoutlong$fit[, 1]
new.data_iceoutlong$se <- pred_iceoutlong$se.fit[, 1]

### combine predictions to single data frame
new.data_duration$variable <- "Duration"
new.data_iceon$variable <- "Ice formation"
new.data_iceout$variable <- "Ice breakup"
new.data_iceoutlong$variable <- "Ice breakup (all data)"

new.data <- rbind(new.data_duration, 
                            new.data_iceon, 
                            new.data_iceout, 
                            new.data_iceoutlong)

head(new.data)

## Make the predictions spatial ----

MN.shape_agg_centroids <- st_centroid(MN.shape_agg)

MN.shape_agg_centroids$X <- st_coordinates(MN.shape_agg_centroids$geometry)[, "X"]
MN.shape_agg_centroids$Y <- st_coordinates(MN.shape_agg_centroids$geometry)[, "Y"]

new.data_spatial <- inner_join(new.data, MN.shape_agg_centroids) %>% 
  rename(Area_acres=AggDOW_Area_acres) %>% 
  dplyr::select(-dowlknum)

new.data_spatial_sf <- st_as_sf(new.data_spatial)
new.data_spatial_sf <- st_transform(new.data_spatial_sf, crs = st_crs(MN.shape_agg))

head(new.data_spatial_sf)

## Spatial analysis for GAM predictions ----

vg_predicted_duration <- Variogram(object = new.data_spatial_sf$pred[new.data_spatial_sf$winter.year == 2022 & new.data_spatial_sf$variable == "Duration"],
                                   distance = dist(data.frame(X = new.data_spatial_sf$X[new.data_spatial_sf$winter.year == 2022 & new.data_spatial_sf$variable == "Duration"],
                                                              Y = new.data_spatial_sf$Y[new.data_spatial_sf$winter.year == 2022 & new.data_spatial_sf$variable == "Duration"])))

vg_predicted_iceon <- Variogram(object = new.data_spatial_sf$pred[new.data_spatial_sf$winter.year == 2022 & new.data_spatial_sf$variable == "Ice formation"],
                                distance = dist(data.frame(X = new.data_spatial_sf$X[new.data_spatial_sf$winter.year == 2022 & new.data_spatial_sf$variable == "Ice formation"],
                                                           Y = new.data_spatial_sf$Y[new.data_spatial_sf$winter.year == 2022 & new.data_spatial_sf$variable == "Ice formation"])))

vg_predicted_iceoff <- Variogram(object = new.data_spatial_sf$pred[new.data_spatial_sf$winter.year == 2022 & new.data_spatial_sf$variable == "Ice breakup"],
                                 distance = dist(data.frame(X = new.data_spatial_sf$X[new.data_spatial_sf$winter.year == 2022 & new.data_spatial_sf$variable == "Ice breakup"],
                                                            Y = new.data_spatial_sf$Y[new.data_spatial_sf$winter.year == 2022 & new.data_spatial_sf$variable == "Ice breakup"])))

vg_predicted_duration$variable <- "Duration"
vg_predicted_iceon$variable <- "Ice formation"
vg_predicted_iceoff$variable <- "Ice breakup"

vg_predicted_allvars <- rbind(vg_predicted_duration, vg_predicted_iceon, vg_predicted_iceoff) %>% mutate(variable = factor(variable, levels = c("Duration", "Ice formation", "Ice breakup")))

p_predicted_variogram_duration <- ggplot(vg_predicted_allvars %>% filter(variable == "Duration"), 
                                         aes(x = dist/1000, y = variog)) + 
  facet_wrap(~variable, nrow = 1) + 
  geom_point(cex = 0.8, alpha = 0.5, col = 'grey') + 
  geom_smooth(method = 'gam', formula = y~s(x), col = 'black', fill = 'black', alpha = 0.3, lwd = 1) + 
  scale_y_continuous(name = expression(gamma)) + 
  scale_x_continuous(name = "Distance (km)") + 
  ggtitle("Semivariograms for each ice variable prediction") + 
  theme_cowplot(8) + theme(legend.position = 'none', axis.title.x = element_blank(),
                           plot.title = element_text(size = 8, face = "plain"))

p_predicted_variogram_onoff <- ggplot(vg_predicted_allvars %>% filter(variable != "Duration"), 
                                      aes(x = dist/1000, y = variog)) + 
  facet_wrap(~variable, nrow = 1) + 
  geom_point(cex = 0.8, alpha = 0.5, col = 'grey') + 
  geom_smooth(method = 'gam', formula = y~s(x), col = 'black', fill = 'black', alpha = 0.3, lwd = 1) + 
  scale_y_continuous(name = expression(gamma)) + 
  scale_x_continuous(name = "Distance (km)") + 
  theme_cowplot(8)

p_predicted_variogram_allvars <- plot_grid(p_predicted_variogram_duration, 
                                           p_predicted_variogram_onoff, 
                                           ncol = 1)

## FIGURE 3--Map + semivariogram plot ----

### Maps ----
p_pred_duration <- ggplot(data = new.data_spatial_sf %>% filter(variable == "Duration"), 
                                    aes(col = pred)) + 
  facet_grid(~winter.year) + 
  geom_sf(data = MN_outline, 
          fill = 'grey', col = NA) + 
  geom_sf(cex = 0.8) + 
  scale_color_gradient2(breaks = c(100, 120, 140, 160, 180),
                        low = "red", mid = "white", high = "blue",
                        midpoint = median(new.data_duration$pred)) +
  scale_x_continuous(breaks = c(-96, -93, -90)) + 
  scale_y_continuous(breaks = c(44, 46, 48)) +
  ggtitle("Ice cover duration (days)") + 
  theme_cowplot(9) + theme(plot.title = element_text(face = "plain"),
                           legend.title = element_blank(),
                           legend.key.height = unit(0.04, units = 'npc'),
                           axis.title = element_blank())

p_pred_iceon <- ggplot(data = new.data_spatial_sf %>% filter(variable == "Ice formation"), 
                                 aes(col = pred)) + 
  facet_grid(~winter.year) + 
  geom_sf(data = MN_outline, 
          fill = 'grey', col = NA) + 
  geom_sf(cex = 0.8) + 
  scale_color_gradient2(breaks = seq(312, 361, 7),
                        labels = c("Nov 8", "Nov 15", "Nov 22", "Nov 29", "Dec 6", "Dec 13", "Dec 20", "Dec 27"),
                        low = "blue", mid = "white", high = "red",
                        midpoint = median(new.data_iceon$pred)) +
  scale_x_continuous(breaks = c(-96, -93, -90)) + 
  scale_y_continuous(breaks = c(44, 46, 48)) +
  ggtitle("Day of ice formation") + 
  theme_cowplot(9) + theme(plot.title = element_text(face = "plain"),
                           legend.title = element_blank(),
                           legend.key.height = unit(0.04, units = 'npc'),
                           axis.title = element_blank())

p_pred_iceout <- ggplot(data = new.data_spatial_sf %>% filter(variable == "Ice breakup"), 
                                  aes(col = pred)) + 
  facet_grid(~winter.year) + 
  geom_sf(data = MN_outline, 
          fill = 'grey', col = NA) + 
  geom_sf(cex = 0.8) + 
  scale_color_gradient2(breaks = seq(83, 132, 7),
                        labels = c("Mar 24", "Mar 31", "Apr 7", "Apr 14", "Apr 21", "Apr 28", "May 5", "May 12"),
                        low = "red", mid = "white", high = "blue",
                        midpoint = median(new.data_iceout$pred)) +
  scale_x_continuous(breaks = c(-96, -93, -90)) + 
  scale_y_continuous(breaks = c(44, 46, 48)) +
  ggtitle("Day of ice breakup") + 
  theme_cowplot(9) + theme(plot.title = element_text(face = "plain"),
                           legend.title = element_blank(),
                           legend.key.height = unit(0.04, units = 'npc'),
                           axis.title = element_blank())


### Variograms ----

p_predicted_variogram_duration2 <- ggplot(vg_predicted_allvars %>% filter(variable == "Duration"), 
                                          aes(x = dist/1000, y = variog)) + 
  geom_point(cex = 0.8, alpha = 0.5, col = 'grey') + 
  geom_smooth(method = 'gam', formula = y~s(x), col = 'black', fill = 'black', alpha = 0.3, lwd = 1) + 
  scale_y_continuous(name = expression(gamma)) + 
  scale_x_continuous(name = "") + 
  theme_cowplot(8) + theme(legend.position = 'none',
                           plot.title = element_text(size = 8, face = "plain"))

p_predicted_variogram_iceon <- ggplot(vg_predicted_allvars %>% filter(variable == "Ice formation"), 
                                   aes(x = dist/1000, y = variog)) + 
  geom_point(cex = 0.8, alpha = 0.5, col = 'grey') + 
  geom_smooth(method = 'gam', formula = y~s(x), col = 'black', fill = 'black', alpha = 0.3, lwd = 1) + 
  scale_y_continuous(name = expression(gamma)) + 
  scale_x_continuous(name = "") + 
  theme_cowplot(8) + theme(legend.position = 'none',
                           plot.title = element_text(size = 8, face = "plain"))

p_predicted_variogram_iceout <- ggplot(vg_predicted_allvars %>% filter(variable == "Ice breakup"), 
                                    aes(x = dist/1000, y = variog)) + 
  geom_point(cex = 0.8, alpha = 0.5, col = 'grey') + 
  geom_smooth(method = 'gam', formula = y~s(x), col = 'black', fill = 'black', alpha = 0.3, lwd = 1) + 
  scale_y_continuous(name = expression(gamma)) + 
  scale_x_continuous(name = "Distance (km)") + 
  theme_cowplot(8) + theme(legend.position = 'none',
                           plot.title = element_text(size = 8, face = "plain"))

p_mapvariogram <- plot_grid(p_pred_duration, p_predicted_variogram_duration2 + ggtitle(""),
          p_pred_iceon, p_predicted_variogram_iceon + ggtitle(""),
          p_pred_iceout, p_predicted_variogram_iceout + ggtitle(""),
          ncol = 2, rel_widths = c(1, 0.5, 1, 0.5, 1, 0.5), 
          labels = "auto", label_size = 8)

png("figures/Figure3_Map&Semivariogram.png", width = 5, height = 5.5, units='in', res=1200)
p_mapvariogram
dev.off()

# Sensitivity to time frame of analysis ----

## set up windows (time frames) ----
window_df2 <- data.frame(length=74, start=1949, end=2022)

for(i in c(seq(10, 70, 10), 74)){
  
  for(j in seq(1949, 2014, 1)){
    
    if(i < 74){
      
      evenly_i <- evenly(j:2022, by=i-1)
      
      window_start_i <- evenly_i[-length(evenly_i)]
      window_end_i <- evenly_i[-1]
      
      window_df2_i <- data.frame(length=rep(i, length(window_start_i)),
                                 start=window_start_i,
                                 end=window_end_i)
      
      window_df2 <- rbind(window_df2, window_df2_i)
      
    }
    
  }
  
}

window_df2 <- window_df2 %>%
  arrange(length, start)

# I'm not sure why I'm getting duplicates here, I'm sure 
# my code above is wrong, but this gets us what we want...
window_df2 <- window_df2 %>% distinct()


## fit models and pull trends and change estimates ----

window_fit_df2 <- data.frame(length=numeric(0),
                             start=numeric(0),
                             end=numeric(0),
                             winter.year=numeric(0),
                             est=numeric(0),
                             lower=numeric(0),
                             upper=numeric(0))

window_fit_df2_change <- data.frame(length=numeric(0),
                                    start=numeric(0),
                                    end=numeric(0),
                                    change=numeric(0),
                                    change_lower=numeric(0),
                                    change_upper=numeric(0),
                                    change_century=numeric(0),
                                    change_century_lower=numeric(0),
                                    change_century_upper=numeric(0))

for(i in 1:dim(window_df2)[1]){
#for(i in 1){
  
  gam_fwinter.year_i <- bam(min_duration ~ s(winter.year, k=4) + s(fwinter.year, bs='re') + s(lnArea_acres, k=5) + s(ID, bs='re') + s(x, y, k=5),
                            data = ice_all_condense %>% filter(winter.year >= window_df2$start[i] & 
                                                                  winter.year <= window_df2$end[i]),
                            nthreads = detectCores() - 1,
                            #method='REML', 
                            select=TRUE)
  
  # Estimating Trends
  year_seq_i <- window_df2$start[i]:window_df2$end[i]
  # Intercept included in trends
  confint_fwinter.year_i <- confint(object=gam_fwinter.year_i, parm="s(winter.year)", n=length(year_seq_i), type='simultaneous', shift=TRUE)
  
  window_fit_df2_i <- data.frame(length=window_df2$length[i],
                                 start=window_df2$start[i],
                                 end=window_df2$end[i],
                                 winter.year=year_seq_i,
                                 est=confint_fwinter.year_i$est,
                                 lower=confint_fwinter.year_i$lower,
                                 upper=confint_fwinter.year_i$upper)
  
  window_fit_df2 <- rbind(window_fit_df2, window_fit_df2_i)
  
  # Estimating change
  early_i <- data.frame(winter.year=year_seq_i[year_seq_i<=window_df2$start[i]+4])
  late_i <- data.frame(winter.year=year_seq_i[year_seq_i>=window_df2$end[i]-4])
  
  pred_early <- predicted_samples(gam_fwinter.year_i, n=1000, data=early_i, exclude=c("s(fwinter.year)", "s(lnArea_acres)", "s(ID)", "s(x,y)"), newdata.guaranteed=TRUE)
  pred_late <- predicted_samples(gam_fwinter.year_i, n=1000, data=late_i, exclude=c("s(fwinter.year)", "s(lnArea_acres)", "s(ID)", "s(x,y)"), newdata.guaranteed=TRUE)
  
  pred_early <- pred_early %>%
    group_by(draw) %>%
    summarize(response=mean(response))
  
  pred_late <- pred_late %>%
    group_by(draw) %>%
    summarize(response=mean(response))
  
  pred_change <- pred_late$response - pred_early$response
  dt <- window_df2$end[i] - window_df2$start[i] - 5
  
  window_fit_df2_change_i <- data.frame(length=window_df2$length[i],
                                 start=window_df2$start[i],
                                 end=window_df2$end[i],
                                 change=mean(pred_change),
                                 change_lower=quantile(pred_change, probs=0.025),
                                 change_upper=quantile(pred_change, probs=0.975),
                                 change_century=(mean(pred_change)/(dt))*100,
                                 change_century_lower=(quantile(pred_change, probs=0.025)/(dt))*100,
                                 change_century_upper=(quantile(pred_change, probs=0.975)/(dt))*100)
  
  window_fit_df2_change <- rbind(window_fit_df2_change, window_fit_df2_change_i)
  
  print(window_df2[i, ])
  
}

window_fit_df2_change <- window_fit_df2_change %>%
  mutate(Sig=ifelse(change_lower*change_upper <=0, "NS", "Sig."),
         Dir=ifelse(change < 0, "Loss", "Gain"))

summary(factor(paste(window_fit_df2_change$Sig, window_fit_df2_change$Dir)))
summary(factor(paste(window_fit_df2$Sig, window_fit_df2$Dir)))

# add info about significant change to window_fit_df2
window_fit_df2 <- left_join(window_fit_df2,
                            window_fit_df2_change %>%
                              dplyr::select(length, start, end, Sig, Dir),
                            by=join_by(length, start, end))

# summarize change by window length
window_fit_df2_change_sum <- window_fit_df2_change %>%
  group_by(length) %>%
  summarize(mean=mean(change_century),
            sd=sd(change_century),
            lower=quantile(change_century, probs=0.025),
            upper=quantile(change_century, probs=0.975),
            min=min(change_century),
            max=max(change_century),
            n=n()) %>%
  mutate(n=n/length) %>%
  ungroup()



## plot trends and change estimates ----

p_windows2_trends <- grid.arrange(ggplot(data=window_fit_df2, 
                                         aes(x=winter.year, y=est, col=Dir, alpha=Sig, lty=factor(start))) + 
                                    facet_grid(length~.) + 
                                    geom_line() + 
                                    scale_color_manual(values=c("dodgerblue", "red4")) + 
                                    scale_alpha_manual(values=c(0.2, 1)) +
                                    scale_linetype_manual(values=rep(1, length(unique(window_fit_df2$start)))) + 
                                    scale_y_continuous("Estimated trend over time frame") + 
                                    scale_x_continuous("Winter Year") + 
                                    theme_bw(8) + theme(legend.position='none'))

p_windows2_change <- grid.arrange(ggplot(data=window_fit_df2_change, aes(change_century)) + 
                                    facet_grid(length ~ ., scales='free') + 
                                    geom_vline(xintercept=0) + 
                                    geom_histogram(fill='skyblue') + 
                                    geom_point(data=window_fit_df2_change_sum,
                                               aes(x=mean, y=0),
                                               inherit.aes=FALSE,
                                               col="orange", cex=1.5) + 
                                    geom_errorbarh(data=window_fit_df2_change_sum,
                                                   aes(y=0, xmin=mean-sd, xmax=mean+sd),
                                                   inherit.aes=FALSE,
                                                   col='orange', height=0) + 
                                    scale_x_continuous("Estimated change over time frame (days/century)") + 
                                    scale_y_continuous("Number of estimates", breaks=pretty_breaks(n=3)) +
                                    theme_bw(8),
                                  right=textGrob("Time frame duration (years)", rot=-90, gp=gpar(fontsize=8)))

png("figures/Fig_WindowStudy2.png", width=4, height=4, units='in', res=600)
plot_grid(p_windows2_trends, p_windows2_change,
          nrow=1, labels="auto", label_size=9)
dev.off()

p_windows2_biplots <- grid.arrange(ggplot(data=window_fit_df2_change, aes(x=end, y=change_century, pch=Sig, alpha=Sig, col=Dir)) + 
                                     facet_grid(length ~ ., scales="free_y") + 
                                     geom_hline(yintercept=0) + 
                                     geom_hline(yintercept=window_fit_df2_change$change_century[window_fit_df2_change$length==74], lty=3) + 
                                     #geom_errorbar(aes(ymin=change_century_lower, ymax=change_century_upper), width=0, lwd=1) + 
                                     geom_point(cex=0.7) + 
                                     scale_shape_manual(values=c(1, 19)) + 
                                     scale_color_manual(values=c("dodgerblue", "red4")) +
                                     scale_alpha_manual(values=c(0.2, 1)) +
                                     scale_x_continuous("Last year in time frame") + 
                                     scale_y_continuous("Estimated change over time frame (days/century)", breaks=pretty_breaks(n=3)) +
                                     theme_bw(8) + theme(legend.position='none'),
                                   right=textGrob("Time frame duration (years)", rot=-90, gp=gpar(fontsize=8)))
p_windows2_biplots_werror <- grid.arrange(ggplot(data=window_fit_df2_change, aes(x=end, y=change_century, pch=Sig, alpha=Sig, col=Dir)) + 
                                     facet_grid(length ~ ., scales="free_y") + 
                                     geom_hline(yintercept=0) + 
                                     geom_hline(yintercept=window_fit_df2_change$change_century[window_fit_df2_change$length==74], lty=3) + 
                                     geom_errorbar(aes(ymin=change_century_lower, ymax=change_century_upper), width=0, lwd=0.7) + 
                                     geom_point(cex=0.7, pch=21, color='black') + 
                                     scale_shape_manual(values=c(1, 19)) + 
                                     scale_color_manual(values=c("dodgerblue", "red4")) +
                                     scale_alpha_manual(values=c(0.2, 1)) +
                                     scale_x_continuous("Last year in time frame") + 
                                     scale_y_continuous("Estimated change over time frame (days/century)", breaks=pretty_breaks(n=3)) +
                                     theme_bw(8) + theme(legend.position='none'),
                                   right=textGrob("Time frame duration (years)", rot=-90, gp=gpar(fontsize=8)))

## FIGURE 4--Sensitivty to time frame of analysis ----

png("figures/Figure4_WindowStudy2_Trend&Biplot_wError.png", width=3.5, height=3.5, units='in', res=1200)
plot_grid(p_windows2_trends, p_windows2_biplots_werror,
          nrow=1, labels="auto", label_size=9)
dev.off()