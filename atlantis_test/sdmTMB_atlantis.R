require(ggplot2)
require(dplyr)
require(sdmTMB)
require(sf)
require(gridExtra)

##### explore basic functionality of sdmTMB for arrowtooth flounder biomass estimates in the GOA from RACE bottom trawl survey: 
#### build models, predict onto grid, and calculate biomass index

##### this data set has been processed to contain zeros for stations where hauls had no biomass recorded for this species, and only contains data from May-August in RACE survey years
atf <- read.csv("atf_clean.csv")

##### take a quick look at the data spatially (plot takes a while to load)
load('goa_coast.RData')

ggplot()+
  geom_point(data = atf, aes(lon, lat, colour = log1p(biom_kgkm2)), size = 1.5)+
  scale_colour_viridis_c()+
  geom_polygon(data = coast, aes(long, lat, group = group), colour = "black", fill = "grey80")+
  theme_minimal()+
  facet_wrap(~year)

#### scale the depth and bottom temp variables for numerical efficiency
atf$depth_scale <- scale(atf$depth)
atf$btemp_scale <- scale(atf$btemp)
atf$stemp_scale <- scale(atf$stemp)

### this is the spde mesh that the sdmTMB algorithm uses to estimate spatial autocorrelation
##### the speed of model running is highly dependent on number of knots you use to build the mesh. 100 is quite low, you'd want to make sure it was robust by checking multiple resolutions when you have a model you want to actually use (people use like 350-450 or something)
atf_spde <- make_mesh(atf, c("lon", "lat"), n_knots = 100)
plot(atf_spde)

##### check out the distribution of the biomass density response variable. Obviously very right-skewed
hist(atf$biom_kgkm2, breaks = 30)
hist(log(atf$biom_kgkm2), breaks = 30)

#### how many zeros in the data
table(atf$biom_kgkm2 == 0)

##### run a model with just space and time using the Tweedie distribution. This explicitly accounts for zeros in the data within a single distribution framework
#### the 0 in the formula is to predict the mean estimate for each time period if you want to do index standardisation. Not necessary for Atlantis giving we're interested in biomass distribution but I've included it here. 
#### note fancy options like anisotropy, spatial trend are turned off here (but see sdmTMB documentation for full options).
#### increasing the loops and steps arguments increases run time but can help with convergence issues 
start.time <- Sys.time()
m_atf <- sdmTMB(
  data = atf, 
  formula = biom_kgkm2 ~ 0 + as.factor(year),
  time = "year", 
  spde = atf_spde, 
  reml = TRUE,
  anisotropy = FALSE,
  spatial_trend = FALSE, 
  spatial_only = FALSE,
  silent = FALSE,
  control = sdmTMBcontrol(),
  nlminb_loops = 3,
  newton_steps = 10,
  family = tweedie(link = "log"))
end.time <- Sys.time()
time.taken_m_atf <- end.time - start.time
time.taken_m_atf

#### check out model residuals - these look good
atf$resids <- residuals(m_atf) # randomized quantile residuals, should look normal
hist(atf$resids)

### could be much worse 
qqnorm(atf$resids)
abline(a = 0, b = 1)

##### try running a slightly more complex model with smooth term for depth (sdmTMB runs GAMs in mgcv under the hood for the covariates in addition to the spatiotemporal model in TMB). You can also use linear terms if you want. 

#### the s() terms are smooths, k is the smoothing parameter that controls the number of "knots" in the splines and therefore the wiggliness of the response curves. I would go lower (3-5) if youre trying to understand general responses to environment and use the models for e.g. forecasting out of sample, maybe higher if you only care about making accurate predictions back onto your data, but good practice to check sensitivity of predictions to these values and optimise depending on model purpose. See mgcv for more. 

start.time <- Sys.time()
m_atf_depth <- sdmTMB(
  data = atf, 
  formula = biom_kgkm2 ~ 0 + s(depth_scale, k = 5) + as.factor(year),
  time = "year", 
  spde = atf_spde, 
  reml = TRUE,
  anisotropy = FALSE,
  spatial_trend = FALSE, 
  spatial_only = FALSE,
  silent = FALSE,
  control = sdmTMBcontrol(),
  nlminb_loops = 3,
  newton_steps = 10,
  family = tweedie(link = "log"))
end.time <- Sys.time()
time.taken_m_atf_depth <- end.time - start.time
time.taken_m_atf_depth

#### check out model residuals - look good
atf$resids <- residuals(m_atf_depth) 
hist(atf$resids)

##### better fit than the previous model
qqnorm(atf$resids)
abline(a = 0, b = 1)

#### you can plot the response curve from the depth smooth term in the mgcv model - arrowtooth do not like the deepest waters  
plot(m_atf_depth$mgcv_mod, rug = TRUE)

#### you can also interrogate the GAM to get the standard summary output (this is only the covariate part so not describing the performance of the model as a whole, but maybe useful to check contribution of covariates etc)
summary(m_atf_depth$mgcv_mod)

##### try an even more complex model with climate covariates (bottom temperature and surface temperature) 
#### the package doesn't seem to have a built in way of dealing with missing values (there are NAs in the temperature columns). Need to build a new mesh that's the correct length to account for for the missing values and run the model on the data with NAs removed (using data = na.exclude(atf)), just in order to inspect the model residuals 
atf_spde_na <- make_mesh(na.exclude(atf), c("lon", "lat"), n_knots = 100)

start.time <- Sys.time()
m_atf_clim <- sdmTMB(
  data = na.exclude(atf), 
  formula = biom_kgkm2 ~ 0 + s(depth_scale, k = 5) + s(stemp_scale, k = 5) + s(btemp_scale, k = 5) + as.factor(year),
  time = "year", 
  spde = atf_spde_na, 
  reml = TRUE,
  anisotropy = FALSE,
  spatial_trend = FALSE, 
  spatial_only = FALSE,
  silent = FALSE,
  control = sdmTMBcontrol(),
  nlminb_loops = 3,
  newton_steps = 10,
  family = tweedie(link = "log"))
end.time <- Sys.time()
time.taken_m_atf_clim <- end.time - start.time
time.taken_m_atf_clim    

resids <- residuals(m_atf_clim)
hist(resids)

# this is looking even better!
qqnorm(resids)
abline(a = 0, b = 1)

#### plot the response curves 
#### temperature partial effects not as strong as depth, seems like arrowtooth maybe prefer intermediate surface temps and maybe more extreme bottom temps
plot(m_atf_clim$mgcv_mod, rug = TRUE)

##### which of these models performed better according to AIC?
#### looks like the one with the climate covariates did substantially better (comes out the same if you run the climate model with NAs), followed by the one with depth and then the one with just space time. Let's see what they look like spatially.
AIC(m_atf, m_atf_depth, m_atf_clim)

#### read in Atlantis prediction grids modified in script Atlantis_grid_covars.R to add environmental covariates 
load("atlantis_grid_bgm.RData")
atlantis_grid_df <- read.csv('atlantis_grid_df.csv')

### scale covariates to match what we did for the arrowtooth data
atlantis_grid_df$depth_scale <- scale(atlantis_grid_df$depth)
atlantis_grid_df$btemp_scale <- scale(atlantis_grid_df$btemp)
atlantis_grid_df$stemp_scale <- scale(atlantis_grid_df$stemp)

#### make SDM predictions onto new data from climate atf model
#### I am using the dataframe version to predict onto because sdmTMB doesn't seem to like the format of the other grid
#### predict for all years
predictions_atf_clim <- predict(m_atf_clim, newdata = atlantis_grid_df, return_tmb_object = TRUE)
atlantis_grid_df$atf_clim <- predictions_atf_clim$data$est

#### add 1990 predictions to Atlantis polygons for plotting 
atlantis_grid_bgm$atf_clim <- subset(atlantis_grid_df$atf_clim, atlantis_grid_df$year == 1990)

#### not plotting most of Canada because the predictions look terrible (due to not having temperature or biomass data from there in this model)
coast <- st_as_sf(coast)
coast_utm <- st_transform(coast, crs = st_crs(atlantis_grid_bgm))

p1 <- ggplot(subset(atlantis_grid_bgm, box_id < 92))+
  geom_sf(aes(fill=atf_clim))+
  geom_sf(data = coast_utm, colour = "black", fill = "grey80")+
  scale_fill_viridis_c(limits = c(-10.25, 9.4))+
  theme_minimal()

p1

##### check against patterns in raw data for 1990
#### this is very quick and dirty 
atf1990 <- atf %>% filter(year == 1990)

### note thatI'm making the colour scales the same on all these plots for easier comparison
p2 <- ggplot()+
  geom_point(data = atf1990, aes(lon, lat, colour = log(biom_kgkm2)), size = 3)+
  scale_colour_viridis_c(limits = c(-10.25, 9.4), na.value = 0.0001)+
  geom_sf(data = coast, colour = "black", fill = "grey80")+
  theme_minimal()

p2 

##### you can see that the spatial patterns predicted on the Atlantis grid recreate patterns in the trawl data pretty well  - it's getting the absence off shelf in the deep bits (or at least negligible biomass predictions), and it's getting that biomass is generally higher in the central gulf around Kodiak and somewhat lower in the far west and southeast. 

##### for comparison, let's look at the model with depth
predictions_atf_depth <- predict(m_atf_depth, newdata = atlantis_grid_df, return_tmb_object = TRUE)
atlantis_grid_df$atf_depth <- predictions_atf_depth$data$est
atlantis_grid_bgm$atf_depth <- subset(atlantis_grid_df$atf_depth, atlantis_grid_df$year == 1990)

#### this also looks pretty good, it's predicting slightly higher values to the west of Kodiak
p3 <- ggplot(subset(atlantis_grid_bgm, box_id < 92))+
  geom_sf(aes(fill=atf_depth))+
  geom_sf(data = coast_utm, colour = "black", fill = "grey80")+
  scale_fill_viridis_c(limits = c(-10.25, 9.4))+
  theme_minimal()

p3

##### and now look at the model with only space and time covariates
predictions_atf <- predict(m_atf, newdata = atlantis_grid_df, return_tmb_object = TRUE)
atlantis_grid_df$atf <- predictions_atf$data$est
atlantis_grid_bgm$atf <- subset(atlantis_grid_df$atf, atlantis_grid_df$year == 1990)

##### this model doesn't predict the offshore absences because there was no input data from there and it doesn't have depth as a covariate to inform it. The shelf part of the model looks fine though.
p4 <- ggplot(subset(atlantis_grid_bgm, box_id < 92))+
  geom_sf(aes(fill=atf))+
  geom_sf(data = coast_utm, colour = "black", fill = "grey80")+
  scale_fill_viridis_c(limits = c(-10.25, 9.4))+
  theme_minimal()

p4

##### take a look at them all together
grid.arrange(p2, p1, p3, p4, ncol = 2)

# use the inbuilt function to generate a relative biomass index
#### bias correct argument turned off for speed
ind <- get_index(predictions_atf, bias_correct = FALSE)
ind_depth <- get_index(predictions_atf_depth, bias_correct = FALSE)
ind_clim <- get_index(predictions_atf_clim, bias_correct = FALSE)

#### scale by area and convert to metric tonnes (I'm not sure the scaling is correct here because of the irregular Atlantis grid, but the trends should be)
scale <- sum(atlantis_grid_bgm$area)/109/1000

# scale the biomass by the area and units for plotting
### black is base model, purple is depth and red is climate
ggplot()+
  geom_line(data = ind, aes(year, est*scale))+
  geom_ribbon(data = ind, aes(year, est*scale, ymin = lwr*scale, ymax = upr*scale), alpha = 0.4) +
  geom_line(data = ind_depth, aes(year, est*scale), colour = "purple")+
  geom_ribbon(data = ind_depth, aes(year, est*scale, ymin = lwr*scale, ymax = upr*scale), colour = "purple", alpha = 0.4) +
  geom_line(data = ind_clim, aes(year, est*scale), colour = "red")+
  geom_ribbon(data = ind_clim, aes(year, est*scale, ymin = lwr*scale, ymax = upr*scale), colour = "red", alpha = 0.4) +
  xlab('Year') + ylab('Biomass estimate (metric tonnes)')+
  theme_minimal()

