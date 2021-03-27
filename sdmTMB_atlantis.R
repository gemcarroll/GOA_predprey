library(ggplot2)
library(dplyr)
library(sdmTMB)
require(sf)

##### explore functionality of sdmTMB for arrowtooth flounder biomass estimates in the GOA from RACE bottom trawl survey
load("~/Google Drive/NOAA/AFSC/GOA food web modelling/GOA_RACE_srvy/atf.cpue.goa.Rdata")
atf <- atf.cpue.goa$data

names(atf) <- c("vessel", "cruise", "haul", "hauljoin", "year", "lat", "lon", "month", "day", "stratum", "station", "btemp", "stemp", "depth", "depthr", "tempr", "bin", "biom_kgkm2", "num_km2", "species_code", "sn", "cn", "region", "LW_a", "LW_b", "srvy_nm")

#### restrict to summer survey months
atf <- atf %>% filter(month == c(6,7,8))

##### take a quick look at the data 
#### arrowtooth appear to be increasing in biomass since 1984, they seem to prefer shallower water
load('~/Google Drive/NOAA/AFSC/GOA food web modelling/GOA_predprey/goa_coast.RData')

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

### this is the mesh that the sdmTMB algorithm uses to estimate spatial autocorrelation
##### the speed is highly dependent on number of knots. 100 is quite low, you'd want to make sure it was robust by checking multiple resolutions when you have a model you want to actually use (up to like 350-450 or something)
atf_spde <- make_mesh(atf, c("lon", "lat"), n_knots = 100)
plot(atf_spde)

##### check out the distribution of the biomass density response variable. Obviously very right-skewed
hist(atf$biom_kgkm2, bins = 30)

#### it actually doesn't look terrible log transformed so let's try log transforming the response and using a normal distribution
hist(log(atf$biom_kgkm2))

#### but are there 0s? yes, some but not actually very many 
length(which(atf$biom_kgkm2 == 0))/nrow(atf)

#### because there are relatively few zeros (0.8%), transforming the response variable  will be much more efficient/appropriate than trying to fit e.g. Tweedie distribution or delta approach (delta = modeling presence/absence as binomial, then modeling biomass as , and multiplying output together). For different data sets, I would definitely reccommend exploring the second approach/es, but this will do for working example.  
### (using log1p to deal with the zeros)
hist(log1p(atf$biom_kgkm2))

##### run a model with just space and time
start.time <- Sys.time()
m_atf <- sdmTMB(
  data = atf, 
  formula = log1p(biom_kgkm2) ~ as.factor(year),
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
  family = gaussian(link = "identity"))
end.time <- Sys.time()
time.taken_m_atf <- end.time - start.time
time.taken_m_atf

#### check out model residuals
atf$resids <- residuals(m_atf) # randomized quantile residuals
hist(atf$resids)

### not amazing, but could be worse! #modelmantra
qqnorm(atf$resids)
abline(a = 0, b = 1)

##### try a slightly more complex model with smooth term for depth  
#### the s() terms are smooths, k is the smoothing parameter that controls the number of "knots" in the splines and therefore the wiggliness of the response curves
start.time <- Sys.time()
m_atf_depth <- sdmTMB(
  data = atf, 
  formula = log1p(biom_kgkm2) ~ s(depth_scale, k = 10) + as.factor(year),
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
  family = gaussian(link = "identity"))
end.time <- Sys.time()
time.taken_m_atf_depth <- end.time - start.time
time.taken_m_atf_depth

#### check out model residuals
atf$resids <- residuals(m_atf_depth) # randomized quantile residuals
hist(atf$resids)

#### you can plot the response curve from the depth smooth term - arrowtooth seem to not like deep waters
plot(m_atf_depth$mgcv_mod)

##### try an even more complex model with climate covariates (bottom temperature and surface temperature) to compare
start.time <- Sys.time()
m_atf_clim <- sdmTMB(
  data = atf, 
  formula = log1p(biom_kgkm2) ~ s(depth_scale, k = 5) + s(stemp_scale, k = 5) + s(btemp_scale, k = 5) + as.factor(year),
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
  family = gaussian(link = "identity"))
end.time <- Sys.time()
time.taken_m_atf_clim <- end.time - start.time
time.taken_m_atf_clim    

resids <- residuals(m_atf_clim$mgcv_mod) # randomized quantile residuals
hist(resids)

#### plot the response curves 
plot(m_atf_clim$mgcv_mod)

##### which of these models performed better according to AIC?
#### looks like the one with the enviro covariates did substantially better, followed by the one with depth and then the one with just space time. Let's see what they looks like spatially.
AIC(m_atf, m_atf_depth, m_atf_clim)

#### read in Atlantis prediction grids modified in Atlantis_grid_covars.R
load("/Users/gemmacarroll/Google Drive/NOAA/AFSC/GOA food web modelling/GOA_predprey/atlantis_grid_bgm.RData")
atlantis_grid_df <- read.csv('/Users/gemmacarroll/Google Drive/NOAA/AFSC/GOA food web modelling/GOA_predprey/atlantis_grid_df.csv')

atlantis_grid_df$depth_scale <- scale(atlantis_grid_df$depth)
atlantis_grid_df$btemp_scale <- scale(atlantis_grid_df$btemp)
atlantis_grid_df$stemp_scale <- scale(atlantis_grid_df$stemp)

#### make SDM predictions onto new data from climate atf model
#### I am using a dataframe version to predict onto because sdmTMB doesn't seem to support sf polygon objects
#### all years
predictions_atf_clim <- predict(m_atf_clim, newdata = atlantis_grid_df, return_tmb_object = TRUE)
atlantis_grid_df$atf_clim <- predictions_atf_clim$data$est
predictions_atf_clim1990 <- predict(m_atf_clim, newdata = subset(atlantis_grid_df, year == 1990), return_tmb_object = TRUE)
atlantis_grid_bgm$atf_clim <- predictions_atf_clim1990$data$est

###### I chose 1990 as a year to make predictions for?
#### plot predictions for 1990
#### not plotting most of Canada because the predictions look terrible (due to the temperatures being whack and not having biomass data from there)
coast <- st_as_sf(coast)
coast_utm <- st_transform(coast, crs = crs(atlantis_grid_bgm))

ggplot(subset(atlantis_grid_bgm, box_id < 92))+
  geom_sf(aes(fill=atf_clim))+
  geom_sf(data = coast_utm, colour = "black", fill = "grey80")+
  scale_fill_viridis_c()+
  theme_minimal()

##### check against patterns in raw data for 1990
#### this is very quick and dirty 
atf1990 <- atf %>% filter(year == 1990)

ggplot()+
  geom_point(data = atf1990, aes(lon, lat, colour = log1p(biom_kgkm2)), size = 3)+
  scale_colour_viridis_c()+
  geom_sf(data = coast, colour = "black", fill = "grey80")+
  theme_minimal()

#### check also depth distributions
ggplot(data = atf1990, aes(depth))+
  geom_histogram(colour = "black", fill = 'grey80')+
  theme_minimal()
  
##### you can see that the spatial patterns predicted on the Atlantis grid are alright - it's getting the absence off shelf in the deep bits (or at least negligible biomass predictions), and it's getting that biomass is generally higher in the central gulf around Kodiak and somewhat lower in the far west and southeast. The value of most of the cells (3-5ish) is similar to what you'd expect from the raw data

ggplot(data = atf1990, aes(log1p(biom_kgkm2)))+
  geom_histogram(colour = "black", fill = 'grey80')+
  theme_minimal()

ggplot(data = atlantis_grid_bgm, aes(atf_clim))+
  geom_histogram(colour = "black", fill = 'grey80')+
  theme_minimal()

##### for comparison, let's look at the model with depth
#### this also looks pretty good
predictions_atf_depth1990 <- predict(m_atf_depth, newdata = subset(atlantis_grid_df, year == 1990), return_tmb_object = TRUE)
atlantis_grid_bgm$atf_depth <- predictions_atf_depth1990$data$est

ggplot(subset(atlantis_grid_bgm, box_id < 92))+
  geom_sf(aes(fill=atf_depth))+
  geom_sf(data = coast_utm, colour = "black", fill = "grey80")+
  scale_fill_viridis_c()+
  theme_minimal()

##### and now look at the model with only space and time covariates
predictions_atf1990 <- predict(m_atf, newdata = subset(atlantis_grid_df, year == 1990), return_tmb_object = TRUE)
atlantis_grid_bgm$atf <- predictions_atf1990$data$est

ggplot(subset(atlantis_grid_bgm, box_id < 92))+
  geom_sf(aes(fill=atf))+
  geom_sf(data = coast_utm, colour = "black", fill = "grey80")+
  scale_fill_viridis_c()+
  theme_minimal()

###### seems like the climate model is most accurately reconstructing patterns in the raw data

