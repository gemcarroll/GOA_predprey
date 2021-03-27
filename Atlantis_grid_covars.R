require(sf)
require(ggplot2)
require(sdmTMB)
require(rbgm)
require(dplyr)

#### code to add enviro covars to the Atlantis prediction grid for SDM predicitons (these will come from ROMS in practice, I am just interpolating them to the Atlantis grid scale from the bottom trawl survey measurements using sdmTMB)

#### read in bottom trawl data (arrowtooth)
load("~/Google Drive/NOAA/AFSC/GOA food web modelling/GOA_RACE_srvy/atf.cpue.goa.Rdata")
atf <- atf.cpue.goa$data

#### restrict to summer survey months
atf <- atf %>% filter(MONTH == c(6,7,8))

#### read in Atlantis prediction grid 
atlantis_grid_bgm <- read_bgm("~/Google Drive/NOAA/AFSC/GOA food web modelling/Atlantis_grid/GOA_BGM_NAD83.bgm")
atlantis_grid_bgm <- atlantis_grid_bgm %>% box_sf()

ggplot(atlantis_grid_bgm)+
  geom_sf(aes(fill=botz))+
  scale_fill_viridis_c()+
  theme_minimal()

names(atlantis_grid_bgm)[3] <- "depth"

###### create a prediction grid with all years of bottom trawl data 
year <- rep(unique(atf$YEAR), each = nrow(atlantis_grid_bgm))
x <- rep(atlantis_grid_bgm$insideX, times = length(unique(atf$YEAR)))
y <- rep(atlantis_grid_bgm$insideY, times = length(unique(atf$YEAR)))
depth <- rep(atlantis_grid_bgm$depth, times = length(unique(atf$YEAR)))

atlantis_grid_df <- as.data.frame(cbind(year, x, y, depth))
atlantis_grid_df$year<- as.integer(atlantis_grid_df$year)  

#### doing this because I'm too lazy to convert the SDM code to UTM, but that should prob be in UTM too
coords <- st_as_sf(atlantis_grid_df, coords = c("x", "y"), crs = crs(atlantis_grid_bgm))
coords_latlon <- st_transform(coords, crs = 4326)

latlon <- data.frame(do.call(rbind, st_geometry(coords_latlon)))
names(latlon) <- c("lon", "lat")
atlantis_grid_df <- cbind(latlon, atlantis_grid_df)

#### add lat lon to bgm grid
atlantis_grid_bgm$lon <- latlon[1:109, 1]
atlantis_grid_bgm$lat <- latlon[1:109, 2]

##### interpolate temperature values from bottom trawl to the atlantis grid
names(atf)[names(atf) == 'YEAR'] <- 'year'
names(atf)[names(atf) == 'LON'] <- 'lon'
names(atf)[names(atf) == 'LAT'] <- 'lat'

atf_spde <- make_mesh(atf, c("lon", "lat"), n_knots = 100)
plot(atf_spde)

##### run spatio-temporal model of bottom temperature
m_btemp <- sdmTMB(
  data = atf, 
  formula = TEMP ~ 0 + as.factor(year),
  time = "year", spde = atf_spde, family = gaussian(link = "identity"))

#### make predictions onto the GOA grid
### couldn't get these predictions to work directly on the bgm grid, so using a dataframe intermediate step
#### predicting for 1990 (all covars in model + lat, lon, and year must be in the 'newdata' object)
predictions_btemp <- predict(m_btemp, newdata = atlantis_grid_df, return_tmb_object = TRUE)
atlantis_grid_df$btemp <- predictions_btemp$data$est
predictions_btemp1990 <- predict(m_btemp, newdata = subset(atlantis_grid_df, year == 1990), return_tmb_object = TRUE)
atlantis_grid_bgm$btemp <- predictions_btemp1990$data$est

#### here's interpolated bottom temp predicted for 1990 from the bottom trawls
#### you can clearly see the west to east progression in sampling... this will do as an example
ggplot(subset(atlantis_grid_bgm, box_id < 92))+
  geom_sf()+
  geom_sf(aes(fill=btemp))+
  scale_fill_viridis_c()+
  theme_minimal()

#### do same thing for surface temperature
m_stemp <- sdmTMB(
  data = atf, 
  formula = SST ~ 0 + as.factor(year),
  time = "year", spde = atf_spde, family = gaussian(link = "identity"))

#### predicting surface temp
predictions_stemp <- predict(m_stemp, newdata = atlantis_grid_df, return_tmb_object = TRUE)
atlantis_grid_df$stemp <- predictions_stemp$data$est
predictions_stemp1990 <- predict(m_stemp, newdata = subset(atlantis_grid_df, year == 1990), return_tmb_object = TRUE)
atlantis_grid_bgm$stemp <- predictions_stemp1990$data$est

ggplot(subset(atlantis_grid_bgm, box_id < 92))+
  geom_sf()+
  geom_sf(aes(fill=stemp))+
  scale_fill_viridis_c()+
  theme_minimal()

#### save these grid objects to make example SDM predictions onto 
write.csv(atlantis_grid_df, file = '/Users/gemmacarroll/Google Drive/NOAA/AFSC/GOA food web modelling/GOA_predprey/atlantis_grid_df.csv', row.names = F)
save(atlantis_grid_bgm, file = "/Users/gemmacarroll/Google Drive/NOAA/AFSC/GOA food web modelling/GOA_predprey/atlantis_grid_bgm.RData")
