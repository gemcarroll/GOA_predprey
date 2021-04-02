require(sf)
require(ggplot2)
require(sdmTMB)
require(rbgm)
require(dplyr)

#### code to add enviro covars to the Atlantis prediction grid for SDM testing (these will come from ROMS in practice, I am just interpolating them to the Atlantis grid scale from the bottom trawl survey measurements using sdmTMB)

#### read in bottom trawl data (arrowtooth)
load("atlantis_test/atf.cpue.goa.Rdata")
atf <- atf.cpue.goa$data

#### restrict to summer survey months
atf <- atf %>% filter(MONTH == c(6,7,8))

#### read in Atlantis prediction grid from Albi
atlantis_grid_bgm <- read_bgm("atlantis_test/GOA_BGM_NAD83.bgm")
atlantis_grid_bgm <- atlantis_grid_bgm %>% box_sf()

#### plot depth and rename
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
area <- rep(atlantis_grid_bgm$area, times = length(unique(atf$YEAR)))

atlantis_grid_df <- as.data.frame(cbind(year, x, y, depth, area))
atlantis_grid_df$year<- as.integer(atlantis_grid_df$year)  

#### give lat lon to the grid for making predictions because that's the way I have it set up (ultimately should probably all be in UTM though)
coords <- st_as_sf(atlantis_grid_df, coords = c("x", "y"), crs = st_crs(atlantis_grid_bgm))
coords_latlon <- st_transform(coords, crs = 4326)

latlon <- data.frame(do.call(rbind, st_geometry(coords_latlon)))
names(latlon) <- c("lon", "lat")
atlantis_grid_df <- cbind(latlon, atlantis_grid_df)

#### now have a prediction grid with all grid locations and years of data
head(atlantis_grid_df)

##### interpolate temperature values from bottom trawl to the atlantis grid using sdmTMB
names(atf)[names(atf) == 'YEAR'] <- 'year'
names(atf)[names(atf) == 'LON'] <- 'lon'
names(atf)[names(atf) == 'LAT'] <- 'lat'

atf_spde <- make_mesh(atf, c("lon", "lat"), n_knots = 100)
plot(atf_spde)

##### run spatio-temporal model of bottom temperature from trawl
m_btemp <- sdmTMB(
  data = atf, 
  formula = TEMP ~ as.factor(year),
  time = "year", spde = atf_spde, family = gaussian(link = "identity"))

#### make predictions onto the grid
### couldn't get these predictions to work directly on the bgm grid, so using a dataframe intermediate step
#### predicting for all years onto the full grid
predictions_btemp <- predict(m_btemp, newdata = atlantis_grid_df, return_tmb_object = TRUE)
atlantis_grid_df$btemp <- predictions_btemp$data$est

#### Adding predictions for 1990 to Atlantis grid
atlantis_grid_bgm$btemp <- subset(atlantis_grid_df$btemp, atlantis_grid_df$year == 1990)

#### here's interpolated bottom temp predicted for 1990 from the bottom trawls
#### you can clearly see the west to east seasonal progression in sampling...! anyway, this will do as an example
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

#### predicting surface temp onto grid
predictions_stemp <- predict(m_stemp, newdata = atlantis_grid_df, return_tmb_object = TRUE)
atlantis_grid_df$stemp <- predictions_stemp$data$est

#### Adding predictions for 1990 to Atlantis grid
atlantis_grid_bgm$stemp <- subset(atlantis_grid_df$stemp, atlantis_grid_df$year == 1990)

ggplot(subset(atlantis_grid_bgm, box_id < 92))+
  geom_sf()+
  geom_sf(aes(fill=stemp))+
  scale_fill_viridis_c()+
  theme_minimal()

#### save these grid objects to make example SDM predictions onto 
write.csv(atlantis_grid_df, file = 'atlantis_test/atlantis_grid_df.csv', row.names = F)
save(atlantis_grid_bgm, file = "atlantis_test/atlantis_grid_bgm.RData")
