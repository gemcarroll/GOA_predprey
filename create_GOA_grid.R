require(sf)
require(raster)
require(rgeos)
require(ggplot2)
require(sdmTMB)

#### read in GOA RACE survye data to get dimensions
load("~/Google Drive/NOAA/AFSC/GOA food web modelling/GOA_RACE_srvy/atf.cpue.goa.Rdata")
atf <- atf.cpue.goa$data

#### read in shapefile for GOA
coast <- st_read('~/Google Drive/NOAA/AFSC/GOA food web modelling/coastfiles/GSHHS_shp/h/GSHHS_h_L1.shp')

### trim coast to relevant spatial extent
min.lo <- min(atf$LON)
max.lo <- max(atf$LON)
min.la <- min(atf$LAT)
max.la <- max(atf$LAT)

coast.sp <- as(coast, "Spatial")

coast <- crop(coast.sp, extent(min.lo+1, max.lo+1, min.la+1, max.la+1))
plot(coast)

##### read in kml file I created in Google Earth to cut off the coastal complexity so we don't make predictions in fjords & inlets etc
goa <- st_read("~/Google Drive/NOAA/AFSC/GOA food web modelling/Enviro_data/GOA.kml")
goa_poly <- st_collection_extract(
       goa,
       type = c("POLYGON"),
       warn = FALSE)
goa_poly2 <- goa_poly$geometry
goa_poly2 <- st_zm(goa_poly2)
goa.sp <- as(goa_poly2, "Spatial")
goa.sp <- crop(goa.sp, extent(min.lo+1, max.lo+1, min.la+1, max.la+1))

### create a regular grid for spatial analyses
grid <- spsample(goa.sp, cellsize = c(0.5, 0.5), type = "regular")
gridded(grid) = TRUE
plot(grid)
plot(coast, add = TRUE)

goa.grid.df <- data.frame(grid)

###### Get depth at the grid scale to subset out

bathy <- raster(x = "~/Google Drive/NOAA/AFSC/GOA food web modelling/Enviro_data/ETOPO1_Bed_g_gdal.grd")

bathy <- crop(bathy, extent(min.lo+1, max.lo+1, min.la+1, max.la+1))
plot(bathy)
plot(coast, add = T)

rasValue=extract(bathy, grid)
goa.grid.df$depth <- rasValue
names(goa.grid.df) <- c("lon", "lat", "depth")

#### remove cells that are not well sampled by the bottom trawl survey
goa.grid.df <- goa.grid.df[goa.grid.df$depth > -1000,]
goa.grid.df <- goa.grid.df[!goa.grid.df$depth > -1,]
#goa.grid.df <- goa.grid.df[!is.na(goa.grid.df$depth),]

ggplot()+
  geom_tile(data = goa.grid.df, aes(lon, lat, fill = depth))+
  geom_polygon(data = coast, aes(long, lat, group = group), colour = "black", fill = "grey80")

###### add surface and bottom temperature

year <- rep(unique(atf$YEAR), each = nrow(goa.grid.df))
lon <- rep(goa.grid.df$lon, times = length(unique(atf$YEAR)))
lat <- rep(goa.grid.df$lat, times = length(unique(atf$YEAR)))
depth <- rep(goa.grid.df$lat, times = length(unique(atf$YEAR)))

goa.grid.df <- as.data.frame(cbind(year, lon, lat, depth))
goa.grid.df$year<- as.integer(goa.grid.df$year)  

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
predictions_btemp <- predict(m_btemp, newdata = goa.grid.df, return_tmb_object = TRUE)

plot_map <- function(dat, column) {
  ggplot(dat, aes_string("lon", "lat", fill = column)) +
    geom_raster() +
    facet_wrap(~year) +
    coord_fixed()
}

plot_map(predictions_btemp$data, "est") +
  #scale_fill_viridis() +
  ggtitle("Prediction (fixed effects + all random effects)")

#### plot some individual years to get a better idea of spatial variability
plot_map(subset(predictions_btemp$data, year == 2015), "est") +
  #scale_fill_viridis() +
  ggtitle("Prediction (fixed effects + all random effects)")+
  geom_polygon(data = coast, aes(long, lat, group = group), colour = "black", fill = "grey80")

#### same thing for surface temperature
m_stemp <- sdmTMB(
  data = atf, 
  formula = SST ~ 0 + as.factor(year),
  time = "year", spde = atf_spde, family = gaussian(link = "identity"))

#### make predictions onto the GOA grid
predictions_stemp <- predict(m_stemp, newdata = goa.grid.df, return_tmb_object = TRUE)

plot_map(predictions_stemp$data, "est") +
  #scale_fill_viridis() +
  ggtitle("Prediction (fixed effects + all random effects)")+
  geom_polygon(data = coast, aes(long, lat, group = group), colour = "black", fill = "grey80")

#### add predictions to the grid
goa.grid.df$btemp <- predictions_btemp$data$est
goa.grid.df$stemp <- predictions_stemp$data$est

#### save prediction grid
write.csv(goa.grid.df, file = '~/Google Drive/NOAA/AFSC/GOA food web modelling/goa_prediction_grid.csv', row.names = F)
