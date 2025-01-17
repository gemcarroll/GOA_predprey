require(dplyr)
require(stringr)

##### code to combine biomass cpue data for 6 predators and add length bins prior to building species distribution models. 

##### create size-structured data set for all predators to facilitate efficient model building
load("~/Google Drive/NOAA/AFSC/GOA food web modelling/GOA_RACE_srvy/atf.cpue.goa.Rdata")
atf <- atf.cpue.goa$data
load("~/Google Drive/NOAA/AFSC/GOA food web modelling/GOA_RACE_srvy/sablefish.cpue.goa.Rdata")
sablefish <- sablefish.cpue.goa$data
load("~/Google Drive/NOAA/AFSC/GOA food web modelling/GOA_RACE_srvy/plk.cpue.goa.Rdata")
pollock <- plk.cpue.goa$data
load("~/Google Drive/NOAA/AFSC/GOA food web modelling/GOA_RACE_srvy/halibut.cpue.goa.Rdata")
halibut <- halibut.cpue.goa$data
load("~/Google Drive/NOAA/AFSC/GOA food web modelling/GOA_RACE_srvy/pcod.cpue.goa.Rdata")
pcod <- pcod.cpue.goa$data
pop <- read.csv("/Users/gemmacarroll/Google Drive/NOAA/AFSC/GOA food web modelling/GOA_RACE_srvy/POP_data_survey47/data.csv")
pop <- pop[, 2:length(pop)]

#### bind datasets together
preds <- rbind(atf, sablefish, pollock, halibut, pcod, pop)

names(preds) <- c("vessel", "cruise", "haul", "hauljoin", "year", "lat", "lon", "month", "day", "stratum", "station", "btemp", "stemp", "depth", "depthr", "tempr", "bin", "biom_kgkm2", "num_km2", "species_code", "species_name", "pred", "region", "LW_a", "LW_b", "srvy_nm")

##### add zone information 
strata <- read.csv('/Users/gemmacarroll/Google Drive/NOAA/AFSC/GOA food web modelling/Diet data/table_strata.csv')

### reduce to just column of interest for joining (IPHC)
strata <- strata %>% filter(Region == "GOA") %>% select(Stratum, ECOPATH_PP, area, mindepth, maxdepth) %>% rename(stratum = Stratum, zone = ECOPATH_PP)
preds <- left_join(preds, strata, by = "stratum")
preds$zone2 <- str_sub(preds$zone, 1, -6)
preds$zone2[is.na(preds$zone2)] <- "Southeast"
preds$zone2[preds$zone2 == "Southeast" & preds$lon < -140] <- NA

ggplot(data = preds, aes(lon, lat, colour = zone2))+
  geom_point()

##### add size class identifiers based on explorations conducted in 'size_classes.R'
preds$length_class <- NA

#### arrowtooth 
preds$length_class[preds$pred == "arrowtooth flounder" & preds$bin > 99 & preds$bin < 300] <- "small"
preds$length_class[preds$pred == "arrowtooth flounder" & preds$bin > 299 & preds$bin < 600] <- "medium"  
preds$length_class[preds$pred == "arrowtooth flounder" & preds$bin > 599 & preds$bin < 800] <- "large"    
  
#### sablefish
#### small = 10-20cm; medium = 20-30cm, 30-40cm; large = 40-50cm, 50-60cm; x.large = 60-70cm, 70-80cm 
preds$length_class[preds$pred == "sablefish" & preds$bin > 99 & preds$bin < 200] <- "small"
preds$length_class[preds$pred == "sablefish" & preds$bin > 199 & preds$bin < 400] <- "medium"  
preds$length_class[preds$pred == "sablefish" & preds$bin > 399 & preds$bin < 600] <- "large"   
preds$length_class[preds$pred == "sablefish" & preds$bin > 599 & preds$bin < 800] <- "x.large"     

#### walleye pollock
#### small = 10-20cm; medium = 20-30cm, 30-40cm, 40-50cm, 50-60cm, 60-70cm; large = 70-80cm
preds$length_class[preds$pred == "walleye pollock" & preds$bin > 99 & preds$bin < 200] <- "small"
preds$length_class[preds$pred == "walleye pollock" & preds$bin > 199 & preds$bin < 700] <- "medium"  
preds$length_class[preds$pred == "walleye pollock" & preds$bin > 699 & preds$bin < 800] <- "large"   

#### halibut
#### small = 10-20cm; medium = 20-30cm, 30-40cm, 40-50cm, 50-60cm; large = 60-70cm, 70-80cm; x.large = 80-90cm, 90-100cm
preds$length_class[preds$pred == "Pacific halibut" & preds$bin > 99 & preds$bin < 200] <- "small"
preds$length_class[preds$pred == "Pacific halibut" & preds$bin > 199 & preds$bin < 600] <- "medium"  
preds$length_class[preds$pred == "Pacific halibut" & preds$bin > 599 & preds$bin < 800] <- "large"   
preds$length_class[preds$pred == "Pacific halibut" & preds$bin > 799 & preds$bin < 1000] <- "x.large" 

#### Pacific cod
#### small = 10-20cm; medium = 20-40cm, large = 40-70cm; x.large = 70-80cm, 80-90cm
preds$length_class[preds$pred == "Pacific cod" & preds$bin > 99 & preds$bin < 200] <- "small"
preds$length_class[preds$pred == "Pacific cod" & preds$bin > 199 & preds$bin < 400] <- "medium"  
preds$length_class[preds$pred == "Pacific cod" & preds$bin > 399 & preds$bin < 700] <- "large"   
preds$length_class[preds$pred == "Pacific cod" & preds$bin > 699 & preds$bin < 900] <- "x.large"

#### Pacific ocean perch
#### small = 10-20cm; medium = 20-40cm; large = 40-50cm
preds$length_class[preds$pred == "Pacific ocean perch" & preds$bin > 99 & preds$bin < 200] <- "small"
preds$length_class[preds$pred == "Pacific ocean perch" & preds$bin > 199 & preds$bin < 400] <- "medium"  
preds$length_class[preds$pred == "Pacific ocean perch" & preds$bin > 399 & preds$bin < 500] <- "large"  

#### remove length bins that are not represented well by the diet data
preds <- preds %>% filter(! is.na(preds$length_class))

#### use length bins based on ontogenetic shifts in diet to identify classes for modelling
preds$species_length_class <- NA
preds$species_length_class[preds$pred == 'arrowtooth flounder'] <- "atf"
preds$species_length_class[preds$pred == 'Pacific cod'] <- "cod"
preds$species_length_class[preds$pred == 'Pacific ocean perch'] <- "pop"
preds$species_length_class[preds$pred == 'Pacific halibut'] <- "hal"
preds$species_length_class[preds$pred == 'sablefish'] <- "sab"
preds$species_length_class[preds$pred == 'walleye pollock'] <- "pol"
preds$species_length_class <- paste(preds$species_length_class, preds$length_class, sep = "_")

#### remove small sablefish class because there are only 23 observations
preds <- preds %>% filter(! species_length_class == "sab_small")

##### reduce to just summer survey
preds <- preds %>% filter(month == c(6,7,8))

##### get total biomass cpue for each size class 
preds_sum <- preds %>% group_by(year, station, pred, species_length_class) %>% summarise(biom_kgkm2 = sum(biom_kgkm2))

###### add 0s into the data
station <- preds %>% distinct(year, station)
species <- preds %>% distinct(year, pred, species_length_class)
lonlat <- preds %>% distinct(year, station, lon, lat, month, depth, btemp, stemp)

pred_grid <- left_join(station, species)
preds_sum2 <- left_join(pred_grid, lonlat)
preds_sum3 <- left_join(preds_sum2, preds_sum)

preds_sum3$biom_kgkm2[is.na(preds_sum3$biom_kgkm2)] <- 0

###### order by year and station
pred <- preds_sum3[with(preds_sum3, order(year, station, species_length_class)),]

##### save file 
write.csv(pred, file = "/Users/gemmacarroll/Google Drive/NOAA/AFSC/GOA food web modelling/GOA_predprey/goa_biomass_clean.csv", row.names = F)

