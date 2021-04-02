require(dplyr)
require(stringr)
require(ggplot2)
require(mapdata)

#### code to tidy up Groundfish diet data in Gulf of Alaska

### read in data files (GOA diet data, NODC codes & higher level classifications)
diet <- read.csv('/Users/gemmacarroll/Google Drive/NOAA/AFSC/GOA food web modelling/Diet data/Data_March2021/PredPreyForWeb.csv')
table_nodc <- read.csv('/Users/gemmacarroll/Google Drive/NOAA/AFSC/GOA food web modelling/Diet data/table_nodc.csv')
lookup <- read.csv('/Users/gemmacarroll/Google Drive/NOAA/AFSC/GOA food web modelling/Diet data/GOA_med_diet_lookup.csv')
strata <- read.csv('/Users/gemmacarroll/Google Drive/NOAA/AFSC/GOA food web modelling/Diet data/table_strata.csv')

#### read in file for coast
goa.map <- map('worldHires', fill=T, col='black', plot=F, xlim = c(-170, -130), ylim = c(53,60))

#### reduce diet data set to key variables and rename
diet <- diet %>% select(-VESSEL, -HAUL, -CRUISE, -PREY_PARTS, -PRED_DIG, -REGION, -TYPE) %>% rename(pred = Pred_common, pred_species = Pred_Species, predhaul_join = HAULJOIN, pred_nodc = PRED_NODC, pred_specimen = PRED_SPECN, prey_nodc = PREY_NODC, prey_cnt = PREY_CNT, prey_twt = PREY_TWT, predjoin = PREDJOIN, pred_stomwt = PRED_STOMWT, pred_len = PRED_LEN, pred_full = PRED_FULL, pred_wt = PRED_WT, year = Year, month = Month, prey_group_name = Prey_group_Name, prey_name = Prey_Name, pred_lh = PRED_LH, prey_lh = PREY_LH, pred_sex = PRED_SEX, day = Day, lat = RLAT, lon = RLONG, gear_depth = GEAR_DEPTH, bottom_depth = BOTTOM_DEPTH, start_hour = START_HOUR, stemp = SURFACE_TEMP, btemp = GEAR_TEMP, inpfc_area = INPFC_AREA, station_id = STATIONID, start_date = START_DATE, stratum = STRATUM, cruise_type = CRUISE_TYPE) 

#### reduce data set to key predators
preds <- c("Walleye pollock", "Pacific cod", "Pacific halibut", "Pacific ocean perch", "Arrowtooth flounder", "Sablefish")
diet <- diet %>% filter(pred %in% preds) 

#### reduce to just column of interest for joining (ECOPATH_Prey)
table_nodc <- table_nodc[, c(1, 9)]

### reduce to just column of interest for joining (IPHC)
lookup <- lookup[, c(1, 4)]
names(lookup) <- c("PREYNAME", "prey_iphc")

### reduce to just column of interest for joining (IPHC)
strata <- strata %>% filter(Region == "GOA") %>% select(Stratum, ECOPATH_PP, area, mindepth, maxdepth) %>% rename(stratum = Stratum, zone = ECOPATH_PP)

### join data by prey nodc cod  
diet <- left_join(diet, table_nodc, by = c("prey_nodc" = "nodc"))
diet <- left_join(diet, lookup, by = c("ECOPATH_Prey" = "PREYNAME"))
diet <- left_join(diet, strata, by = "stratum")

diet$prey_iphc <- as.character(diet$prey_iphc)
diet$prey_iphc[which(is.na(diet$prey_iphc))] <- "Empty"
diet$zone2 <- str_sub(diet$zone, 1, -6)
diet$zone2[is.na(diet$zone2)] <- "Southeast"
diet$zone2[diet$zone2 == "Southeast" & diet$lon < -140] <- NA

#### reorder columns
diet <- diet %>% select(year, month, day, lat, lon, pred, pred_species, pred_len, pred_wt, pred_sex, pred_stomwt, pred_full, prey_name, prey_group_name, ECOPATH_Prey, prey_iphc, prey_cnt, prey_twt, everything())

### only include RACE groundfish survey years
survey <- diet %>% group_by(year) %>% summarise(race = any(cruise_type == "Race_Groundfish")) %>% filter(race == "TRUE")

diet <- diet %>% filter(year %in% survey$year)

#### plot sampling coverage by year
ggplot(data = diet, aes(lon, lat))+
  geom_point()+
  facet_wrap(~year)

#### plot sampling coverage by year
ggplot(data = subset(diet, pred == "Walleye pollock"), aes(lon, lat))+
  geom_point()+
  facet_wrap(~year)

#### remove years before 1990
diet <- diet %>%  filter(year > 1989)

#### Look at spatial sampling effort by month
ggplot()+
  geom_path(data = goa.map, aes(long, lat, group = group))+
  geom_point(data = diet, aes(lon, lat), col = "red")+
  xlim(c(-170,-125)) +ylim(c(51, 62))+
  facet_wrap(~month)+
  theme_classic()+
  theme(strip.text = element_text(size=16))

### subset to May-August samples (other time periods represent observer data and out-of-season data may not be comparable in terms of ecology)
diet <- diet[diet$month > 4 & diet$month < 9,]

#### add 10cm length bins
diet$pred_len_bin <- substring(diet$pred_len, 1, 1)
diet$pred_len_bin[diet$pred_len < 10] <- 0

#### plot all sample locations for each species 
ggplot(data = diet, aes(lon, lat))+
  geom_point()+
  facet_wrap(~pred)

#### plot all sample locations for each species by zone
ggplot(data = diet, aes(lon, lat))+
  geom_point(aes(colour = zone2))+
  facet_wrap(~pred)

write.csv(diet, file = "/Users/gemmacarroll/Google Drive/NOAA/AFSC/GOA food web modelling/GOA_predprey/goa_diet_clean.csv", row.names = F)

