require(ggplot2)
require(dplyr)
require(mapdata)
require(tidyr)

#### code to look for size-based ontogenetic breakpoints in predation by groundfish in the Gulf of Alaska

diet <- read.csv("/Users/gemmacarroll/Google Drive/NOAA/AFSC/GOA food web modelling/GOA_predprey/goa_diet_clean.csv")

#### read in file for coast
goa.map <- map('worldHires', fill=T, col='black', plot=F, xlim = c(-170, -130), ylim = c(53,60))

#### look at size variability in predators 
ggplot()+
  geom_histogram(data = diet, aes(pred_len, fill = pred))+
  facet_wrap(~pred)+
  theme_bw()

#### reduce to prey items that are most frequently consumed
mainprey <- c("Bairdi", "Motile.epifauna", "W.Pollock", "Shelf.forage", "Capelin", "Pandalidae", "Copepods", "Pred.zoop")
diet$prey_iphc <- as.character(diet$prey_iphc)
diet$prey_iphc[!diet$prey_iphc %in% mainprey] <- 'Other'

diet_mainprey <- diet %>% filter(prey_iphc %in% mainprey)

ggplot()+
  geom_histogram(data = diet_mainprey, aes(pred_len_bin, fill = pred), stat = "count")+
  facet_wrap(~pred)+
  theme_bw()

#### determine useful size bins based on ontogenetic shifts
##### prey type breakdown by predator and 10cm length bin

contrib <- diet_mainprey %>% group_by(pred, pred_len_bin, prey_iphc) %>% summarise(preycount = sum(prey_cnt, na.rm = T), preywt = sum(prey_twt)) %>% mutate(sumwt = sum(preywt), sumcnt = sum(preycount)) %>% group_by(pred, pred_len_bin, prey_iphc) %>% summarise(propwt = preywt/sumwt, propcnt = preycount/sumcnt) %>% arrange(pred, desc(propwt)) 

#### shifts in diet preference by arrowtooth flounder
#### only consider size bins with > 100 samples - this number is arbitrary, but helps to at least remove size bins that are not well-sampled

contrib_arrowtooth <- contrib %>% filter(pred == "Arrowtooth flounder") %>% as_tibble(rownames = "pred_len_bin") %>% select(-pred)

sample <- diet_mainprey %>% filter(pred == "Arrowtooth flounder") %>% group_by(pred, pred_len_bin) %>% summarise(count = length(lat)) %>% filter(count > 100)

contrib_arrowtooth <- contrib_arrowtooth %>% filter(pred_len_bin %in% sample$pred_len_bin) %>% select(-propcnt)
contrib_arrowtooth$pred_len_bin <- factor(contrib_arrowtooth$pred_len_bin, order = T)
contrib_arrowtooth <- contrib_arrowtooth[with(contrib_arrowtooth, order(pred_len_bin, prey_iphc)),]

ggplot()+
  geom_bar(data = contrib_arrowtooth, aes(pred_len_bin, propwt, fill = prey_iphc), colour = "black", size = 0.2, position = "stack", stat = "identity")+
  xlab("Size bin")+ ylab("Proportion total annual weight")+ labs(fill = "Prey group")+
  theme_classic()+
  theme(axis.title = element_text(size=30), legend.text=element_text(size=24), legend.title=element_text(size=30), axis.text = element_text(size = 18))

ggplot()+
  geom_path(data = contrib_arrowtooth, aes(pred_len_bin, propwt, colour = prey_iphc, group = prey_iphc))+
  geom_point(data = contrib_arrowtooth, aes(pred_len_bin, propwt, colour = prey_iphc, group = prey_iphc), size = 2)+
  xlab("Size bin")+ ylab("Proportion total annual weight")+ labs(fill = "Prey group")+
  theme_classic()+
  theme(axis.title = element_text(size=30), legend.text=element_text(size=24), legend.title=element_text(size=30), axis.text = element_text(size = 18))

#### use kmeans clustering as a supporting metric
contrib_arrowtooth_wide  <- contrib_arrowtooth %>%
  pivot_wider(names_from = prey_iphc, values_from = propwt) %>% select(-pred_len_bin) 

contrib_arrowtooth_wide[is.na(contrib_arrowtooth_wide)] <- 0

##### apply kmeans clustering to length bins 
kmeans(contrib_arrowtooth_wide, centers = 3)

#### arrowtooth size bins: 
#### small = 10-20cm, 20-30cm; medium = 30-40cm, 40-50cm, 50-60cm; large = 60-70cm, 70-80cm 

###################
#### shifts in diet preference by pacific cod
contrib_cod <- contrib %>% filter(pred == "Pacific cod") %>% as_tibble(rownames = "pred_len_bin") %>% select(-pred)

sample <- diet_mainprey %>% filter(pred == "Pacific cod") %>% group_by(pred, pred_len_bin) %>% summarise(count = length(lat)) %>% filter(count > 100)

contrib_cod <- contrib_cod %>% filter(pred_len_bin %in% sample$pred_len_bin) %>% select(-propcnt)
contrib_cod$pred_len_bin <- factor(contrib_cod$pred_len_bin, order = T)
contrib_cod <- contrib_cod[with(contrib_cod, order(pred_len_bin, prey_iphc)),]

ggplot()+
  geom_bar(data = contrib_cod, aes(pred_len_bin, propwt, fill = prey_iphc), colour = "black", size = 0.2, position = "stack", stat = "identity")+
  xlab("Size bin")+ ylab("Proportion total annual weight")+ labs(fill = "Prey group")+
  theme_classic()+
  theme(axis.title = element_text(size=30), legend.text=element_text(size=24), legend.title=element_text(size=30), axis.text = element_text(size = 18))

contrib_cod_wide  <- contrib_cod %>%
  pivot_wider(names_from = prey_iphc, values_from = propwt) %>% select(-pred_len_bin) 

contrib_cod_wide[is.na(contrib_cod_wide)] <- 0

##### apply kmeans clustering to length bins 
kmeans(contrib_cod_wide, centers = 3)

##### Pacific cod size bins:
#### small = 10-20cm; medium = 20-30cm, 30-40cm, large = 40-50cm, 50-60cm, 60-70cm; x.large = 70-80cm, 80-90cm

#### shifts in diet preference by walleye pollock
contrib_pollock <- contrib %>% filter(pred == "Walleye pollock") %>% as_tibble(rownames = "pred_len_bin") %>% select(-pred)

sample <- diet_mainprey %>% filter(pred == "Walleye pollock") %>% group_by(pred, pred_len_bin) %>% summarise(count = length(lat)) %>% filter(count > 100)

contrib_pollock <- contrib_pollock %>% filter(pred_len_bin %in% sample$pred_len_bin) %>% select(-propcnt)
contrib_pollock$pred_len_bin <- factor(contrib_pollock$pred_len_bin, order = T)
contrib_pollock <- contrib_pollock[with(contrib_pollock, order(pred_len_bin, prey_iphc)),]

ggplot()+
  geom_bar(data = contrib_pollock, aes(pred_len_bin, propwt, fill = prey_iphc), colour = "black", size = 0.2, position = "stack", stat = "identity")+
  xlab("Size bin")+ ylab("Proportion total annual weight")+ labs(fill = "Prey group")+
  theme_classic()+
  theme(axis.title = element_text(size=30), legend.text=element_text(size=24), legend.title=element_text(size=30), axis.text = element_text(size = 18))

ggsave(file = 'Diet data/Figures/preyprop_size_pollock.jpg', width = 15, height = 10)

contrib_pollock_wide  <- contrib_pollock %>%
  pivot_wider(names_from = prey_iphc, values_from = propwt) %>% select(-pred_len_bin) 

contrib_pollock_wide[is.na(contrib_pollock_wide)] <- 0

##### apply kmeans clustering to length bins 
kmeans(contrib_pollock_wide, centers = 3)

##### size bins for pollock
#### small = 10-20cm; medium = 20-30cm, 30-40cm, 40-50cm, 50-60cm, 60-70cm; large = 70-80cm

#### shifts in diet preference by halibut
contrib_halibut <- contrib %>% filter(pred == "Pacific halibut") %>% as_tibble(rownames = "pred_len_bin") %>% select(-pred)

sample <- diet_mainprey %>% filter(pred == "Pacific halibut") %>% group_by(pred, pred_len_bin) %>% summarise(count = length(lat)) %>% filter(count > 100)

contrib_halibut <- contrib_halibut %>% filter(pred_len_bin %in% sample$pred_len_bin) %>% select(-propcnt)
contrib_halibut$pred_len_bin <- factor(contrib_halibut$pred_len_bin, order = T)
contrib_halibut <- contrib_halibut[with(contrib_halibut, order(pred_len_bin, prey_iphc)),]

ggplot()+
  geom_bar(data = contrib_halibut, aes(pred_len_bin, propwt, fill = prey_iphc), colour = "black", size = 0.2, position = "stack", stat = "identity")+
  xlab("Size bin")+ ylab("Proportion total annual weight")+ labs(fill = "Prey group")+
  theme_classic()+
  theme(axis.title = element_text(size=30), legend.text=element_text(size=24), legend.title=element_text(size=30), axis.text = element_text(size = 18))

contrib_hal_wide  <- contrib_halibut %>%
  pivot_wider(names_from = prey_iphc, values_from = propwt) %>% select(-pred_len_bin) 

contrib_hal_wide[is.na(contrib_hal_wide)] <- 0

##### apply kmeans clustering to length bins 
kmeans(contrib_hal_wide, centers = 3)

##### size bins for halibut
#### small = 10-20cm; medium = 20-30cm, 30-40cm, 40-50cm, 50-60cm; large = 60-70cm, 70-80cm; x.large = 80-90cm, 90-100cm

#### ontogenetic shifts in diet preference by Pacific ocean perch
contrib_pop <- contrib %>% filter(pred == "Pacific ocean perch") %>% as_tibble(rownames = "pred_len_bin") %>% select(-pred)

sample <- diet_mainprey %>% filter(pred == "Pacific ocean perch") %>% group_by(pred, pred_len_bin) %>% summarise(count = length(lat)) %>% filter(count > 100)

contrib_pop <- contrib_pop %>% filter(pred_len_bin %in% sample$pred_len_bin) %>% select(-propcnt)
contrib_pop$pred_len_bin <- factor(contrib_pop$pred_len_bin, order = T)
contrib_pop <- contrib_pop[with(contrib_pop, order(pred_len_bin, prey_iphc)),]


ggplot()+
  geom_bar(data = contrib_pop, aes(pred_len_bin, propwt, fill = prey_iphc), colour = "black", size = 0.2, position = "stack", stat = "identity")+
  xlab("Size bin")+ ylab("Proportion total annual weight")+ labs(fill = "Prey group")+
  theme_classic()+
  theme(axis.title = element_text(size=30), legend.text=element_text(size=24), legend.title=element_text(size=30), axis.text = element_text(size = 18))

##### size bins for POP
#### small = 10-50cm

contrib_sab <- contrib %>% filter(pred == "Sablefish") %>% as_tibble(rownames = "pred_len_bin") %>% select(-pred)

sample <- diet_mainprey %>% filter(pred == "Sablefish") %>% group_by(pred, pred_len_bin) %>% summarise(count = length(lat)) %>% filter(count > 100)

contrib_sab <- contrib_sab %>% filter(pred_len_bin %in% sample$pred_len_bin) %>% select(-propcnt)
contrib_sab$pred_len_bin <- factor(contrib_sab$pred_len_bin, order = T)
contrib_sab <- contrib_sab[with(contrib_sab, order(pred_len_bin, prey_iphc)),]


ggplot()+
  geom_bar(data = contrib_sab, aes(pred_len_bin, propwt, fill = prey_iphc), colour = "black", size = 0.2, position = "stack", stat = "identity")+
  xlab("Size bin")+ ylab("Proportion total annual weight")+ labs(fill = "Prey group")+
  theme_classic()+
  theme(axis.title = element_text(size=30), legend.text=element_text(size=24), legend.title=element_text(size=30), axis.text = element_text(size = 18))

contrib_sab_wide  <- contrib_sab %>%
  pivot_wider(names_from = prey_iphc, values_from = propwt) %>% select(-pred_len_bin) 

contrib_sab_wide[is.na(contrib_sab_wide)] <- 0

##### apply kmeans clustering to length bins 
kmeans(contrib_sab_wide, centers = 4)

##### size bins for sablefish
#### small = 10-20cm; medium = 20-30cm, 30-40cm; large = 40-50cm, 50-60cm; x.large = 60-70cm, 70-80cm 
