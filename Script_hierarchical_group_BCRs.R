####

# script to run a multispecies analysis of estimated trends on BBS routes
# using various route-level covariates to predict the trends

# trends are produced from a spatial iCAR model

# requires:
 ## data frame with columns for species, route, slope, sd-slope
 ## matrix of covariates nrows = ncovs, ncols = nspecies
 ## vector of species-level non-breeding region covariates (change in footprint)
 ## possible addition vector of species-level breeding region covariates (change in footprint)


library(tidyverse)
library(rstan)
library(tidybayes)

rstan_options(auto_write = TRUE,javascript = FALSE)

library(shinystan)

# data --------------------------------------------------------------------

# trend estimates
all = read.csv("data/combined_2004_2019_Canadian_trends_and_intercepts2.csv")
#weird header missmatch fix
names(all) <- c(names(all)[-1],"geometry2")


#drop the non-trend relevant columns
all <- all %>% select(route,strat,b,sd,species) %>% 
  mutate(bcr = str_extract(strat,pattern = "[:digit:]+"))


# species available
sps_all <- unique(all$species)


# load species groupings from State of Canada's Birds
socb = read.csv("data/SOCB data supplement.csv")

groups = c("Grassland.birds",
           "Forest.birds",
           "Other.birds",
           "Aerial.insectivores",
           "suburban",
           "other.wetland.birds")



# group loop --------------------------------------------------------------

g_sel = groups[1]


# species selection of trend data -----------------------------------------


socb_c = which(grepl(names(socb),pattern = g_sel))

sp_sel1 = socb[which(socb[,socb_c] == "Included in group"),"species"]

sp_sel2 = sp_sel1[which(sp_sel1 %in% sps_all)]


# covariate data ----------------------------------------------------------

load("data/fp_change.RData")

# non-breeding range change in footprint by species
sp_w_fp_ch_dat <- fp_ch_non_breeding$english
#some species don't have range data on footprint
sp_sel = sp_sel2[which(sp_sel2 %in% sp_w_fp_ch_dat)]


dat = all %>% filter(species %in% sp_sel)

rts_inc = unique(dat$route)


load("data/compiled_footprint_data.RData")


buf_sel = buffer_sizes[5] #selecting the 4 km buffer

# fp_components
# [1] "cumulative"                   "built"                       
# [3] "crop"                         "dam_and_associated_reservoir"
# [5] "forestry_harvest"             "mines"                       
# [7] "nav_water"                    "night_lights"                
# [9] "oil_gas"                      "pasture"                     
# [11] "population_density"           "rail"                        
# [13] "roads" 

preds <- fp_components[13]

cls_sel <- paste(preds,buf_sel,"mean",sep = "_")
cls_sel_i <- c("rt.uni",cls_sel)

fp_can_sel <- fp_can_by_route %>% select(all_of(cls_sel_i)) %>% 
  filter(rt.uni %in% rts_inc)



rts_w_preds <- unique(fp_can_sel$rt.uni)


# filter out routes with no predictor data --------------------------------

dat <- dat %>% filter(route %in% rts_w_preds) %>% 
  mutate(routeF = as.integer(factor(route)),
         speciesF = as.integer(factor(species)),
         stratF = as.integer(factor(bcr)))

strata_df <- unique(dat[,c("bcr","stratF")])
strata_df <- arrange(strata_df,stratF)
nstrata <- max(strata_df$stratF)


routes_df <- unique(dat[,c("route","routeF")])
routes_df <- arrange(routes_df,routeF)
nroutes <- max(routes_df$routeF)

species_df <- unique(dat[,c("species","speciesF")])
species_df <- arrange(species_df,speciesF)

nspecies <- max(species_df$speciesF)

# align predictor route indicators ----------------------------------------

fp_can_sel = inner_join(routes_df,fp_can_sel,by = c("route" = "rt.uni")) %>% 
  arrange(routeF)

for(cl in cls_sel){
  fp_can_sel[,paste(cl,"standard",sep = "_")] <- scale(fp_can_sel[,cl])
}

#matrix input 
covr <- as.numeric(unlist(fp_can_sel[,paste(cls_sel,"standard",sep = "_")]))

n <- nrow(dat)



# non-breeding range covariates -------------------------------------------

sp_covs <- fp_ch_non_breeding %>% filter(english %in% species_df$species) %>% 
  select(english,
         dif_mean) %>% 
  inner_join(.,species_df,by = c("english" = "species")) %>% 
  arrange(speciesF) %>% 
  mutate(pred_scale = scale(dif_mean))

cov_nb <- as.numeric(unlist(sp_covs$pred_scale))

npredictions <- 100
cov_pred <- seq(min(covr),max(covr),length = npredictions)



stan_data <- list(n = n,
                  nspecies = nspecies,
                  nroutes = nroutes,
                  npredictions = npredictions,
                  nstrata = nstrata,
                  
                  slopes = dat$b,
                  sd_slopes = dat$sd,
                  species = dat$speciesF,
                  route = dat$routeF,
                  strata = dat$stratF,
                  
                  cov_nb = cov_nb,
                  covr = covr,
                  cov_pred = cov_pred)

parms = c("b_nb",
          "b",
          "b_species",
          "b_group",
          #"sdrte",
          #"rte",
          "preds",
          "sd_b_species",
          "sd_b_group",
          "sps",
          "SPS",
          "epsilon")


mod.file = "model/One_covariate_on_slopes_strata2.stan"

## compile model
model = stan_model(file=mod.file)

## run sampler on model, data
stanfit <- sampling(model,
                    data=stan_data,
                    verbose=TRUE, refresh=50,
                    chains=4, iter=2000,
                    warmup=1000,
                    cores = 4,
                    pars = parms,
                    control = list(adapt_delta = 0.8,
                                   max_treedepth = 15))






save(list = c("stanfit","stan_data","fp_can_sel","routes_df","species_df","dat",
              "cls_sel","strata_df"),
     file = paste0("output/",g_sel,"_",preds,"_covariate_slope_output.RData"))



launch_shinystan(stanfit) 






load(paste0("output/",g_sel,"_covariate_slope_output.RData"))



# plot trends vs covariate ------------------------------------------------

# t_pch <- function(x){
#   (exp(x)-1)*100
# }

rescale <- sd(fp_can_sel[,cls_sel])
recenter <- mean(fp_can_sel[,cls_sel])

b_species_samp <- spread_draws(stanfit,b_species[speciesF,stratF]) %>% 
  mutate(b_species = b_species/rescale) %>% 
  left_join(.,species_df,by = "speciesF") %>%
  left_join(.,strata_df,by = "stratF") %>%
  ungroup() %>% 
  group_by(species,bcr) %>% 
  summarise(mean = mean(b_species),
            lci = quantile(b_species,0.025),
            uci = quantile(b_species,0.975))

b_species_plot <- ggplot(data = b_species_samp,aes(x = species,y = mean))+
  geom_point()+
  geom_errorbar(aes(ymin = lci,ymax = uci),alpha = 0.2, width = 0)+
  geom_abline(intercept = 0,slope = 0)+
  coord_flip()+
  facet_wrap(~bcr,nrow = 5,ncol = 3)

print(b_species_plot)  



b_group_samp <- spread_draws(stanfit,b_group[stratF,speciesF]) %>% 
  mutate(b_group = b_group/rescale)


b_samp <- spread_draws(stanfit,b[stratF,speciesF])  %>% 
  mutate(b = b/rescale) %>% 
  left_join(.,species_df,by = "speciesF") %>%
  left_join(.,strata_df,by = "stratF") %>%
  ungroup() %>% 
  group_by(species,bcr) %>% 
  summarise(mean = mean(b),
         lci = quantile(b,0.025),
         uci = quantile(b,0.975))

b_plot <- ggplot(data = b_samp,aes(x = species,y = mean))+
  geom_point()+
  geom_errorbar(aes(ymin = lci,ymax = uci),alpha = 0.2, width = 0)+
  geom_abline(intercept = 0,slope = 0)+
  coord_flip()+
  facet_wrap(~bcr,nrow = 5,ncol = 3)
  
print(b_plot)  

cov_resc <- data.frame(p = 1:length(stan_data$cov_pred),
                       cov_pred = stan_data$cov_pred,
                       cov_orig = (stan_data$cov_pred*rescale)+recenter)

b_preds <- spread_draws(stanfit,preds[speciesF,p])  %>% 
  left_join(.,species_df,by = "speciesF") %>%
  ungroup() %>% 
  group_by(species,p) %>% 
  summarise(mean = mean(preds),
            lci = quantile(preds,0.025),
            uci = quantile(preds,0.975)) %>% 
left_join(.,cov_resc,by = "p")

names(fp_can_sel)[3] <- "cov_orig"


dat <- left_join(dat,fp_can_sel,by = c("route","routeF")) 


sl_plot = ggplot(data = b_preds,aes(x = cov_orig,y = mean))+
  geom_point(data = dat,aes(x = cov_orig,y = b,size = 1/sd),alpha = 0.1)+
  geom_line(aes(colour = species))+
  geom_ribbon(aes(ymin = lci,ymax = uci,fill = species),alpha = 0.1)+
  facet_wrap(~species,nrow = ceiling(sqrt(stan_data$nspecies)))
print(sl_plot)
