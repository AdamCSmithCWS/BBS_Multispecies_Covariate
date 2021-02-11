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

rstan_options(auto_write = TRUE,javascript = FALSE)

library(shinystan)

# data --------------------------------------------------------------------

# trend estimates
all = read.csv("data/combined_2004_2019_Canadian_trends_and_intercepts.csv")
#weird header missmatch fix
names(all) <- c(names(all)[-1],"geometry2")


#drop the non-trend relevant columns
all <- all %>% select(route,strat,b,sd,species) 

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


buf_sel = buffer_sizes[5] #selecting the 1 km buffer

# fp_components
# [1] "cumulative"                   "built"                       
# [3] "crop"                         "dam_and_associated_reservoir"
# [5] "forestry_harvest"             "mines"                       
# [7] "nav_water"                    "night_lights"                
# [9] "oil_gas"                      "pasture"                     
# [11] "population_density"           "rail"                        
# [13] "roads" 

preds <- fp_components[c(3,8,10)]

cls_sel <- paste(preds,buf_sel,"mean",sep = "_")
cls_sel_i <- c("rt.uni",cls_sel)

fp_can_sel <- fp_can_by_route %>% select(all_of(cls_sel_i)) %>% 
  filter(rt.uni %in% rts_inc)



rts_w_preds <- unique(fp_can_sel$rt.uni)


# filter out routes with no predictor data --------------------------------

dat <- dat %>% filter(route %in% rts_w_preds) %>% 
  mutate(routeF = as.integer(factor(route)),
         speciesF = as.integer(factor(species)))

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

#matrix input (transposed to fit)
covr <- t(as.matrix(fp_can_sel[,paste(cls_sel,"standard",sep = "_")]))
ncovs <- nrow(covr)

n <- nrow(dat)


# non-breeding range covariates -------------------------------------------

sp_covs <- fp_ch_non_breeding %>% filter(english %in% species_df$species) %>% 
  select(english,
         dif_mean) %>% 
  inner_join(.,species_df,by = c("english" = "species")) %>% 
  arrange(speciesF) %>% 
  mutate(pred_scale = scale(dif_mean))

cov_nb <- as.numeric(unlist(sp_covs$pred_scale))


stan_data <- list(n = n,
                  nspecies = nspecies,
                  nroutes = nroutes,
                  ncovs = ncovs,
                  
                  slopes = dat$b,
                  sd_slopes = dat$sd,
                  species = dat$speciesF,
                  route = dat$routeF,
                  
                  cov_nb = cov_nb,
                  covr = covr)

parms = c("B",
          "b_nb",
          "b",
          "sdrte",
          "sd_b",
          "sps")


mod.file = "model/covariates_on_slopes.stan"

## compile model
model = stan_model(file=mod.file)

## run sampler on model, data
stanfit <- sampling(model,
                    data=stan_data,
                    verbose=TRUE, refresh=50,
                    chains=4, iter=2200,
                    warmup=1500,
                    cores = 4,
                    pars = parms,
                    control = list(adapt_delta = 0.95,
                                   max_treedepth = 15))






save(list = c("stanfit","stan_data","fp_can_sel","routes_df","species_df"),
     file = paste0("output/",g_sel,"_covariate_slope_output.RData"))



launch_shinystan(stanfit) 





