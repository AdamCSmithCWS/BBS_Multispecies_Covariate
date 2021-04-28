//
// modify this to estimate covariate effects separately by strata - random effects

// 
data {


// int<lower=0,upper=1> I[nspecies,nroutes]; // indicator matrix for species by route inclusions

int<lower=0> n;//number of slope estimates 
int<lower=0> nspecies;//number of species 
int<lower=0> nroutes;//number of routes 
int<lower=0> npredictions;//number of values in prediction smooth
int<lower=0> nstrata;//number of strata


real slopes[n];//vector of slope estimates
real sd_slopes[n];// vector of observed sd of slope estimates
int species[n];// species indicators
int route[n];// route indicators
int strata[n];// strata indicators

real cov_nb[nspecies];// species-level covariates - change in footprint on wintering grounds

vector[nroutes] covr; // route-level predictors
 
 
real cov_pred[npredictions]; // vector covering range of predictors
 
 
}

// 
parameters {

//vector[nroutes] rte_raw;
//real<lower=0> sdrte;
real SPS[nspecies]; //species hyperparameters on trend intercepts
matrix[nspecies,nstrata] sps_raw; //species intercept on trends by strata - random
real<lower=0> sd_sps[nspecies];

matrix[nstrata,nspecies] b_raw_group;
matrix[nspecies,nstrata] b_raw_species;
vector[nstrata] B_group_raw; //group level hyper parameter in a stratum
vector[nspecies] B_species_raw; // species level hyperparameter - covariate effect
real B_GROUP; //Group level mean effect of covariate across all strata
//real B_SPECIES; // species level mean effect of covariate across all strata


real<lower=0> sd_B_group; //sd on 
real<lower=0> sd_B_species;
real<lower=0> sd_b_group[nstrata];
real<lower=0> sd_b_species[nspecies];
real B;
real b_nb;
real<lower=0> epsilon;
vector[n] noise_raw;

}

transformed parameters{

  //vector[nroutes] rte;

  real mu_slopes[n];
  matrix[nstrata,nspecies] b; //combined species covariate response
  vector[nstrata] B_group; // group hyper-hyper param
  matrix[nstrata,nspecies] b_group; // group-component of species resonse
  
  matrix[nspecies,nstrata] b_species;// species-component of species resonse
  vector[nspecies] B_species; // species hyper-hyper param
  matrix[nspecies,nstrata] sps; //species intercept on trends by strata - random
  vector[n] noise;
  real true_slopes[n];

  
  noise = epsilon*noise_raw;
 

 //group level component of the response to covariate - varies among species
 // centered on a group-strata hyperparameter, centered on a group- hyper-hyper param
B_group  = sd_B_group*B_group_raw + B_GROUP;
 for(t in 1:nstrata){
 b_group[t,] = sd_b_group[t]*b_raw_group[t,] + B_group[t];
}
   
   //species level component in the response to covariate - same as above but
   // varies among strata for a given species
 B_species = sd_B_species*B_species_raw;
for(s in 1:nspecies){
    sps[s,] = sd_sps[s]*sps_raw[s,] + SPS[s]; //species intercept on the trends

 b_species[s,] = sd_b_species[s]*b_raw_species[s,] + B_species[s];
}


for(t in 1:nstrata){
  for(s in 1:nspecies){
    b[t,s] = b_group[t,s] + b_species[s,t];
  }}


  
for(i in 1:n){
mu_slopes[i] = (sps[species[i],strata[i]] + b_nb*cov_nb[species[i]] + b[strata[i],species[i]]*covr[route[i]]);// + rte[route[i]]);

   true_slopes[i] = mu_slopes[i] + noise[i];

}



}
// 
model {


// random route effects on trends to account for repeated samples at routes
//rte_raw ~ normal(0,1);
//sum(rte_raw) ~ normal(0,0.0001*nroutes);//sum to zero soft constraint on route level random effects

//sdrte ~ gamma(2,0.1);//boundary avoiding prior on // ~ normal(0,0.01);

// Group level effects
B_group_raw ~ std_normal();
sum(B_group_raw) ~ normal(0,0.0001*nstrata);
B_GROUP ~ normal(0,0.1);
sd_B_group ~ normal(0,0.1);
sd_b_group ~ normal(0,0.1);
for(t in 1:nstrata){
    b_raw_group[t,] ~ std_normal(); // species random effects of covariates
  sum(b_raw_group[t,]) ~ normal(0,0.0001*nstrata);//sum to zero soft constraint on species level random effects of covariates

}



// species effects
SPS ~ normal(0,0.1);
sd_sps ~ normal(0,0.1);

B_species_raw ~ std_normal();
sum(B_species_raw) ~ normal(0,0.0001*nspecies);

sd_B_species ~ normal(0,0.1);
sd_b_species ~ normal(0,0.1);
for(s in 1:nspecies){
  sps_raw[s,] ~ std_normal();
    b_raw_species[s,] ~ std_normal(); // species random effects of covariates
  sum(b_raw_species[s,]) ~ normal(0,0.0001*nspecies);//sum to zero soft constraint on species level random effects of covariates

}


//espilon
epsilon ~ gamma(2,0.1);//boundary avoiding prior // ~ normal(0,0.01);// process error

// non-breeding range footprint effect
b_nb ~ normal(0,0.01);

// likelihood
// process model
   noise_raw ~ std_normal();//normal(mu_slopes,epsilon);
// observation model   
for(i in 1:n){
   slopes[i] ~ normal(true_slopes[i],sd_slopes[i]);
    } 
  
}

generated quantities {

  real preds[nstrata,nspecies,npredictions];

for(s in 1:nspecies){
  for(t in 1:nstrata){
for(j in 1:npredictions){
preds[t,s,j] = (sps[s,t] + b_nb*cov_nb[s] + b[t,s]*cov_pred[j]);// + rte[route[i]]);
}
}
}


}



