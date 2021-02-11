//
//

// 
data {


// int<lower=0,upper=1> I[nspecies,nroutes]; // indicator matrix for species by route inclusions

int<lower=0> n;//number of slope estimates 
int<lower=0> nspecies;//number of species 
int<lower=0> nroutes;//number of routes 
int<lower=0> ncovs;//number of route-level covariates 

real slopes[n];//vector of slope estimates
real sd_slopes[n];// vector of observed sd of slope estimates
int species[n];// species indicators
int route[n];// route indicators

real cov_nb[nspecies];// species-level covariates - change in footprint on wintering grounds

matrix[ncovs,nroutes] covr; // route-level predictors
 
 
 
}

// 
parameters {

vector[nroutes] rte_raw;
real<lower=0> sdrte;
real sps[nspecies];
matrix[nspecies,ncovs] b_raw;
real<lower=0> sd_b[ncovs];
real B[ncovs];
real b_nb;
real<lower=0> epsilon;
real true_slopes[n];

}

transformed parameters{

  vector[nroutes] rte;

  real mu_slopes[n];
  matrix[nspecies,nroutes] cov_adj;
  matrix[nspecies,ncovs] b;

 
 rte = sdrte*rte_raw;
 
   for(v in 1:ncovs){
 b[,v] = sd_b[v]*b_raw[,v] + B[v];
   }
   

 cov_adj = b*covr; // matrix multiplication of covariate effects
  //cov_adj[s,r] = b[s,v]*covr[v,r];// result of above
  
  
for(i in 1:n){
mu_slopes[i] = (sps[species[i]] + b_nb*cov_nb[species[i]] + cov_adj[species[i],route[i]] + rte[route[i]]);
}

  
}
// 
model {


// random route effects on trends to account for repeated samples at routes
rte_raw ~ normal(0,1);
sum(rte_raw) ~ normal(0,0.001*nroutes);//sum to zero soft constraint on route level random effects

sdrte ~ normal(0,0.1);

// species mean trends - fixed effects
sps ~ normal(0,0.1);


// species effects of covariates
  for(v in 1:ncovs){
  b_raw[,v] ~ normal(0,0.1); // species random effects of covariates
  sum(b_raw[,v]) ~ normal(0,0.001*nspecies);//sum to zero soft constraint on species level random effects of covariates

  sd_b[v] ~ normal(0,0.1); // covariate level variance across species
  B[v] ~ student_t(4,0,0.05); // covariate hyperparameter - regularizing hyperparameters
  } 

//espilon
epsilon ~ normal(0,0.1);// process error

// non-breeding range footprint effect
b_nb ~ normal(0,0.1);

// likelihood
// process model
   true_slopes ~ normal(mu_slopes,epsilon);
// observation model   
for(i in 1:n){
   slopes[i] ~ normal(true_slopes[i],sd_slopes[i]);
    } 
  
}



