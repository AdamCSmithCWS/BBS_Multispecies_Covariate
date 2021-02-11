//
//

// 
data {


// int<lower=0,upper=1> I[nspecies,nroutes]; // indicator matrix for species by route inclusions

int<lower=0> n;//number of slope estimates 
int<lower=0> nspecies;//number of species 
int<lower=0> nroutes;//number of routes 
int<lower=0> ncovs;//number of route-level covariates 
real  nB; //expected number of non-zero effects
 

real slopes[n];//vector of slope estimates
real sd_slopes[n];// vector of observed sd of slope estimates
int species[n];// species indicators
int route[n];// route indicators

real cov_nb[nspecies];// species-level covariates - change in footprint on wintering grounds

matrix[ncovs,nroutes] covr; // route-level predictors
 
 
 
}

transformed data {
 real slab_scale = 0.005; //scale for large slopes = ~0.5% trend change for 1 SD of predictor
  real slab_scale2 = square(slab_scale); //
  real slab_df = 25; //effective df for t-distribution on large slopes
  real half_slab_df = 0.5 * slab_df;
}

// 
parameters {

vector[nroutes] rte_raw;
real<lower=0> sdrte;
real sps[nspecies];
real<lower=0> sd_b[ncovs];
real b_nb;
real<lower=0> epsilon;
real true_slopes[n];

matrix[nspecies,ncovs] b_raw;

// Finnish horeshoe parameter estimates

vector[ncovs] B_tilde;
vector<lower=0>[ncovs] lambda; 
real<lower=0> c2_tilde;
real<lower=0> tau_tilde;

}

transformed parameters{

vector[ncovs] B;
  vector[nroutes] rte;

  real mu_slopes[n];
  matrix[nspecies,nroutes] cov_adj;
  matrix[nspecies,ncovs] b;

  {
    real tau0 = (nB / (ncovs - nB)) * (epsilon / sqrt(1.0 * nroutes));//adapted from original to use nroutes instead of n
    real tau = tau0 * tau_tilde; // tau ~ cauchy(0, tau0)

    // c2 ~ inv_gamma(half_slab_df, half_slab_df * slab_scale2)
    // Implies that marginally B ~ student_t(slab_df, 0, slab_scale)
    real c2 = slab_scale2 * c2_tilde;

    vector[ncovs] lambda_tilde =
      sqrt( c2 * square(lambda) ./ (c2 + square(tau) * square(lambda)) );

    // B ~ normal(0, tau * lambda_tilde)
    B = tau * lambda_tilde .* B_tilde;
  }


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

//Finnish Horseshoe
B_tilde ~ normal(0,1);
lambda ~ cauchy(0,1);
tau_tilde ~ cauchy(0,1);
c2_tilde ~ inv_gamma(half_slab_df,half_slab_df);

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
  //B[v] ~ student_t(4,0,0.05); // covariate hyperparameter - regularizing hyperparameters
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



