data {
  int<lower=1> N_reactor; // Number of reactors
  int<lower=1> N_timepoint; // Number of timepoints
  int<lower=1> N_biomass_measurement; // Number of biomass measurements
  array[N_biomass_measurement] int<lower=1,upper=N_reactor> reactor_ybiomass;
  array[N_biomass_measurement] int<lower=1,upper=N_timepoint> timepoint_ybiomass;
  vector<lower=0>[N_timepoint] timepoint_time;
  vector<lower=0>[N_biomass_measurement] ybiomass;
  int<lower=0,upper=1> likelihood;
}
parameters {
  vector<lower=0>[N_reactor] biomass_0;
  // growth rate effects
  real amu_intercept;
  vector[N_reactor] amu_reactor; // Reactor specific growth rate
  real<lower=0> tmu_reactor; // Biological variation of growth rate
  // measurement error
  real<lower=0.001, upper=50> sigma_ybiomass;
}
transformed parameters {
  vector[N_reactor] mu = amu_intercept + amu_reactor;
  matrix[N_timepoint, N_reactor] biomass;
  matrix[N_timepoint, N_reactor] dbiomass;
  for (tix in 1:N_timepoint){
    for (r in 1:N_reactor){
      biomass[tix, r] = fmin(1e9, fmax(1e-7, biomass_0[r] .* exp(mu[r] * timepoint_time[tix])));
    }
    dbiomass[tix] = biomass[tix] - biomass_0';
  }
}
model {
  biomass_0 ~ lognormal(0, 4);
  amu_intercept ~ normal(0, 1);
  amu_reactor ~ normal(0, tmu_reactor);
  tmu_reactor ~ normal(0, 1);
  sigma_ybiomass ~ lognormal(-1, 1);
  if (likelihood){
    for (n in 1:N_biomass_measurement){
      real biomass_hat = biomass[timepoint_ybiomass[n], reactor_ybiomass[n]];
      ybiomass[n] ~ lognormal(log(biomass_hat), sigma_ybiomass);
    }
  }
}

