data {
  int<lower=0> N_species; // Number of species
  int<lower=1> N_reactor; // Number of reactors
  int<lower=1> N_timepoint; // Number of timepoints
  int<lower=1> N_conc_measurement; // Number of conc measurements
  int<lower=1> N_biomass_measurement; // Number of biomass measurements
  array[N_species] int<lower=0,upper=1> is_substrate;
  array[N_conc_measurement] int<lower=1,upper=N_reactor> reactor_yconc;
  array[N_conc_measurement] int<lower=1,upper=N_species> species_yconc;
  array[N_conc_measurement] int<lower=1,upper=N_timepoint> timepoint_yconc;
  array[N_biomass_measurement] int<lower=1,upper=N_reactor> reactor_ybiomass;
  array[N_biomass_measurement] int<lower=1,upper=N_timepoint> timepoint_ybiomass;
  vector<lower=0>[N_timepoint] timepoint_time;
  vector<lower=0>[N_conc_measurement] yconc;
  vector<lower=0>[N_biomass_measurement] ybiomass;
  int<lower=0,upper=1> likelihood;
}
transformed data {
  vector[N_species] yield_sign;
  for (s in 1:N_species) yield_sign[s] = is_substrate[s] ? -1.0 : 1.0;
}
parameters {
  vector<lower=0>[N_reactor] biomass_0;
  array[N_reactor] vector<lower=0>[N_species] conc_0;
  // growth rate effects
  real amu_intercept;
  vector[N_reactor] amu_reactor; // Reactor specific growth rate
  real<lower=0> tmu_reactor; // Biological variation of growth rate
  // yield effects
  vector[N_species] ayd_species;
  matrix[N_reactor, N_species] ayd_reactor;
  vector<lower=0>[N_species] tyd_reactor; // Biological variation of growth rate
  // measurement error
  real<lower=0.001, upper=50> sigma_ybiomass;
  vector<lower=0.001, upper=50>[N_species] sigma_yconc;
}
transformed parameters {
  vector[N_reactor] mu = amu_intercept + amu_reactor;
  matrix[N_timepoint, N_reactor] biomass;
  matrix[N_timepoint, N_reactor] dbiomass;
  array[N_reactor] vector[N_species] yield;
  array[N_timepoint, N_reactor] vector[N_species] conc;
  for (tix in 1:N_timepoint){
    for (r in 1:N_reactor){
      biomass[tix, r] = fmin(1e9, fmax(1e-7, biomass_0[r] .* exp(mu[r] * timepoint_time[tix])));
    }
    dbiomass[tix] = biomass[tix] - biomass_0';
  }
  for (r in 1:N_reactor){
    yield[r] = exp(ayd_species + ayd_reactor[r]');
    for (tix in 1:N_timepoint){
      for (s in 1:N_species)
        conc[tix, r][s] = fmin(1e9, fmax(1e-7, conc_0[r][s] + yield_sign[s] .* (yield[r][s] * dbiomass[tix, r])));
    }
  }
}
model {
  biomass_0 ~ lognormal(0, 4);
  amu_intercept ~ normal(0, 1);
  amu_reactor ~ normal(0, tmu_reactor);
  ayd_species ~ normal(0, 3);
  tmu_reactor ~ normal(0, 1);
  tyd_reactor ~ normal(0, 1);
  sigma_yconc ~ lognormal(-1, 1);
  sigma_ybiomass ~ lognormal(-1, 1);
  for (r in 1:N_reactor)
    conc_0[r] ~ lognormal(0, 3);
  for (s in 1:N_species)
    ayd_reactor[,s] ~ normal(0, tyd_reactor[s]);
  if (likelihood){
    for (n in 1:N_conc_measurement){
      real conc_hat = conc[timepoint_yconc[n], reactor_yconc[n], species_yconc[n]];
      yconc[n] ~ lognormal(log(conc_hat), sigma_yconc[species_yconc[n]]);
    }
    for (n in 1:N_biomass_measurement){
      real biomass_hat = biomass[timepoint_ybiomass[n], reactor_ybiomass[n]];
      ybiomass[n] ~ lognormal(log(biomass_hat), sigma_ybiomass);
    }
  }
}
