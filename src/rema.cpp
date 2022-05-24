// REMA

// Model author: Pete Hulson
// Adapted to TMB by Jane Sullivan in May 2022
// Hulson, P.-J. F., K. B. Echave, P. D. Spencer, and J. N. Ianelli. 2021. Using multiple Indices for
// biomass and apportionment estimation of Alaska groundfish stocks. U.S. Dep. Commer., NOAA
// Tech. Memo. NMFS-AFSC-414, 28 p.

// Code summary: A multivariate random effects model that accepts an additional
// relative (i.e. cpue) survey index. If the multi-survey mode is off, REMA runs
// the same as the univariate (RE) and multivariate (REM) versions of the randon
// effects model. The biomass and cpue surveys can have multiple strata, which
// can be uniquely defined by survey. Because the cpue survey is intended to inform the
// trend of the biomass, the number of cpue strata cannot exceed the number of
// biomass strata.


#define TMB_LIB_INIT R_init_wham
#include <TMB.hpp>
#include <iostream>

#define see(object) std::cout << #object ":\n" << object << "\n";

template <class Type>
Type square(Type x){return x * x;}

template<class Type>
bool isNA(Type x){return R_IsNA(asDouble(x));}

template<class Type>
Type objective_function<Type>::operator() ()
{
  // 0 = no, 1 = yes (estimate scaling parameters)
  DATA_INTEGER(multi_survey_mode);

  // model dimensions
  DATA_IVECTOR(model_yrs);

  // survey biomass (absolute biomass gives scale to the model)
  // DATA_IVECTOR(biomass_yrs);
  DATA_MATRIX(biomass_obs);
  DATA_MATRIX(biomass_cv);
  DATA_SCALAR(wt_biomass);
  DATA_IVECTOR(pointer_PE_biomass); // length = ncol biomass obs (# strata), unique values = number of PE parameters
  DATA_IVECTOR(pointer_q_biomass);  // length = ncol biomass obs (# strata), unique values = number of q parameters

  // survey cpue (relative index that can inform trend in the model)
  // DATA_IVECTOR(cpue_yrs);
  DATA_MATRIX(cpue_obs);
  DATA_MATRIX(cpue_cv);
  DATA_SCALAR(wt_cpue);
  // DATA_IVECTOR(pointer_PE_cpue);
  DATA_IVECTOR(pointer_q_cpue);

  // parameter section
  PARAMETER(dummy);  // dummy var for troubleshooting
  PARAMETER_VECTOR(log_PE);
  PARAMETER_VECTOR(log_q);
  PARAMETER_MATRIX(log_biomass_pred);

  // negative log likelihood
  vector<Type> jnll(3); // random walk, biomass obs, cpue obs
  jnll.setZero();
  Type nll = 0;

  // derived quantities - model dimensions
  int nyrs = model_yrs.size();
  int n_strata_biomass = biomass_obs.cols();
  int n_strata_cpue = cpue_obs.cols();

  // // model predictions
  matrix<Type> biomass_pred(nyrs,n_strata_biomass);
  matrix<Type> cpue_pred(nyrs,n_strata_cpue);
  // predicted biomass summed to the same strata level as the cpue strata
  matrix<Type> biomass_pred_cpue_strata(nyrs, n_strata_cpue);
  biomass_pred_cpue_strata.setZero();

  // derived quantities - biomass survey observations
  matrix<Type> log_biomass_obs(nyrs,n_strata_biomass);
  log_biomass_obs = log(biomass_obs.array());
  matrix<Type> biomass_sd(nyrs,n_strata_biomass);
  biomass_sd = biomass_cv.array() * biomass_cv.array() + Type(1.0);
  biomass_sd = sqrt(log(biomass_sd.array()));

  // derived quantities - cpue survey observations
  matrix<Type> log_cpue_obs(nyrs,n_strata_cpue);
  log_cpue_obs = log(cpue_obs.array());
  matrix<Type> cpue_sd(nyrs,n_strata_cpue);
  cpue_sd = cpue_cv.array() * cpue_cv.array() + Type(1.0);
  cpue_sd = sqrt(log(cpue_sd.array()));
  matrix<Type> log_cpue_pred(nyrs,n_strata_cpue);
  log_cpue_pred.setZero();

  // process errors
  vector<Type> PE_var(log_PE.size());
  PE_var = exp(Type(2.0) * log_PE.array()); // FLAG - should this by x^2 not 2*x?

  // random effects contribution to likelihood
  for(int i = 1; i < nyrs; i++) {
    for(int j = 0; j < n_strata_biomass; j++) {
      jnll(0) -= dnorm(log_biomass_pred(i-1,j), log_biomass_pred(i,j), exp(log_PE(pointer_PE_biomass(j))), true);
    }
  }
  biomass_pred = exp(log_biomass_pred.array());

  // data likelihood for biomass survey observations
  for(int i = 0; i < nyrs; i++) {
    for(int j = 0; j < n_strata_biomass; j++) {
      if(biomass_obs(i,j) >= 0) {
        jnll(1) -= dnorm(log_biomass_pred(i,j), log(biomass_obs(i,j) + 0.0001), biomass_sd(i,j), true);
      }
    }
  }
  jnll(1) = jnll(1) * wt_biomass;

  // If in multi-survey mode (1=on, 0=off), calculate predicted cpue and data
  // likelihood for cpue survey observations
  if(multi_survey_mode == 1) {

    // get predicted biomass at the strata level for cpue (i.e. account for the
    // scenario when you may have biomass at a higher resolution than cpue; e.g.
    // 4 biomass strata that are represented by 2 cpue strata; note: currently
    // you cannot have a scenario where cpue has more strata than biomass)
    for(int i = 0; i < nyrs; i++) {
      for(int j = 0; j < n_strata_biomass; j++) {
        biomass_pred_cpue_strata(i,pointer_q_biomass(j)) += exp(log_biomass_pred(i,j) + 0.0001);
      }

      // get predicted cpue
      for(int j = 0; j < n_strata_cpue; j++) {
        cpue_pred(i,j) = exp(log_q(pointer_q_cpue(j))) * biomass_pred_cpue_strata(i,j);
      }

      // get data likelihood for cpue survey
      for(int j = 0; j < n_strata_cpue; j++) {
        if(cpue_obs(i,j) >= 0) {
          jnll(2) -= dnorm(log(cpue_pred(i,j)), log(cpue_obs(i,j) + 0.0001), cpue_sd(i,j), true);
        }
      }
    }
    log_cpue_pred = log(cpue_pred.array());
    jnll(2) = jnll(2) * wt_cpue;
  }

  vector<Type> tot_biomass_pred;
  tot_biomass_pred = biomass_pred.rowwise().sum();

  vector<Type> tot_cpue_pred;
  tot_cpue_pred = cpue_pred.rowwise().sum();

  ADREPORT(biomass_pred);
  ADREPORT(tot_biomass_pred);
  ADREPORT(cpue_pred);
  ADREPORT(tot_cpue_pred);

  REPORT(nll);

  // jnll = dummy * dummy;        // Uncomment when debugging code
  nll = jnll.sum();
  return nll;

}
