// REMA

// A multivariate random effects model that accepts an additional relative (i.e.
// cpue) survey index. If the multi-survey mode is off, REMA runs the same as
// the univariate (RE) and multivariate (REM) versions of the randon effects
// model. The biomass and cpue surveys can have multiple strata, which can be
// uniquely defined by survey. Because the cpue survey is intended to inform the
// trend of the biomass, the number of cpue strata cannot exceed the number of
// biomass strata.

// Hulson, P.-J. F., K. B. Echave, P. D. Spencer, and J. N. Ianelli. 2021. Using multiple Indices for
// biomass and apportionment estimation of Alaska groundfish stocks. U.S. Dep. Commer., NOAA
// Tech. Memo. NMFS-AFSC-414, 28 p.

#define TMB_LIB_INIT R_init_rema
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
  // 0 = no, 1 = yes; include cpue survey index and estimate scaling
  // parameters?
  DATA_INTEGER(multi_survey);

  // model dimensions
  DATA_IVECTOR(model_yrs);

  // survey biomass (absolute biomass gives scale to the model)
  DATA_MATRIX(biomass_obs);
  DATA_MATRIX(biomass_cv);

  // survey cpue (relative index that can inform trend in the model)
  DATA_MATRIX(cpue_obs);
  DATA_MATRIX(cpue_cv);
  DATA_INTEGER(sum_cpue_index); // 0 = no, 1 = yes (default = 0); can the cpue index be summed across strata (e.g. Relative population numbers = 1, Numbers per hachi = 0)

  // switch for observation error distribution
  DATA_INTEGER(obs_error_type); // 0 = default = lognormal (make data assumptions about zeros); 1 = Tweedie (model zeros)

  // indexing vectors
  DATA_IVECTOR(pointer_PE_biomass); // length = ncol biomass obs (# strata), unique values = number of PE parameters
  DATA_IVECTOR(pointer_biomass_cpue_strata);  // length = ncol biomass obs (# strata), unique values = number of cpue survey strata
  DATA_IVECTOR(pointer_q_cpue); // length = ncol cpue obs (# strata), unique values = number of q parameters

  DATA_SCALAR(wt_biomass); // weight for biomass data likelihood (default = 1)
  DATA_SCALAR(wt_cpue); // weight for cpue data likelihood (default = 1)

  // process error penalty/prior options
  DATA_INTEGER(PE_penalty_type); // 0 = "none", 1 = "wt", 2 = "squared_penalty", 3 = "normal_prior"
  DATA_SCALAR(wt_PE); // weight for the random effects and process error likelihood (default = 1)
  DATA_VECTOR(squared_penalty_PE); // prevents process error from shrinking to zero
  DATA_VECTOR(pmu_log_PE); // normal prior
  DATA_VECTOR(psig_log_PE);
  // scaling parameter penalty/prior options and likelihood weight for CPUE data
  // likelihood
  DATA_INTEGER(q_penalty_type); // 0 = "none", 1 = "normal_prior"
  DATA_VECTOR(pmu_log_q); // normal prior
  DATA_VECTOR(psig_log_q);

  // parameter section

  // PARAMETER(dummy);  // dummy var for troubleshooting
  PARAMETER_VECTOR(log_PE); // process errors
  PARAMETER_VECTOR(log_q); // scaling parameters

  PARAMETER_VECTOR(log_tweedie_dispersion); // tweedie dispersion (one for biomass survey and optional one for cpue survey)
  PARAMETER_VECTOR(logit_tweedie_p); // tweedie power parameter (one for biomass survey and optional one for cpue survey)

  vector<Type> tweedie_dispersion = exp(log_tweedie_dispersion);
  vector<Type> tweedie_p = invlogit(logit_tweedie_p) + Type(1.0);

  PARAMETER_MATRIX(log_biomass_pred); // random effects of predicted biomass

  // negative log likelihood
  vector<Type> jnll(3); // random walk, biomass obs, cpue obs
  jnll.setZero();
  Type nll = 0;

  // derived quantities - model dimensions
  int nyrs = model_yrs.size();
  int n_strata_biomass = biomass_obs.cols();
  int n_strata_cpue = cpue_obs.cols();

  // model predictions

  // predicted biomass and cpue on the natural and log scale
  matrix<Type> biomass_pred(nyrs, n_strata_biomass);
  biomass_pred.setZero();
  matrix<Type> cpue_pred(nyrs, n_strata_cpue);
  cpue_pred.setZero();
  matrix<Type> log_cpue_pred(nyrs, n_strata_cpue);
  log_cpue_pred.setZero();

  // predicted biomass summed to the same strata level as the cpue strata
  // (useful when n_strata_biomass > n_strata_cpue)
  matrix<Type> biomass_pred_cpue_strata(nyrs, n_strata_cpue);
  biomass_pred_cpue_strata.setZero();
  matrix<Type> log_biomass_pred_cpue_strata(nyrs, n_strata_cpue);
  log_biomass_pred_cpue_strata.setZero();

  // derived quantities - biomass survey observations
  matrix<Type> log_biomass_obs(nyrs, n_strata_biomass);
  log_biomass_obs = log(biomass_obs.array());
  matrix<Type> log_biomass_sd(nyrs, n_strata_biomass);
  log_biomass_sd = biomass_cv.array() * biomass_cv.array() + Type(1.0);
  log_biomass_sd = sqrt(log(log_biomass_sd.array()));
  matrix<Type> biomass_sd(nyrs, n_strata_biomass);
  biomass_sd = biomass_obs.array() * biomass_cv.array();
  matrix<Type> biomass_dispersion(nyrs, n_strata_biomass);
  biomass_dispersion.setZero();

  // derived quantities - cpue survey observations
  matrix<Type> log_cpue_obs(nyrs, n_strata_cpue);
  log_cpue_obs = log(cpue_obs.array());
  matrix<Type> log_cpue_sd(nyrs, n_strata_cpue);
  log_cpue_sd = cpue_cv.array() * cpue_cv.array() + Type(1.0);
  log_cpue_sd = sqrt(log(log_cpue_sd.array()));
  matrix<Type> cpue_sd(nyrs, n_strata_cpue);
  cpue_sd = cpue_obs.array() * cpue_cv.array();

  // random effects contribution to likelihood
  for(int i = 1; i < nyrs; i++) {
    for(int j = 0; j < n_strata_biomass; j++) {
      jnll(0) -= dnorm(log_biomass_pred(i-1,j), log_biomass_pred(i,j), exp(log_PE(pointer_PE_biomass(j))), 1);
    }
  }
  switch(PE_penalty_type) {

  case 0: // no penalty
    jnll(0) = jnll(0);
    break;

  case 1: // weight
    jnll(0) = wt_PE * jnll(0);
    break;

  case 2: // squared penalty
    for(int k = 0; k < log_PE.size(); k++) {
      jnll(0) += pow(log_PE(k) + squared_penalty_PE(k), 2);
    }
    break;

  case 3: // normal prior
    for(int k = 0; k < log_PE.size(); k++) {
      jnll(0) -= dnorm(log_PE(k), pmu_log_PE(k), psig_log_PE(k), 1);
    }
    break;
  }

  // likelihood for observation error of biomass survey data
  switch(obs_error_type) {

  case 0: // lognormal
    for(int i = 0; i < nyrs; i++) {
      for(int j = 0; j < n_strata_biomass; j++) {
        if(biomass_obs(i,j) > 0) {
          jnll(1) -= dnorm(log_biomass_pred(i,j), log_biomass_obs(i,j), log_biomass_sd(i,j), 1);
        }
        biomass_pred(i,j) = exp(log_biomass_pred(i,j));
      }
    }
    break;

  case 1: // Tweedie

    Type tmp_biomass_sd = 100.0;

    for(int i = 0; i < nyrs; i++) {
      for(int j = 0; j < n_strata_biomass; j++) {

        biomass_pred(i,j) = exp(log_biomass_pred(i,j));

        // if(biomass_sd(i,j) == 0) {
        //   tmp_biomass_sd = 100.0;
        // }
        // if(biomass_sd(i,j) > 0) {
        //   tmp_biomass_sd = biomass_sd(i,j);
        // }

        if(biomass_obs(i,j) == 0) {
          biomass_dispersion(i,j) = Type(100.0);
        }
        if(biomass_obs(i,j) > 0) {
          biomass_dispersion(i,j) = (biomass_sd(i,j) * biomass_sd(i,j)) / pow(biomass_pred(i,j), tweedie_p(0));
        }

        if(biomass_obs(i,j) >= 0) {
          // tweedie_dispersion(0) = (tmp_biomass_sd * tmp_biomass_sd) / pow(biomass_pred(i,j), tweedie_p(0));
          jnll(1) -= dtweedie(biomass_pred(i,j), biomass_obs(i,j), tweedie_dispersion(0), tweedie_p(0), 1);
          // jnll(1) -= dtweedie(biomass_pred(i,j), biomass_obs(i,j), biomass_dispersion(i,j), tweedie_p(0), 1);
        }

      }
    }
    break;
  }

  jnll(1) = jnll(1) * wt_biomass;

  // If in multi-survey mode (1=on, 0=off), calculate predicted cpue and data
  // likelihood for cpue survey observations
  if(multi_survey == 1) {

    // get predicted biomass at the strata level for cpue (i.e. account for the
    // scenario when you may have biomass at a higher resolution than cpue; e.g.
    // 4 biomass strata that are represented by 2 cpue strata; note: currently
    // you cannot have a scenario where survey cpue has more strata than biomass
    // survey)
    for(int i = 0; i < nyrs; i++) {

      for(int j = 0; j < n_strata_biomass; j++) {
        biomass_pred_cpue_strata(i,pointer_biomass_cpue_strata(j)) += biomass_pred(i,j);
      }

      // get predicted cpue, log-transform for variance estimates
      for(int j = 0; j < n_strata_cpue; j++) {
        cpue_pred(i,j) = exp(log_q(pointer_q_cpue(j))) * biomass_pred_cpue_strata(i,j);
        log_cpue_pred(i,j) = log(cpue_pred(i,j));
        log_biomass_pred_cpue_strata(i,j) = log(biomass_pred_cpue_strata(i,j));
      }
    }


    // get data likelihood for cpue survey
    switch(obs_error_type) {

    case 0: // lognormal

      for(int i = 0; i < nyrs; i++) {
        for(int j = 0; j < n_strata_cpue; j++) {
          if(cpue_obs(i,j) > 0) {
            jnll(2) -= dnorm(log_cpue_pred(i,j), log_cpue_obs(i,j), log_cpue_sd(i,j), 1);
          }
        }
      }
      break;

    case 1: // Tweedie
      for(int i = 0; i < nyrs; i++) {
        for(int j = 0; j < n_strata_cpue; j++) {

          if(cpue_obs(i,j) >= 0) {
            jnll(2) -= dtweedie(cpue_pred(i,j), cpue_obs(i,j), tweedie_dispersion(1), tweedie_p(1), 1);
            // tweedie_dispersion(1) =  (cpue_sd(i,j) * cpue_sd(i,j)) / pow(cpue_pred(i,j), tweedie_p(1));

          }
        }
      }
      break;
    }

    jnll(2) = jnll(2) * wt_cpue;

    // optional penalty or prior on log_q
    switch(q_penalty_type) {

    case 0: // no penalty
      jnll(2) = jnll(2);
      break;

    case 1: // normal prior
      for(int k = 0; k < log_q.size(); k++) {
        jnll(2) -= dnorm(log_q(k), pmu_log_q(k), psig_log_q(k), 1);
      }
      break;
    }
  }

  // report section
  ADREPORT(log_biomass_pred);

  if(n_strata_biomass > 1) {
    vector<Type> tot_biomass_pred;
    tot_biomass_pred = biomass_pred.rowwise().sum();
    vector<Type> log_tot_biomass_pred;
    log_tot_biomass_pred = log(tot_biomass_pred);
    ADREPORT(log_tot_biomass_pred);
  }

  if(multi_survey == 1) {
    ADREPORT(log_cpue_pred);

    // only sum cpue index if its appropriate for that index (e.g. appropriate
    // for relative popn numbers, not appropriate for nominal cpue).
    if(sum_cpue_index == 1 & n_strata_cpue > 1){

      vector<Type> tot_cpue_pred;
      tot_cpue_pred = cpue_pred.rowwise().sum();
      vector<Type> log_tot_cpue_pred;
      log_tot_cpue_pred = log(tot_cpue_pred);
      REPORT(log_tot_cpue_pred);
      ADREPORT(log_tot_cpue_pred);

    }
    // report biomass at the cpue strata level (and get SDs) if they have
    // different strata definitions
    if(n_strata_biomass > n_strata_cpue) {
      REPORT(log_biomass_pred_cpue_strata);
      ADREPORT(log_biomass_pred_cpue_strata);
    }
  }

  REPORT(log_biomass_pred);
  // REPORT(log_cpue_pred);
  REPORT(tweedie_dispersion);
  REPORT(tweedie_p);
  REPORT(biomass_pred);
  REPORT(biomass_sd);
  REPORT(jnll);

  // jnll = dummy * dummy;        // Uncomment when debugging code
  nll = jnll.sum();
  return nll;

}
