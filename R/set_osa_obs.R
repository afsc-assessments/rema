# Add data objects needed to calculate OSA residuals. Code modified from
# https://github.com/timjmiller/wham/blob/28cab35ede66506dcaa73f7de8fc0334870893fc/R/set_osa_obs.R

set_osa_obs <- function(input) {

  data = input$data

  obs_colnames <- c('year', 'survey', 'strata', 'value')
  obs <- data.frame(matrix(ncol = length(obs_colnames), nrow = 0))
  colnames(obs) <- obs_colnames

  # 1) observations for biomass survey
  x <- data$biomass_obs
  for(i in 1:length(data$model_yrs)) {
    for(j in 1:ncol(data$biomass_obs)) {

      # obs_error_type = 0 = lognormal, take log, otherwise obs_error_type =
      # tweedie; unname takes away stratum name
      if(!is.na(x[i,j])) {
        if(data$obs_error_type == 0) {
          val <- unname(log(x[i,j]))
        } else {
          val <- unname(x[i,j])
        }

        tmp = data.frame(year = data$model_yrs[i], survey = 'Biomass survey', strata = colnames(x)[j], value = val)
        obs <- rbind(obs, tmp[, obs_colnames])
      }
    }
  }

  # 2) observations for cpue survey
  x <- data$cpue_obs
  for(i in 1:length(data$model_yrs)) {
    for(j in 1:ncol(data$cpue_obs)) {

      # obs_error_type = 0 = lognormal, take log, otherwise obs_error_type =
      # tweedie; unname takes away stratum name
      if(!is.na(x[i,j])) {
        if(data$obs_error_type == 0) {
          val <- unname(log(x[i,j]))
        } else {
          val <- unname(x[i,j])
        }

        tmp = data.frame(year = data$model_yrs[i], survey = 'CPUE survey', strata = colnames(x)[j], value = val)
        obs <- rbind(obs, tmp[, obs_colnames])
      }
    }
  }

  # index for each observation, starting at 0 for TMB
  obs$ind <- 1:dim(obs)[1] - 1

  # set up data-specific keep matrices that will hold the index for each
  # observations in the final observation vector

  keep_biomass_obs <- matrix(NA, nrow = nrow(data$biomass_obs), ncol = ncol(data$biomass_obs))
  for(i in 1:nrow(data$biomass_obs)) {
    for(j in 1:ncol(data$biomass_obs)) {
      if(!is.na(data$biomass_obs[i,j])) {
        keep_biomass_obs[i,j] <- subset(obs, year == data$model_yrs[i] & survey == 'Biomass survey' & strata == colnames(data$biomass_obs)[j])$ind
      }
    }
  }

  keep_cpue_obs <- matrix(NA, nrow = nrow(data$cpue_obs), ncol = ncol(data$cpue_obs))
  for(i in 1:nrow(data$cpue_obs)) {
    for(j in 1:ncol(data$cpue_obs)) {
      if(!is.na(data$cpue_obs[i,j])) {
        keep_cpue_obs[i,j] <- subset(obs, year == data$model_yrs[i] & survey == 'CPUE survey' & strata == colnames(data$cpue_obs)[j])$ind
      }
    }
  }

  data$obsvec <- obs$value
  data$keep_biomass_obs <- keep_biomass_obs
  data$keep_cpue_obs <- keep_cpue_obs

  input$data <- data
  input$osa <- obs
  return(input)

}

