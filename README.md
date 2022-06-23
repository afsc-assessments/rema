# REMA
A generalized **R**andom **E**ffects model for fitting survey biomass estimates with the option of including **M**ultiple survey strata and an **A**dditional survey index

The random effects (RE) model was deveped by the North Pacific Fisheries Management Council (NPFMC) Groundfish Plan Team's Survey Averaging working group and has been used at the Alaska Fisheries Science Center (AFSC) since [2013](https://github.com/afsc-assessments/SurveyAverageRandomEffects/blob/013c9a937fa0133f594c7d66248677685ae77010/code/re.tpl) to estimate biomass in data-limited groundfish and crab stock assessments, and to apportion estimates of Acceptable Biological Catch by area. The RE model uses a Kalman filter approach, where the biomass is estimated as a series of random effects and the underlying state dynamics are modeled as a random walk ([Oct 2013 Joint GPT minutes](https://meetings.npfmc.org/CommentReview/DownloadFile?p=11009549-068b-40cf-903d-67f90686db60.pdf&fileName=C4%20c1%20Joint%20Plan%20Team%20Minutes.pdf)). There are several versions of the original single-strata (univariate) RE and the single-survey, multiple-strata (multivariate; REM) models currently in use at the AFSC. Though the versions share the same underlying state-space dynamics, [Monnahan et al. (2021)](https://meetings.npfmc.org/CommentReview/DownloadFile?p=86098951-a0ed-4021-a4e1-95abe5a357fe.pdf&fileName=Tiers%204%20and%205%20assessment%20considerations.pdf) found inconsistencies in the treatment of zero biomass observations and use of a prior or penalty on the process error parameter. Additionally, [Hulson et al. (2021)](https://repository.library.noaa.gov/view/noaa/28174) developed a third version of the model (REMA) that includes a second relative survey index (e.g., NMFS longline survey Relative Population Numbers) and estimates additional scaling parameters. The purpose of the `rema` R library is to develop a unified random effects model that is flexible enough to accommodate the multitude of use cases at the AFSC. As presented here, the REMA model is generalized to include the RE and REM. The variable names have been updated from their original versions to increase interpretability and transparency, and the model has been rewritten in Template Model Builder (TMB; [Kristensen et al. 2016](https://www.jstatsoft.org/article/view/v070i05)).

The structure, naming conventions, functions, and documentation in this library were inspired by and modeled after the Woods Hole Assessment Model ([WHAM](https://timjmiller.github.io/wham/); Miller and Stock 2020), an open-source, state-space age-structured assessment model and R package.

This library is under development and has not been vetted for use in stock assessments.

## Installation

Use the `devtools` package to install the `rema` package from github. If you do not have `devtools` installed, you must do that first.

```
# install.packages("devtools")
devtools::install_github("JaneSullivan-NOAA/rema", dependencies = TRUE)

# Example R scripts are downloaded when `rema` is installed. Locate them on your computer by running the following commands:
(rema_path <- find.package('rema'))
(rema_examples <- file.path(rema_path, 'example_scripts'))
list.files(rema_examples)

```

## The `rema` worflow

1.  Load `rema` and data. The user can read biomass or other abundance index data from file (e.g. csv files), or they can use the `rwout.rep` report file from the ADMB version of the RE model using `read_admb_re()`.

2.  Define the model structure and assumptions using `prepare_rema_input()`. This function allows users to quickly transition from a single to two survey model, specify alternative process error structures, add likelihood penalties or priors on parameters, and evaluate alternative assumptions about zero biomass observations.

3.  Fit the specified REMA model using `fit_rema()` and determine whether the model has met basic convergence criteria using `check_convergence()`. 

4.  Extract `rema` model output into clean, consistently formatted data frames using `tidy_rema()`. The user can visualize this model output using `plot_rema()`, or quickly format it into tables for a report.

5.  Compare alternative REMA models and conduct model selection using `compare_rema_models()`. Output from this function includes a table of Akaike Information Criteria (AIC) when appropriate, figures, and tidied data frames. This function also accepts model output from the ADMB version of the RE model for easy comparison to past models.

Taken together, these functions allow R users to quickly fit and interrogate a suite of simple statistical models in TMB without needing software-specific training or expertise.

![standards](https://user-images.githubusercontent.com/68303089/174915022-e9b89cb4-79ed-4dcc-af93-3cd199b5f82c.png)
