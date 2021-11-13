library(argparse)
# TODO: Vectorize computations in this function (current inherited version is inefficient)
calculate_ascvd_risk <- function(X){
  # Estimate ASCVD risk conditioned on relevant health information
  # 
  # Check against https://tools.acc.org/ldl/ascvd_risk_estimator/index.html#!/calulate/estimator/ as a reference
  # and to get an intuition for the model's inputs and behavior.
  # 
  #   Args:
  #       X:  An [n_samples x 17] numpy array in which each row represents a single research subject and each
  #           column represents a particular health-related variable. Columns are as follows:
  #               Column 0 represents age,
  #               Column 1 represents gender (1 is female, 0 is male)
  #               Column 2 represents race (1 is black, 0 is white)
  #               Column 4 is systolic blood pressure (mm Hg)
  #               Column 4 is the HDL-cholesterol (mg/dL)
  #               Column 7 is smoking status (1 if smoker, 0 otherwise)
  #               Column 12 is total cholesterol (mg/dL)
  #               Column 16 is the subjects' diabetes status (1 if px is diabetic, 0 otherwise)
  # (for reference, this column ordering is aligned with SPRINT/ACCORD-BP preprocessing as given in
  #   https://github.com/som-shahlab/risk-vs-hte/tree/master/data)
  # Returns:
  #   Calculated 10-year ASCVD Risk based on the "ASCVD Risk Estimator" (for pxs not on hypertension treatment)
  # 
  # To use:
  #   >>> # 55yo female, black, with SBP 120 mm Hg, current smoker, HDL-Cholesterol 80 mm Hg, Total cholesterol 150 mm Hg
  #   >>> X = np.array([[55, 1, 1, 1, 120, 80, 5, 1, 1, 1, 1, 150, 150, 80, 80, 25, 0]])
  # >>> calculate_ascvd_risk(X)
  # array([0.02735308])
  N <- dim(X)[1]
  risks <- rep(0,N)
  coeffs <- list(c(12.344, 11.853, -2.664, -7.99, 1.769, 1.764, 7.837,-1.795, 0.658),
                 c(2.469, 0.302, 0, -0.307, 0, 1.809, 0.549, 0, 0.645),
                 c(-29.799, 4.884, 13.540, -3.114, -13.578, 3.149, 1.957, 0, 7.574, -1.665, 0.661),
                 c(17.114, 0, 0.940, 0, -18.920, 4.475, 27.820, -6.087, 0.691, 0, 0.874))
  
  for (i in 1:N){
    
    female <- ifelse(X[i,2]==1,1,0)
    black <- ifelse(X[i,3]==1,1,0)
    age <- log(X[i,1])
    age_sq <- age ** 2
    chol <- log(X[i,13])
    age_chol <- age * chol
    hdl <- log(X[i,14])
    age_hdl <- age * hdl
    sbp <- log(X[i,5])
    age_sbp <- age * sbp
    smk <- X[i,8]
    age_smk <- age * smk
    dia <- X[i,17]
    
    # how to deal with missing values in female & black??
    if(female==0 & black==0){
      vec <- c(age, chol, age_chol, hdl, age_hdl, sbp,smk, age_smk, dia)
      risks[i] <- 1-0.9144**exp(vec%*%coeffs[[1]]-61.18)
    }
    if(female==0 & black==1){
      vec <- c(age, chol, age_chol, hdl,age_hdl, sbp,smk, age_smk, dia)
      risks[i] <- 1-0.8954**exp(vec%*%coeffs[[2]]-19.54)
    }
    if(female==1 & black==0){
      vec <- c(age, age_sq, chol, age_chol, hdl, age_hdl, sbp, age_sbp, smk, age_smk, dia)
      risks[i] <- 1-0.9665**exp(vec%*%coeffs[[3]]+29.18)
    }
    if(female==1 & black==1){
      vec <- c(age, age_sq, chol, age_chol, hdl, age_hdl, sbp, age_sbp, smk, age_smk, dia)
      risks[i] <- 1-0.9533**exp(vec%*%coeffs[[4]]-86.61)
    }
  }
  return(risks)
}


estimate_rmst <- function(S, unique_times, end_time){
  # Estimate the Restricted Mean Survival Time using survival curves and a pre-specified end time
  # 
  # Args:
  #   S: A [num_subjects x num_unique_times] array representing the survival at time T for each unique time
  # unique_times: An iterable containing the list of unique times, corresponding to columns in S
  # end_time: The end time of the study (need not be [in fact, should not be] the max of the unique times)
  # 
  # Returns:
  #   The Restricted Mean Survival Time for the given survival curve and end time, in the same units as unique_times
  
  col_indx <- length(which(unique_times<=end_time))
  S <- S[,1:col_indx]
  unique_times <- unique_times[unique_times <= end_time]
  time_diffs <- rep(NA,(length(unique_times)-1))
  for (i in 1:(length(unique_times)-1)){
    time_diffs[i] <- time_diffs[i+1]-time_diffs[i]
  }
  time_diffs <- c(unique_times[1], time_diffs) 
  rmst <- S%*%time_diffs
  return(rmst)
}

estimate_treatment_effect_in_rmst <- function(po_sc_control, po_sc_treat, unique_times, end_time){
  # Estimate the treatment effect in terms of the Restricted Mean Survival Time
  # 
  # Args:
  #   po_sc_control: A [num_subjects x num_unique_times] array representing subject survival at time T,
  # for each unique time, assuming subject is assigned to control arm (i.e., not treated)
  # po_sc_treat: A [num_subjects x num_unique_times] array representing subject survival at time T,
  # for each unique time, assuming subject is assigned to treatment arm
  # unique_times: An array with shape [num_unique_times] containing the times that correspond to the columns
  # of po_sc_control and po_sc_treat (which we require to have the same unique times, in this case)
  # end_time: The end time of the study (need not be [in fact, should not be] the max of the unique times)
  # 
  # Returns:
  #   The difference in potential outcomes (in terms of Restricted Mean Survival Time) between treatment and control.
  
  rmst_under_control   <- estimate_rmst(S=po_sc_control, unique_times=unique_times, end_time=end_time)
  rmst_under_treatment <- estimate_rmst(S=po_sc_treat,   unique_times=unique_times, end_time=end_time)
  tau <- rmst_under_treatment - rmst_under_control
  return(tau)
}

round_to_nearest_base <- function(arr, base){
  # Round each element of an array to the nearest base (e.g., [2.4, 2.6, 5.1] w/ base 5 returns [0, 5, 5])
  return (round(arr / base,0) * base)
}

sample_age <- function(num_samples){
  # Sample subject ages from a truncated normal distribution whose first, second moments match the SPRINT RCT data
  sample_ages <- NULL
  repeat{
    tmp_age <- 68+11*rnorm(1,0,1)
    if(tmp_age >= 48 & tmp_age <= 90){
      sample_ages <- c(sample_ages,tmp_age)
    }else{
      sample_ages <- sample_ages
    }
    if(length(sample_ages)==num_samples){
      break
    }
  }
  return(sample_ages)
}

sample_gender <- function(df){
  # Sample subject gender from a distribution whose first and second moments match SPRINT
  prob_is_female <- 0.19731537 + 0.00230716 * df$Age
  return (as.numeric(runif(dim(df)[1]) < prob_is_female))  # 1 if female, 0 otherwise
}

sample_race <- function(df){
  # Sample subject race (black vs. white) using conditional bernoulli distribution matching SPRINT RCT
  prob_is_black <- 1.1792881 - 0.01355267 * df$Age + 0.15666015 * df$Female
  return (as.numeric(runif(dim(df)[1]) < prob_is_black))  # 1 if black, 0 otherwise
}

sample_hispanic <- function(df){
  # Sample ethnicity (non-hispanic vs. hispanic) using conditional bernoulli distribution matching SPRINT RCT
  # (by self-reported race/ethnicity, px may also be white/black/other if hisp)
  prob_is_hispanic <- 0.45825168 - 0.00497156 * df$Age + 0.07282706 * df$Female - -0.12825677 * df$Black
  return (as.numeric(runif(dim(df)[1]) < prob_is_hispanic))  # 1 if hispanic, 0 otherwise
}

sample_sbp <- function(df){
  # Sample subject baseline (seated) systolic blood pressure using conditional normal distribution matching SPRINT
  pred_sbp <- 129.38647 + 0.13493565 * df$Age + 2.0699322 * df$Female + 0.7088606 * df$Black + 1.4437318 * df$Hispanic
  return (pred_sbp + 15.496434586232125 * rnorm(dim(df)[1],0,1)) 
}

sample_dbp <- function(df){
  # Sample subject baseline (seated) diastolic blood pressure using conditional normal distribution matching SPRINT
  pred_dbp <- 129.38647 + 0.13493565 * df$Age + 2.0699322 * df$Female + 0.7088606 * df$Black + 1.4437318 * df$Hispanic
  return (pred_dbp + 8.428068137949666 * rnorm(dim(df)[1],0,1)) 
}

sample_n_agents <- function(df){
  # Sample subject # of medications prescribed using conditional poisson distribution matching SPRINT RCT
  pred_lam <- exp(1.2004554953534512 + 9.21257164e-04 * df$Age + 3.70157118e-02 * df$Female + 
                    1.36445705e-01 * df$Black - 3.54435776e-02 * df$Hispanic -
                    3.06197274e-05 * df$SBP - 9.08632450e-03 * df$DBP)
  return (rpois(dim(df)[1],pred_lam)) 
}

sample_curr_smoker <- function(df){
  # Sample subject current smoker status (1 = px is current smkr) using conditional bernoulli dist matching SPRINT
  prob_curr_smoker <- 0.70086056 - 0.00966906 * df$Age - 0.00343624 * df$Female + 
                      0.08887907 * df$Black - 0.03297873 * df$Hispanic+ 
                      0.0009654 * df$SBP - 0.00051291 * df$DBP - 0.01574018 * df$N_Medications
  return (as.numeric(runif(dim(df)[1]) < prob_curr_smoker)) 
}

sample_former_smoker <- function(df){
  # Sample subject former smoker status (1 = px is former smkr) using conditional bernoulli dist matching SPRINT
  prob_former_smoker <- 0.47079387 + 2.3620378e-03 * df$Age - 1.5869808e-01 * df$Female - 
                        5.0397918e-02 * df$Black - 1.2431867e-01 * df$Hispanic + 
                        1.5233239e-04 * df$SBP - 1.3951637e-03 * df$DBP + 
                        1.5536858e-02 * df$N_Medications - 4.5159003e-01 * df$CurrentSmoker
  return (as.numeric(runif(dim(df)[1]) < prob_former_smoker)) 
}

sample_aspirin <- function(df){
  # Sample subject daily aspirin use (1 = daily aspirin user) using conditional bernoulli dist matching SPRINT RCT
  prob_aspirin <- 0.3764397 + 0.00573419 * df$Age - 0.08851127 * df$Female - 
                  0.11044942 * df$Black - 0.13806044 * df$Hispanic - 
                  0.00074539 * df$SBP - 0.0022618 * df$DBP + 
                  0.05358611 * df$N_Medications - 0.01058393 * df$CurrentSmoker + 
                  0.02378987 * df$FormerSmoker
  return (as.numeric(runif(dim(df)[1]) < prob_aspirin)) 
}

sample_statin <- function(df){
  # Sample subject 'on any statin' (1=yes) using conditional bernoulli distribution matching SPRINT RCT
  prob_statin <- 0.29434785 + 0.0047816 * df$Age - 0.06067061 * df$Female - 
                 0.07177708 * df$Black + 0.01000513 * df$Hispanic - 
                 0.00152325 * df$SBP - 0.00205615 * df$DBP + 
                 0.05148918 * df$N_Medications + 0.0506965 * df$CurrentSmoker + 
                 0.06907348 * df$FormerSmoker + 0.2043268 * df$Aspirin
  return (as.numeric(runif(dim(df)[1]) < prob_statin)) 
}

sample_creatinine <- function(df){
  # Sample subject serum creatinine (mg/dL) using conditional log-normal distribution matching SPRINT RCT
  pred_log_creatinine <- -0.22361551 + 4.8064794e-03 * df$Age - 2.2061218e-01 * df$Female + 
                         1.0354647e-01 * df$Black - 6.6799179e-02 * df$Hispanic - 
                         7.2427012e-04 * df$SBP + 1.3587828e-04 * df$DBP + 
                         4.2082604e-02 * df$N_Medications - 3.7121084e-02 * df$CurrentSmoker - 
                         1.4266469e-02 * df$FormerSmoker- 4.6004378e-03 * df$Aspirin + 
                         1.6370658e-02 * df$Statin
  return (exp(pred_log_creatinine + 0.28 * rnorm(dim(df)[1],0,1)))
}

sample_cholesterol <- function(df){
  # Sample subject cholesterol (mg/dL) using conditional log-normal distribution matching SPRINT RCT
  pred_log_cholesterol <- 5.2979894 - 0.00220362 * df$Age + 0.12154006 * df$Female - 
                          0.00485418 * df$Black - 0.00789208 * df$Hispanic + 
                          0.00014631 * df$SBP + 0.00139958 * df$DBP - 
                          0.01902662 * df$N_Medications - 0.01929336 * df$CurrentSmoker + 
                          0.00329189 * df$FormerSmoker - 0.02225017 * df$Aspirin - 
                          0.11290579 * df$Statin + 0.00222458 * df$Creatinine
  return (exp(pred_log_cholesterol + 0.19 * rnorm(dim(df)[1],0,1)))
}

sample_hdl <- function(df){
  # Sample subject HDL-cholesterol direct (mg/dL) using conditional log-normal distribution matching SPRINT RCT
  pred_log_hdl <- 3.1693678 + 0.0058122 * df$Age + 0.12116382 * df$Female + 
                  0.07040945 * df$Black - 0.0429091 * df$Hispanic + 
                  0.0007821 * df$SBP - 0.00051581 * df$DBP - 
                  0.01099174 * df$N_Medications - 0.00403379 * df$CurrentSmoker + 
                  0.00821252 * df$FormerSmoker + 0.00567395 * df$Aspirin + 
                  0.02488792 * df$Statin - 0.06565224 * df$Creatinine + 
                  0.0016482 * df$Cholesterol
  return (exp(pred_log_hdl + 0.235 * rnorm(dim(df)[1],0,1)))
}

sample_triglycerides <- function(df){
  # Sample subject Triglycerides (mg/dL) using conditional log-normal distribution matching SPRINT RCT
  pred_log_triglycerides <- 4.5964236 - 0.00241438 * df$Age + 0.08382218 * df$Female - 
                            0.17654029 * df$Black + 0.06732232 * df$Hispanic + 
                            0.00018662 * df$SBP + 0.00050756 * df$DBP + 
                            0.01851496 * df$N_Medications + 0.07510402 * df$CurrentSmoker + 
                            0.0260251 * df$FormerSmoker + 0.00043784 * df$Aspirin + 
                            0.11768849 * df$Statin + 0.08381014 * df$Creatinine + 
                            0.00591391 * df$Cholesterol - 0.02088042 * df$HDL_cholesterol
  return (exp(pred_log_triglycerides + 0.53 * rnorm(dim(df)[1],0,1)))
}

sample_bmi <- function(df){
  # Sample subject Body Mass Index (kg/m^2) using conditional log-normal distribution matching SPRINT RCT
  pred_log_bmi <- 3.7662797 - 3.8123669e-03 * df$Age + 4.1947711e-02 * df$Female + 
                  2.9479891e-02 * df$Black - 2.3104459e-02 * df$Hispanic - 
                  1.2005708e-03 * df$SBP + 2.5209847e-03 * df$DBP + 
                  2.2204543e-02 * df$N_Medications - 9.0576112e-02 * df$CurrentSmoker + 
                  7.7617536e-03 * df$FormerSmoker + 8.1281550e-03 * df$Aspirin + 
                  7.3386673e-03 * df$Statin - 2.4926960e-02 * df$Creatinine +
                  9.5494324e-06 * df$Cholesterol - 3.6268646e-03 * df$HDL_cholesterol - 
                  1.8088980e-05 * df$Triglycerides
  return (exp(pred_log_bmi + 0.171 * rnorm(dim(df)[1],0,1)))
}

get_SPRINT_input_vars_df <- function(n){
  # Generate subject covariates/input variables whose first and second moments match the SPRINT RCT.
  # 
  # Args:
  #   n: An integer representing the number of subjects whose data should be simulated
  # 
  # Returns:
  #   An [n x 17] pandas dataframe where each row is a subject and each column is an input variable.
  # NOT USING Diabetes variable
  df <- data.frame(matrix(NA,n,16))
  names(df) <- c("Age","Female","Black","Hispanic","SBP","DBP","N_Medications","CurrentSmoker","FormerSmoker",
                 "Aspirin","Statin","Creatinine","Cholesterol","HDL_cholesterol","Triglycerides","BMI")
  df$Age <- sample_age(n)
  df$Female <- sample_gender(df)
  df$Black <- sample_race(df)
  df$Hispanic <- sample_hispanic(df)
  df$SBP <- sample_sbp(df)
  df$DBP <- sample_dbp(df)
  df$N_Medications <- sample_n_agents(df)
  df$CurrentSmoker <- sample_curr_smoker(df)
  df$FormerSmoker <- sample_former_smoker(df)
  df$Aspirin <- sample_aspirin(df)
  df$Statin <- sample_statin(df)
  df$Creatinine <- sample_creatinine(df)
  df$Cholesterol <- sample_cholesterol(df)
  df$HDL_cholesterol <- sample_hdl(df)
  df$Triglycerides <- sample_triglycerides(df)
  df$BMI <- sample_bmi(df)
  #df$Diabetes <- rep(0,n)
  
  return (df)
}

SPRINT_main_outcomes_weibull_aft_cumulative_hazard <- function(df,
                                                               W,
                                                               e,
                                                               unique_times,
                                                               fixed_treatment_effect=NULL,
                                                               variable_treatment_effect=NULL,
                                                               event_rate_multiplier=1.0){
  # Simulate potential outcome cumulative hazard curves using a Weibull AFT model based on SPRINT RCT data
  # 
  # Given a set of subjects and their CVD risk-related covariates, estimate the cumulative hazard
  # Parameters for this model were estimated by fitting a Weibull AFT model, as implemented in lifelines 
  #(see https://lifelines.readthedocs.io/en/latest/fitters/regression/WeibullAFTFitter.html), on the actual SPRINT RCT
  # data, extracting the parameters associated with the estimated survival curve, and then plugging those back into
  # the appropriate place in the cumulative hazard function with modifications via the variable_treatment_effect f'n.
  # 
  #   Args:
  #       df: A [num_subjects x num_features] pandas dataframe with at least those columns in `get_SPRINT_input_vars_df`
  #       W: An array with shape [num_subjects] where W[i] = 1 if the i'th subject was treated and 0 otherwise
  #       e: Either an array with shape [n] or float between 0 and 1 representing treatment propensities for each subject
  #       unique_times: An [num_times] numpy array representing the times at which to evaluate the cumulative hazard.
  #       fixed_treatment_effect: A python function that takes no arguments and returns a scalar value indicating the
  #                               log of the factor by which. In practice, the output of `fixed_treatment_effect` is incorporated into
  #                               the cumulative hazard for each subject such that, if the `fixed_treatment_effect` always returned a
  #                               value $\texttt{fte}$, then his/her cumulative hazard would be multiplied by a factor of
  #                                $$\exp \left(\texttt{fte}\right)$$
  #                               at every time point in consideration. For `fixed_treatment_effect=-1`, this corresponds to a
  #                               uniform 37% reduction (log(0.37)=-1) in cumulative hazard at each time point for those who were treated. In theory,
  #                               in the case of the Weibull AFT model, the `fixed_treatment_effect` is redundant with the
  #                              `variable_treatment_effect` option in that one could let
  #                               $$\texttt{vte}(x) = \frac{-\texttt{fte}}{\rho \cdot (W - e)}$$
  #                               for all values of $x$ and this would achieve the same thing as assigning `variable_treatment_effect` to
  #                               `None` and `fixed_treatment_effect` to $\texttt{fte}$. We may eliminate the `fixed_treatment_effect` option
  #                               in the future for this very reason. In the meantime, though, the two offer somewhat distinct intuitions for
  #                               heterogeneous treatment effect functional forms.
  #       variable_treatment_effect: A python function that takes in a [num_subjects x num_features] pandas dataframe
  #                                  like `df` and returns a float vector corresponding to a transformed per-subject treatment effect
  #                                  conditioned on each subject's covariates. In practice, the output of `variable_treatment_effect` is
  #                                  incorporated into the cumulative hazard for each subject such that if `variable_treatment_effect(X[0,:])`
  #                                  returns $\texttt{vte}(x)$ for the 0th subject and that subject had a treatment propensity `e` and
  #                                  treatment assignment `W`, then his/her cumulative hazard would be multiplied by a factor of
  #                                  $$\exp \left(-\rho \cdot (W - e) \cdot \texttt{vte}(x) \right)$$
  #                                  at every time point in consideration. For `variable_treatment_effect=1` with `W=1` and `e=0.5`,
  #                                  this corresponds to a 55% reduction (log(0.55)=1*-\rho \cdot (W - e)) in cumulative hazard at each time point 
  #                                  for those who were treated. For `variable_treatment_effect=-1`, this corresponds to an 83% increase 
  #                                  (log(1.83)=(-1)*-\rho \cdot (W - e)) in cumulative hazard at each time point.
  #       event_rate_multiplier: A scalar indicating (at a high level) by how much to increase/decrease the observed
  #                              event rate in the simulated data. More concretely, when we include an `event_rate_multiplier` of
  #                              size $\texttt{erm}$, then our cumulative hazard at each time point increases by a factor of
  #                              $$\texttt{erm}^{\rho}$$
  #                              and $\rho$ is estimated from the observed RCT SPRINT data to be $e^{0.191244}$.
  #                              In practice, an `event_rate_multiplier` of 2.0 turns out to approximately double the number of
  #                              observed events over the course of the simulated study. This can be helpful if you want assess
  #                              SPRINT-like data under the scenario where the outcome is much more common/prevalent/frequent.
  # 
  #   Returns:
  #       A [num_subjects x num_unique_times] array, representing the cumulative hazard for each subject specified
  #       in the input data frame at each of the time points given, assuming that the functional form mapping subject
  #       input covariates to cumulative hazard at time t is approximately the same as what we observe/estimate from
  #       the control individuals in the original SPRINT RCT.
  
  
  log_lambda <- 15.1 + 
                -0.053699 * df$Age + 
                 0.260077 * df$Female + 
                 0.010973 * df$Black + 
                 0.264252 * df$Hispanic + 
                -0.003996 * df$SBP + 
                -0.002673 * df$DBP + 
                -0.148897 * df$N_Medications + 
                -0.651439 * df$CurrentSmoker + 
                -0.143992 * df$FormerSmoker + 
                -0.069259 * df$Aspirin + 
                -0.065738 * df$Statin + 
                -0.326439 * df$Creatinine + 
                -0.003859 * df$Cholesterol + 
                 0.010695 * df$HDL_cholesterol + 
                -0.000297 * df$Triglycerides + 
                -0.002829 * df$BMI
  
  if(is.null(variable_treatment_effect)==FALSE){
    log_lambda <- log_lambda + (W - e) * variable_treatment_effect(df) # If e -> 1, then tau -> 0 for treated individuals
                                                                       # reprogram the "variable_treatment_effect" python function! 
  }else{
    log_lambda <- log_lambda - log(event_rate_multiplier)  # If event_rate_multiplier is 2, everyone's timeline is accelerated 2x
  }
  rho <- exp(0.180197)
  log_cum_hazard <- rho*(t(matrix(rep(log(unique_times),length(log_lambda)), length(unique_times), length(log_lambda))) -log_lambda) 
  
  if(is.null(fixed_treatment_effect)==FALSE){
    log_cum_hazard <- log_cum_hazard + fixed_treatment_effect # reprogram this function!
  }
  
  cum_hazard <- exp(log_cum_hazard)
  return(cum_hazard)
}

SPRINT_main_outcomes_weibull_aft_survival <- function(df,
                                                      W,
                                                      e,
                                                      unique_times,
                                                      treated=FALSE,
                                                      fixed_treatment_effect=NULL,
                                                      variable_treatment_effect=NULL,
                                                      event_rate_multiplier=1.0){
  # Estimate potential outcomes (as survival curves, estimated at given times) under treatment or control
  # 
  # This function takes in a set of subject covariates, (true) treatment assignments, treatment propensities, and
  # arguments regarding the form of heterogeneity in treatment effect, and returns the potential outcome
  # (in the form of a survival curve estimated at the time points specified by the user) under either treatment or
  # control for each subject in the dataset.
  #
  # Args:
  #   treated: A boolean True/False indicating whether the potential outcomes should be calculated considering the
  #            outcomes under treatment (i.e., `True`) or under control (i.e., `False`)
  # 
  # Returns:
  #   A [num_subjects x num_unique_times] array where the entry in row i and column j represents the probability
  #   that the i^th subject's failure time exceeds unique_times[j] (assuming that subject is treated if treated==True
  #   and that the subject is not treated if treated==False)
  
  if(treated==TRUE){
    return (exp(-SPRINT_main_outcomes_weibull_aft_cumulative_hazard(
                 df,
                 W=rep(1,dim(df)[1]),  # Everyone is treated
                 e=e,
                 unique_times=unique_times,
                 fixed_treatment_effect=fixed_treatment_effect,
                 variable_treatment_effect=variable_treatment_effect,
                 event_rate_multiplier=event_rate_multiplier)))
  }else{
    return (exp(-SPRINT_main_outcomes_weibull_aft_cumulative_hazard(
                 df,
                 W=rep(0,dim(df)[1]),  # Everyone is treated
                 e=e,
                 unique_times=unique_times,
                 fixed_treatment_effect=NULL,
                 variable_treatment_effect=NULL,
                 event_rate_multiplier=event_rate_multiplier)))
  }
}

SPRINT_censoring_weibull_cumulative_hazard <- function(unique_times){
  # Simulate a cumulative hazard curve that looks like SPRINT RCT, estimated at a number of specified time points
  cum_hazard = exp(5.121360 * (log(unique_times) - log(1298.354641)))
  return(cum_hazard)
}

SPRINT_censoring_weibull_survival <- function(df, unique_times){
  # Simulate a survival curve for censoring times that looks like what we see in the SPRINT RCT, using cum. hazard
  # 
  # We assume that the censoring survival curve is the same for all the subjects (i.e., that censoring is not dependent
  # on the covariates). This may be overly simplistic in reality but, in our experience, modeling covariate-dependent
  # censoring increases the variance of downstream estimates (e.g., the CATE) while offering only limited benefits
  # in terms of decreased bias. So we invoke Occam's Razor for the purpose of our simulator.
  # 
  #   Args:
  #       df: A [num_subjects x num_features] pandas dataframe with columns specified as per `get_SPRINT_input_vars_df`
  #       unique_times: unique_times: A [num_unique_times] array indicating the times at which the corresponding
  #           (censoring) survival probabilities should be estimated, such that survival_curve[i, j] represents the
  #           probability that the censoring time for subject i is greater than unique_times[j]).
  # 
  #   Returns:
  #       A [num_subjects x num_unique] times array, sc, where sc[i, j] represents the probability that the censoring time
  #       for subject i is/was greater than unique_times[j].
  
  survival_curve <- exp(-SPRINT_censoring_weibull_cumulative_hazard(unique_times))
  sc <- t(matrix(rep(survival_curve,dim(df)[1]),length(survival_curve),dim(df)[1]))
  return(sc)
}

sample_T_from_sc <- function(sc,unique_times){
  # Sample time-to-event for a number of subjects given their survival curve using inverse transform sampling
  # 
  # We use inverse transform sampling to sample failure times for each subject given their survival curves. Concretely,
  # this is accomplished by (1) drawing a uniform random variable, U; (2) Finding
  # 
  # Args:
  #   sc: A [num_subjects x num_unique_times] data frame representing the survival curves for each subject, such that
  #       sc[i, j] represents the probability that the failure time for subject i is larger than time unique_times[j].
  #   unique_times: A [num_unique_times] array indicating the times corresponding to the survival probabilities
  #                 (e.g., such that sc[i, j] represents the prob. that failure time for subject i is > time unique_times[j]).
  # 
  # Returns:
  #   An array with shape [num_subjects] representing the sampled failure times for each subject.
  # To use:
  #   >>> # 55yo female, black, with SBP 120 mm Hg, current smoker, HDL-Cholesterol 80 mm Hg, Total cholesterol 150 mm Hg
  #   >>> X = np.array([[55, 1, 1, 1, 120, 80, 5, 1, 1, 1, 1, 150, 150, 80, 80, 25, 0]])
  #   >>> calculate_ascvd_risk(X)
  #   >>> sc_np = np.array(
  #     ...     [[0.90, 0.88, 0.85, 0.81],  # Row 0 is Survival curve for subject 0
  #              ...     [0.70, 0.68, 0.65, 0.61],  # Row 1 is Survival curve for subject 1
  #              ...     [0.80, 0.78, 0.75, 0.71],
  #              ...     [0.60, 0.58, 0.55, 0.51]])
  #   >>> unq_times = np.array([1, 2, 3, 4])
  #   >>> sim_results_np = np.array([sample_T_from_sc(sc_np, unq_times) for _ in range(100000)])
  #   >>> np.mean(sim_results_np > 2, axis=0)  # Find the empirical probs w/ which time-to-event is greater than 2
  #   array([0.86336, 0.66507, 0.76552, 0.5654 ])
  #   
  #   TODO: Empirical probabilities of \hat{Pr}(Y > T) = 1/2 * (Pr(Y > T) + Pr(Y > T + 1)) for survival time Y and
  #   TODO: query time T (in other words, the empirical probabilities are not quite aligned with what we'd expect,
  #   TODO: but the difference is of immaterial impact given a reasonable granularity.
  
  # Inverse transform sampling using the CDF induced by the survival curve
  U_vec <- runif(dim(sc)[1])        # Generate as many random uniform r.v.s as there are subjects
  dist_to_U <- abs(sc - U_vec)      # For each subject, find survival prob closest to sampled U and...
  dist_to_U_min <- apply(dist_to_U, 1, FUN=min)
  index <- rep(NA, dim(sc)[1])
  for (i in 1:dim(sc)[1]){
    index[i] <- which(dist_to_U[i,] == dist_to_U_min[i])
  }
  time_to_event <- unique_times[index]     # ...extract time associated with that survival prob
  return(time_to_event)
}

# TODO: Make this function more (time and memory) efficient
simulate_SPRINT_outcomes <- function (df,
                                      e,
                                      unique_times,
                                      treated=FALSE,
                                      fixed_treatment_effect=NULL,
                                      variable_treatment_effect=NULL,
                                      event_rate_multiplier=1.0){
  # Generate potential outcome survival curves and failure times for both under specified treatments
  # 
  # Args:
  #   df: A [num_subjects x num_features] pandas dataframe with columns specified as per `get_SPRINT_input_vars_df`
  #   e: Either an array with shape [n] or float between 0 and 1 representing treatment propensities for each subject
  #   unique_times: An array indicating the unique times at which the survival curve probabilities should be estimated
  #   treated: Whether the survival curves should be generated under treatment vs. control potential outcomes
  #   fixed_treatment_effect: See documentation for SPRINT_main_outcomes_weibull_aft_cumulative_hazard
  #   variable_treatment_effect: See documentation for SPRINT_main_outcomes_weibull_aft_cumulative_hazard
  #   event_rate_multiplier: See documentation for SPRINT_main_outcomes_weibull_aft_cumulative_hazard
  # 
  # Returns:
  #   A tuple with four elements:
  #       po_times: An array with shape [num_unique_times] indicating the sampled times at which a failure
  #                 time/primary outcome was "observed" in the simulated data.
  #       censoring_times: An array with shape [num_unique_times] indicating the sampled times at which the subject's
  #                        censoring was observed. Note that this is sampled independently from the primary outcome/failure time.
  #       po_sc: A [num_subjects x num_unique_times] data frame representing the (failure time/primary outcome)
  #              survival curves for each subject, such that sc[i, j] represents the probability that the failure time
  #              for subject i is larger than time unique_times[j].
  #       censoring_sc: A  [num_subjects x num_unique_times] data frame representing the (censoring)
  #                     survival curves for each subject, such that sc[i, j] represents the probability that the censor time
  #                     for subject i is larger than time unique_times[j].
  
  # Generate survival curves for primary outcomes estimated at unique times given
  po_sc <- SPRINT_main_outcomes_weibull_aft_survival(df,
                                                     e,
                                                     unique_times=unique_times,
                                                     treated=treated,
                                                     fixed_treatment_effect=fixed_treatment_effect,
                                                     variable_treatment_effect=variable_treatment_effect,
                                                     event_rate_multiplier=event_rate_multiplier)
  po_times <- sample_T_from_sc(po_sc, unique_times=unique_times)  # Sample failure times from the survival curve
  
  # Generate survival curves for censoring estimated at the same unique times as above
  censoring_sc <- SPRINT_censoring_weibull_survival(df,
                                                    unique_times=unique_times)
  censoring_times <- sample_T_from_sc(censoring_sc, unique_times=unique_times)  # Sample censoring times from the sc
  
  outcomes <- list(NA,4)
  outcomes$po_sc <- po_sc
  outcomes$po_times <- po_times
  outcomes$censoring_sc <- censoring_sc
  outcomes$censoring_times <- censoring_times
  return(outcomes)
}
  
estimate_treatment_effect_via_monte_carlo <- function(sc_treat, sc_control, unique_times, end_time, n_mc=1000){
  # Estimate the CATE for each subject using Monte Carlo sampling
  # 
  # Args:
  #   sc_treat: A [num_subjects x num_unique_times] survival curve for subjects under treatment
  #   sc_control: A [num_subjects x num_unique_times] survival curve for subjects under control
  #   unique_times: An array with shape [num_unique_times] indicating the times corresponding to columns
  #                 of the survival curves sc_treat and sc_control.
  #   end_time: The effective end time of the study. See `simulate_parametric_SPRINT_data` for details.
  #   n_mc: The number of monte carlo draws to estimate the treatment effect with. Default is 1000.
  # 
  # Returns:
  #   An array with shape [num_subjects] representing the estimated CATE for each subject
  
  ft1 <- rep(NA,n_mc)
  ft0 <- rep(NA,n_mc)
  for (i in 1:n_mc){
    ft1[i] <- sample_T_from_sc(sc_treat, unique_times)
    ft0[i] <- sample_T_from_sc(sc_control, unique_times)
  }
  
  cate <- mean(pmin(ft1,end_time) - pmin(ft0,end_time))  
}

simulate_parametric_SPRINT_data <- function(n,
                                            e,
                                            fixed_treatment_effect_strength=0.0,
                                            variable_treatment_effect_type='ascvd_correlated',
                                            variable_treatment_effect_strength=0.0,
                                            event_rate_multiplier=1.0,
                                            end_time=365.25*3,
                                            base_to_round_to=1,
                                            calculate_mc_cate=FALSE,
                                            n_mc=1000){
  # Simulate input and output/survival outcome data that closely represents the SPRINT Randomized Controlled Trial
  # 
  # Args:
  #   n: Integer representing the total number of subjects whose data we should simulate
  #   e: Either an [n x 1] array-like or float between 0 and 1 representing treatment propensities for each subject
  #   fixed_treatment_effect_strength: as above
  #   variable_treatment_effect_type:
  #   variable_treatment_effect_strength: as above
  #   end_time: The effective end time of the study. This is used for calculating the Restricted Mean Survival Time.
  #             An important assumption for identifying the conditional average treatment effect tau(X) is that there
  #             exists a fixed positive constant M such that the probability of observing an event time past the
  #             maximum follow-up time end_time is at least M (formally, we assume: P(Y >= end_time | X) > M).
  #   base_to_round_to: An integer to which all survival times should be rounded, such that e.g. survival times of
  #                     [102.49, 103.19, 106.13, 107.95] w/ base_to_round_to=5 yields [100, 105, 105, 110]
  #   calculate_mc_cate: Whether to calculate the CATE for each subject using Monte Carlo Sampling (estimating the
  #                      CATE using the RMST derived from the survival curves is performed by default/automatically.
  #   n_mc: The number of monte carlo draws to estimate the treatment effect with. Default is 1000.
  #
  # Returns:
  #   A dictionary with the following keys and associated values:
  #     'X': A [num_unique_subjects x num_features] pandas dataframe where the columns are named in
  #          accordance with what they represent in the original SPRINT RCT.
  #     'Y': An array with shape [num_unique_subjects] representing the observed study time, which is the minimum
  #          of the failure/primary outcome time and censoring time.
  #     'e': An array with shape [num_unique_subjects] representing the probability with which each subject was
  #          assigned to treatment.
  #     'W': An array with shape [num_unique_subjects] s.t. W[i] = 1 if subject i was treated and 0 otherwise
  #     'D': An array with shape [num_unique_subjects] representing the event type such that D[i] = 0 if the
  #          i^th subject was censored and 1 if the i^th had an observed failure/primary outcome
  #     'tau': An array with shape [num_unique_subjects] where tau[i] represents the CATE associated with subject i,
  #            calculated as restricted mean survival time (RMST) under treatment minus RMST under control/no treatment
  #     'tau_sc_diff': An array with shape [num_unique_subjects] where tau_sc_diff[i] represents the survival
  #                    probability under treatment at study end time minus the survival probability under control at study end.
  #     'tau_mc': An array with shape [num_unique_subjects] where tau_mc[i] represents the (Monte Carlo sampled)
  #               difference in time-to-failure under treatment vs. control. Should approach tau asymptotically.
  #     'failure_time_treat':
  #     'failure_time_ctrl':
  #     'failure_time':
  #     'censor_time':
  #     'po_sc_treated': A [num_unique_subjects x num_unique_times] array representing the simulated survival curve
  #                      for each subject under treatment.
  #     'po_sc_control': A [num_unique_subjects x num_unique_times] array representing the simulated survival curve
  #                      for each subject under control.
  #     'end_time': The study end time, used in calculating the Restricted Mean Survival Time.
  
  # TODO: Bake in a decent way of estimating default value (see Steve's code)
  rounded_end_time <- round_to_nearest_base(end_time, base_to_round_to)
  
  # These times are used for the purpose of simulation via the `sample_T_from_sc` function.
  survival_curve_times <- seq(base_to_round_to, rounded_end_time * 2, base_to_round_to)
  
  # A little kludgy right now (needs work). We have the user pass in a treatment effect strength, and then
  # convert that into a function for them. This enables us to pass simple scalar arguments to the program at
  # the command line, but it also means that user's can't pass in a function directly if they so desire.
  # If no "fixed treatment effect" is specified, then we assume that
  fixed_treatment_effect_fn <- fixed_treatment_effect_strength
  
  variable_treatment_effect_fn = NULL
  if(variable_treatment_effect_strength != 0){
    if(variable_treatment_effect_type == "ascvd_correlated"){
      variable_treatment_effect_fn <- function(df){
        variable_treatment_effect_strength * calculate_ascvd_risk(df)
      }}else if(variable_treatment_effect_type == "ascvd_anticorrelated"){
        variable_treatment_effect_fn <- function(df){
          -variable_treatment_effect_strength * calculate_ascvd_risk(df)
        }
      }
  }
  
  X <- get_SPRINT_input_vars_df(n)  # Get study subject covariates
  if(is.integer(e)==FALSE){
    e <- rep(e,dim(X)[1])
  }
  W <- as.numeric(runif(n) <= e)    # Assign each subject to treatment/control based on propensities
  
  # Generate potential outcomes under treatment, where everyone has W=1
  po_treated <- simulate_SPRINT_outcomes(X,
                                         e,
                                         unique_times=survival_curve_times,
                                         treated=TRUE,
                                         fixed_treatment_effect=fixed_treatment_effect_fn,
                                         variable_treatment_effect=variable_treatment_effect_fn,
                                         event_rate_multiplier=event_rate_multiplier)
  
  failure_time_treat <- po_treated$po_times
  censor_time_treat<- po_treated$censoring_times
  po_sc_treat <- po_treated$po_sc
  censor_sc <- po_treated$censoring_sc

  # Generate potential outcomes under treatment, where everyone has W=0
  po_control <- simulate_SPRINT_outcomes(X,
                                         e,
                                         unique_times=survival_curve_times,
                                         treated=FALSE,
                                         fixed_treatment_effect=NULL,
                                         variable_treatment_effect=NULL,
                                         event_rate_multiplier=event_rate_multiplier)
  
  failure_time_ctrl <- po_control$po_times
  censor_time_ctrl <- po_control$censoring_times
  po_sc_control <- po_control$po_sc
  censor_sc <- po_control$censoring_sc
  
  failure_time <- W * failure_time_treat + (1-W) * failure_time_ctrl
  #failure_time <- pmin(failure_time, rounded_end_time)        # Truncate failure times to the end time
  censor_time <- W * censor_time_treat + (1-W) * censor_time_ctrl
  
  is_censored <- as.numeric(censor_time < failure_time)
  D <- as.numeric(failure_time <= censor_time)
  Y <- pmin(failure_time, censor_time)
  
  if(is.null(base_to_round_to)==FALSE){
    failure_time <- round_to_nearest_base(failure_time, base_to_round_to)
    censor_time <- round_to_nearest_base(censor_time, base_to_round_to)
    Y <- round_to_nearest_base(Y, base_to_round_to)
  }
  
  # Estimate "True" treatment effect using RMST (aligns with units of AIPW scores in Yifan's paper)
  cate_rmst <- estimate_treatment_effect_in_rmst(po_sc_control, po_sc_treat, survival_curve_times, rounded_end_time)  
  
  # Estimate "True" treatment effect using monte carlo sampling
  if(calculate_mc_cate==TRUE){
    cate_mc <- estimate_treatment_effect_via_monte_carlo(po_sc_treat,
                                                         po_sc_control,
                                                         survival_curve_times,
                                                         rounded_end_time,
                                                         n_mc)
  }
  
  # Estimate "True" treatment effect by taking difference between survival prob at study end time under treatment
  # and survival prob at study end time under control. (NOT in the same units as RMST, for example)
  end_time_indx <- which(survival_curve_times == rounded_end_time)
  cate_sc_diff <- po_sc_treat[, end_time_indx] - po_sc_control[, end_time_indx]
  
  results_dict <- list(NA,16)
  results_dict$X <- X
  results_dict$Y <- Y
  results_dict$e <- e
  results_dict$W <- W
  results_dict$D <- D
  results_dict$tau <- cate_rmst
  results_dict$tau_sc_diff <- cate_sc_diff
  if(calculate_mc_cate==TRUE){
    results_dict$tau_mc <- cate_mc
  }else{
    results_dict$tau_mc <- NULL
  }
  results_dict$failure_time_treat <- failure_time_treat
  results_dict$failure_time_ctrl <- failure_time_ctrl
  results_dict$failure_time <- failure_time
  results_dict$censor_time <- censor_time
  results_dict$po_sc_treated <- po_sc_treat  # Should have shape (num_samples, num_unique_times)
  results_dict$po_sc_control <- po_sc_control
  results_dict$survival_curve_times <- survival_curve_times
  results_dict$end_time <- rounded_end_time
  
  return(results_dict)
  parser <- ArgumentParser()
  parser$add_argument('--num_samples',
                      type="integer",
                      default=100,
                      help='Number of samples (default 100)')
    
  parser$add_argument('--prob_treat',
                      type="double",
                      default=0.5023,
                      help='Probability of each individual receiving treatment (assumed RCT, default 0.5023 [from SPRINT])')
  parser$add_argument('--out_path',
                      type="character",
                      default='tmp_SPRINT_parametric_out.csv',
                      help='Path indicating where to save the simulated SPRINT data')
    
  parser$add_argument('--fte',
                      type="double",
                      default=0.0,
                      help='Strength of fixed treatment effect (i.e., constant reduction in cum. hazard) for survival outcomes')
    
  parser$add_argument('--vte_type',
                      type="character",
                      default='ascvd_correlated',
                      choices=c('ascvd_correlated', 'ascvd_anticorrelated'),
                      help='Type of variable treatment effect (i.e., variable effect on cum. hazard in AFT model) for survival data')
    
  parser$add_argument('--vte_strength',
                      type="double",
                      default=0.0,
                      help='Strength of variable treatment effect for survival outcomes')
    
  parser$add_argument('--round_times_to_int',
                      type="integer",
                      default=1,
                      help='Whether survival times should be rounded to nearest integer (default True, for memory efficiency)')
    
  parser$add_argument('--verbose',
                      type="integer",
                      default=0,
                      help='Whether to print out a snapshot of the data upon completion')
  
  args <- parser$parse_args()
  tick <- Sys.time()
  Out <- simulate_parametric_SPRINT_data(n=args$num_samples,
                                         e=args$prob_treat,
                                         fixed_treatment_effect=args$fte,
                                         variable_treatment_effect_type=args$vte_type,
                                         variable_treatment_effect_strength=args$vte_strength,
                                         base_to_round_to=args$round_times_to_int)
  
  simulated_df <- data.frame(Out$X, Out$Y, Out$e, Out$W, Out$D, Out$tau, Out$tau_sc_diff, Out$tau_mc, 
                             Out$failure_time_treat,Out$failure_time_ctrl, Out$failure_time, 
                             Out$censor_time, Out$po_sc_treated, Out$po_sc_control, Out$survival_curve_times)
  
  tock <- Sys.time()
  print(paste0("Took ",round((tock-tick)/60,4)," minutes to simulate ",args$num_samples," samples:"))
  print(simulated_df)
  write.csv(args$out_path,row.names = FALSE)  
  print(paste0("Simulated data written to ",args$out_path))
}







