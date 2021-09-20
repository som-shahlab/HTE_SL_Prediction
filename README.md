# HTE_SL_Prediction
We developed the codes here to estimate heterogeneous treatment effects using statistical learning approaches 

This repository contains the following folders: 
- Data: 
-- SPRINT data 
-- ACCORD data 
- Models: Deriving prediction models for risk and HTE estimation
           1) We use a specturm of statistical learning methods including GBM, RSF, BART, Deepsurv
           2) We conduct internal validation using SPRINT data for both the CVD and SAE outcomes
           3) We also conduct internal validation using ACCORD data for the CVD outcome only 
           4) We implement external validation on SPRINT and ACCCORD, with the ACCORD and SPRINT as the training sets, respectively
- Xlearner: Implement the X-learner meta-algorithm using the counterfactual estimates from T-learners 
- Evaluations: Different evaluation metrics are used to assess the risk and HTE model performance 
           1) C-index 
           2) GND test 
           3) RATE metric 
           4) HTE calibration plot
- Summary: codes to generate summarized tables and figures for publication 
- Simulations: 1) Simulations based on the settings in the causal survival forest paper 
