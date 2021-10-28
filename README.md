# HTE_SL_Prediction
We developed the codes here to estimate heterogeneous treatment effects using statistical learning approaches 

This repository contains the following folders: 
- Data: 
  - SPRINT data 
  - ACCORD data 
- Models: Deriving prediction models for risk and HTE estimation
  - We use a specturm of statistical learning methods including GBM, RSF, BART, Deepsurv
  - We conduct internal validation using SPRINT data for both the CVD and SAE outcomes
  - We also conduct internal validation using ACCORD data for the CVD outcome only 
  - We implement external validation on SPRINT and ACCCORD, with the ACCORD and SPRINT as the training sets, respectively
- Xlearner: Implement the X-learner meta-algorithm using the counterfactual estimates from T-learners 
- Evaluations: Different evaluation metrics are used to assess the risk and HTE model performance 
   - C-index 
   - GND test 
   - RATE metric 
   - HTE calibration plot
- Summary: codes to generate summarized tables and figures for publication 
- Simulations: 
  - Simulations based on the settings in the causal survival forest paper 
