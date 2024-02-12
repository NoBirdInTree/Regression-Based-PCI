# Replication Files for Regression-Based Proximal Causal Inference (Liu, Park, Li, Tchetgen Tchetgen, _ARXIV_, 2024+) 

This Github repository contains replication files for Section A.8 Stimulation Study of the Appendix of the paper URegression-Based Proximal Causal Inference (Liu, Park, Li, Tchetgen Tchetgen, 2024)](https://arxiv.org/abs/2402.00335).

## Code

The repository includes the exact data-generation procedure for the synthetic dataset and the model fitting, as well as the comparison of estimates from different estimators mentioned in A.8. The following two files implement these procedures:

* AZU.R: This file generates the i.i.d. vectors (A_i,Z_i,U_i) following the joint probability functions derived in A.8 via the accept-reject sampling algorithm.

* BT_Evaluation.R: This file further generates the i.i.d vectors (Y_i,W_i,A_i,Z_i,U_i) from vectors (A_i,Z_i,U_i) specified above following the joint probability functions derived in A.8 via the accept-reject sampling algorithm. Then, using samples drawn from a large dataset composed of i.i.d. vectors (Y_i,W_i,A_i,Z_i,U_i), finite sample analysis is conducted. Three estimations -- Two-Stage, Naive, Orale -- are evaluated, and their performance at different sample sizes are reported.
