Matlab codes for replicating both the simulation and empirical results of Chen and Lee (2023) on the L0-regularized high dimensional quantile regression methods. Description of the L0-regularized quantile regression estimators and details of its numerical studies can be found in the paper:

Chen, Le-Yu and Lee, Sokbae (2023), "Sparse Quantile Regression".

The paper has been accepted for publication at Journal of Econometrics. The latest working paper version of this work can be found in this repository.

To replicate the simulation and empirical results of the paper, put all program and data files in the same work directory. The Matlab version of the Gurobi solver has to be installed for implementing the various quantile regression estimators studied in this paper. The Gurobi solver is freely available for academic purposes. You can download it via:

https://www.gurobi.com/

To replicate Tables 5.1, 5.2 and 5.3 of the paper, run the Matlab program "Simulation_Results_replication.m" and change the simulation settings therein.
To replicate Empirical results of the paper, run the program "Empirical_Results_replication.m" and change the covariate specification therein accordingly.
