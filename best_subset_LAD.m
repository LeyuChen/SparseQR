
% Code for replicating numerical results reported in the paper:
% Chen Le-Yu and Lee Sokbae (2023), Sparse Quantile Regression, forthcoming at the Journal of
% Econometrics

% This function performs median regression subject to a constraint on the
% maximal number of selected variables. 
% We thank Rahul Mazumder for providing us the code for implementing this function. 

% function input :
% y : the outcome vector
% x : (n by p) matrix of the covariate data 
% q     : the cardinality constraint for the covariate selection 
% T     : the time limit specified for early termination of the MIO solver
% abgap : the absolute gap specified for early termination of the MIO solver
% bnd   : (p by 2) matrix where the first and second columns  
%         respectively store the lower and upper bounds 
%         of the unknown coefficients

% function output :
% bhat  : the estimated coefficients of the L0-constrained median regression
% obj_v : the value of the MIO objective function
% gap   : the MIO optimization gap value in case of early termination
%         gap = 0 ==> optimal solution is found within the time limit
% rtime : the CPU time used by the MIO solver 
% ncount: the number of branch-and-bound nodes used by the MIO solver 


function [bhat,obj_v,gap,rtime,ncount] = best_subset_LAD(y,x,q,T,abgap,bnd)


% run First order method for initialization
bb0 = zeros(size(x,2),1);
maxiter=1000; RUNS=10; tol=1e-6;
[~, ~, bet_best]  = FUN_bs_LP_FO(x,y,q,bb0,maxiter,RUNS,tol); 


%specify starting values and also BigM values:
sideInfo.beta_start=bet_best;
sideInfo.Linfty_bound=bnd;

%specfiy Gurobi params [e.g., TimeLimit]
params.outputflag = 0; 
params.OptimalityTol=tol;
params.FeasibilityTol=tol;
params.IntFeasTol=tol;

if T > 0
params.TimeLimit =T;
end

if abgap > 0
params.MIPGapAbs=abgap;
end

[results2, bhat, ~ ]  = FUNbs_LP_demo(x,y,q,sideInfo,params);

obj_v=results2.objval;
gap=(obj_v-results2.objbound);
rtime=results2.runtime;
ncount=results2.nodecount;
end

