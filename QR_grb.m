% Code for replicating numerical results reported in the paper:
% Chen Le-Yu and Lee Sokbae (2023), Sparse Quantile Regression, forthcoming at the Journal of
% Econometrics

% This function performs the standard quantile regression. 
% function input :
% y : vector of dependent variables
% x : (n by p) matrix of covariate data 
% tau : the quantile regression index
% T     : the time limit specified for early termination of the MIO solver
% bnd   : (p by 2) matrix where the first and second columns  
%         respectively store the lower and upper bounds 
%         of the unknown coefficients

% function output :
% bhat  : the QR estimates for the unknown coefficients
% obj_v : the value of QR objective function
% gap   : the MIO optimization gap value in case of early termination
%         gap = 0 ==> optimal solution is found within the time limit
% rtime : the time used by the MIO solver in the estimation procedure
% ncount : the number of branch-and-bound nodes used by the MIO solver

function [bhat,obj_v,gap,rtime,ncount] = QR_grb(y,x,tau,T,abgap,bnd)

[n,p]=size(x);

bhat=zeros(p,1);

gap=0;
rtime=0;
ncount=0;

model.modelsense = 'min';

model.sense = '=';
model.lb = [bnd(:,1);zeros(2*n,1)];
model.ub = [bnd(:,2);1/eps*ones(2*n,1)];

tol=1e-6;

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

model.obj = [zeros(p,1);tau*ones(n,1);(1-tau)*ones(n,1)];
model.A = sparse([x eye(n) -eye(n)]);
model.rhs = y;

try
    result = gurobi(model, params);
    bhat=result.x(1:p);
    obj_v=result.objval;
    gap=(obj_v-result.objbound);
    rtime=result.runtime;
    ncount=result.nodecount;
    
    %fprintf('Optimization returned status: %s\n', result.status);
  
catch gurobiError
    fprintf('Error reported\n');
end

end

