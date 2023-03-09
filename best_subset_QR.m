
% Code for replicating numerical results reported in the paper:
% Chen Le-Yu and Lee Sokbae (2023), Sparse Quantile Regression, forthcoming at the Journal of
% Econometrics

% This function performs quantile regression subject to a constraint on the
% maximal number of selected variables. 

% function input :
% y : the outcome vector
% x : (n by p) matrix of the covariate data 
% q     : the cardinality constraint for the covariate selection 
% T     : the time limit specified for early termination of the MIO solver
% tau : the targeted quantile level 
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

function [bhat,obj_v,gap,rtime,ncount] = best_subset_QR(y,x,q,T,tau,abgap,bnd)

 
% run First order method for initialization

[bhat,~]  = best_subset_QR_FO(y,x,tau,q,bnd); 

tol=1e-6;

uhat = y-x*bhat;
bp = uhat.*(uhat>=0); 
bn = -uhat.*(uhat<0);
zvec= (abs(bhat)>tol);

[n,p]=size(x);

gap=0;
rtime=0;
ncount=0;

model.modelsense = 'min';
model.start=[bhat;bp;bn;zvec];      
model.sense = [repmat('=',1,n) repmat('<',1,2*p+1)];
model.vtype = char([67*ones(1,p+2*n) 66*ones(1,p)]); 
model.lb = [bnd(:,1);zeros(2*n+p,1)];
model.ub = [bnd(:,2);1/eps*ones(2*n,1);ones(p,1)];

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

QR_obj_v = [tau*ones(n,1);(1-tau)*ones(n,1)]/n;
temp1 = [x eye(n) -eye(n) zeros(n,p)];
temp2 = [-eye(p) zeros(p,2*n) diag(bnd(:,1))];
temp3 = [eye(p) zeros(p,2*n) -diag(bnd(:,2))];
temp4 = [zeros(1,p+2*n) ones(1,p)];
model.obj = [zeros(p,1);QR_obj_v;zeros(p,1)];
model.A = sparse([temp1;temp2;temp3;temp4]);
model.rhs = [y;zeros(2*p,1);q];

try
    result = gurobi(model, params);
    bhat=result.x(1:p);
    obj_v=result.objval;
    gap=(obj_v-result.objbound);
    rtime=result.runtime;
    ncount=result.nodecount;
    
    fprintf('Optimization returned status: %s\n', result.status);
  
catch gurobiError
    fprintf('Error reported\n');
end

end

