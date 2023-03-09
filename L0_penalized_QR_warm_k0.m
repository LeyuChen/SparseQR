
% Code for replicating numerical results reported in the paper:
% Chen Le-Yu and Lee Sokbae (2023), Sparse Quantile Regression, forthcoming at the Journal of
% Econometrics

% This function performs the MIO based L0-penalized quantile regression. 

% function input :
% y : the outcome vector
% x : (n by p) matrix of the covariate data 
% tau : the targeted quantile level 
% lamda : the penalty parameter given by (5.2) of the paper      
% T     : the time limit specified for early termination of the MIO solver
% abgap : the absolute gap specified for early termination of the MIO solver
% bnd   : (p by 2) matrix where the first and second columns  
%         respectively store the lower and upper bounds 
%         of the unknown coefficients
% tuning : the tuning constant c stated in (5.2) of Chen and Lee (2023)
% K0     : maximal number of variables allowed to be selected 

% function output :
% bhat  : the estimated coefficients of the L0-penalized quantile regression
% obj_v : the value of the MIO objective function
% gap   : the MIO optimization gap value in case of early termination
%         gap = 0 ==> optimal solution is found within the time limit
% rtime : the CPU time used by the MIO solver 
% ncount: the number of branch-and-bound nodes used by the MIO solver 
% warm_start_time : the CPU time used to warm start the L0-PQR procedure


function [bhat,obj_v,gap,rtime,ncount,warm_start_time] = L0_penalized_QR_warm_k0(y,x,tau,lamda,T,abgap,bnd,tuning,k0)

[n,p]=size(x);
timerVal = tic;
[bhat,obj_v] = L0_PQR_FO_k0(y,x,tau,lamda,bnd,tuning,k0);
warm_start_time=toc(timerVal);

tol=1e-6;
bhat=bhat';
uhat = y-x*bhat;
bp = uhat.*(uhat>=0); 
bn = -uhat.*(uhat<0);
zvec= (abs(bhat)>tol);
      
gap=0;
rtime=0;
ncount=0;

model.modelsense = 'min';
model.start=[bhat;bp;bn;zvec];      

model.sense = [repmat('=',1,n) repmat('<',1,2*p+1)];
model.lb = [bnd(:,1);zeros(2*n+p,1)];
model.ub = [bnd(:,2);1/eps*ones(2*n,1);ones(p,1)];

% 'B' : int code 66
% 'C' : int code 67
model.vtype = char([67*ones(1,p+2*n) 66*ones(1,p)]); 

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
model.obj = [zeros(p,1);QR_obj_v;lamda*ones(p,1)];
model.A = sparse([temp1;temp2;temp3;temp4]);
model.rhs = [y;zeros(2*p,1);k0];


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

