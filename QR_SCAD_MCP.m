
% Code for replicating numerical results reported in the paper:
% Chen Le-Yu and Lee Sokbae (2023), Sparse Quantile Regression, forthcoming at the Journal of
% Econometrics

% This function performs the non-convex penalized quantile
% regression of Wang, Wu and Li (2012).

% function input :
% y : the outcome vector
% x : (n by p) matrix of the covariate data 
% tau : the targeted quantile level    
% bnd   : (p by 2) matrix where the first and second columns  
%         respectively store the lower and upper bounds 
%         of the unknown coefficients
% tuning : the tuning constant mu stated in (5.6) of Chen and Lee (2023)
% type : the type of penalty weight 
% set type = 1 for SCAD
% set type = 0 for MCP

% function output :
% bhat  : the estimated coefficients of the non-convex penalized quantile regression
% rtime : the CPU time used to finish the estimation procedure

function [bhat,rtime] = QR_SCAD_MCP(y,x,tau,bnd,tuning,type)

tStart = tic;

[n,p]=size(x);

[bhat_ini,~] = L1_penalized_QR(y,x,tau,bnd,tuning);

a = 3.7;
lamda = tuning;
weight = SCAD_MCP(abs(bhat_ini),a,lamda,type); 

bhat=zeros(p,1);

model.modelsense = 'min';
model.sense = repmat('=',1,n+p);
model.lb = [bnd(:,1);zeros(2*(n+p),1)];
model.ub = [bnd(:,2);1/eps*ones(2*n,1);bnd(:,2);bnd(:,2)];

tol=1e-6;

params.outputflag = 0; 
params.OptimalityTol=tol;
params.FeasibilityTol=tol;
params.IntFeasTol=tol;

QR_obj_v = [tau*ones(n,1);(1-tau)*ones(n,1)]/n;
temp1 = [x eye(n) -eye(n) zeros(n,2*p)];
temp2 = [eye(p) zeros(p,2*n) -eye(p) eye(p)];
model.A = sparse([temp1;temp2]);
model.rhs = [y;zeros(p,1)];
obj_temp = [zeros(p,1);QR_obj_v];
maxiter=100; tol = 10^-4;

for j = 1: maxiter    
model.obj = [obj_temp;weight;weight];
    
try
    result = gurobi(model, params);
    bhat=result.x(1:p);    
    %fprintf('Optimization returned status: %s\n', result.status);
catch gurobiError
    fprintf('Error reported\n');
end
old_w = weight;
weight = SCAD_MCP(abs(bhat),a,lamda,type); 
if (norm(old_w-weight) < tol)
disp(['converged in:' num2str(j) ' steps']);    
break;
end 
end
rtime=toc(tStart);
end

