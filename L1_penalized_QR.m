
% Code for replicating numerical results reported in the paper:
% Chen Le-Yu and Lee Sokbae (2023), Sparse Quantile Regression, forthcoming at the Journal of
% Econometrics

% This function performs the L1-penalized quantile regression of Belloni and Chernozhukov (2011). 

% function input :
% y : the outcome vector
% x : (n by p) matrix of the covariate data 
% tau : the targeted quantile level    
% bnd   : (p by 2) matrix where the first and second columns  
%         respectively store the lower and upper bounds 
%         of the unknown coefficients
% lamda : the tuning constant c_BC stated in (5.1) of Chen and Lee (2023) 

% function output :
% bhat  : the estimated coefficients of the L1-penalized quantile regression
% rtime : the CPU time used to finish the estimation procedure


function [bhat,rtime] = L1_penalized_QR(y,x,tau,bnd,lambda)

tStart = tic;

[n,p]=size(x);
NormXX=sqrt(diag(x'*x)/n);

 NumSim = max(5000,n);
    U = rand(n, NumSim );
    lambda_temp = zeros(NumSim,1);

for k = 1 : NumSim     
    lambda_temp(k) = max( abs((x'*( (U(:,k)<tau) - tau ) )./NormXX) )/n;        
end

if (nargin == 4)
   lambda_p = quantile(lambda_temp, 0.9  )*NormXX; 
else
    lambda_p = lambda*quantile(lambda_temp, 0.9  )*NormXX; 
end

bhat=zeros(p,1);

model.modelsense = 'min';
model.sense = repmat('=',1,n+p);
%model.vtype = repmat('C',1,2*n+3*p);
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
model.obj = [zeros(p,1);QR_obj_v;lambda_p;lambda_p];

model.A = sparse([temp1;temp2]);
model.rhs = [y;zeros(p,1)];


try
    result = gurobi(model, params);
    bhat=result.x(1:p);
    %rtime=result.runtime;
    rtime=toc(tStart);   
   % fprintf('Optimization returned status: %s\n', result.status);
  
catch gurobiError
    fprintf('Error reported\n');
end

end


function [ lambda_final ] = L1QR_simulate_lambda ( XX, tau, NumSim, n, NormXX )

NumSim = max(NumSim, n);

[ Numrows, NumColumns ] = size( XX );

U = rand(Numrows, NumSim );
lambda = zeros(NumSim,1);

for k = 1 : NumSim     
    lambda(k) = max( abs((XX'*( (U(:,k)<tau) - tau ) )./NormXX) )/n;        
end

lambda_final = 2*quantile(lambda, 0.9  ); 

end



