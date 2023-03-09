
% Code for replicating numerical results reported in the paper:
% Chen Le-Yu and Lee Sokbae (2023), Sparse Quantile Regression, forthcoming at the Journal of
% Econometrics

% This function performs the FO based L0-penalized quantile regression. 

% function input :
% y : the outcome vector
% x : (n by p) matrix of the covariate data 
% tau : the targeted quantile level 
% lamda : the penalty parameter given by (5.2) of the paper      
% bnd   : (p by 2) matrix where the first and second columns  
%         respectively store the lower and upper bounds 
%         of the unknown coefficients
% tuning : the tuning constant c stated in (5.2) of Chen and Lee (2023)
% K0     : maximal number of variables allowed to be selected 

% function output :
% bhat  : the estimated coefficients of the L0-penalized quantile regression
% obj_v : the value of the L0-PQR objective function
% rtime : the CPU time used to finish the estimation procedure


function [bhat,obj_v,rtime] = L0_PQR_FO_k0(y,x,tau,lamda,bnd,tuning,k0)

tStart = tic;

[n,p]=size(x);

bhat=zeros(p,1);
maxiter=1000; RUNS=50; 
tol = 10^-4;
delta = 2*tol/((max([tau;1-tau]))^2);
objvals=zeros(RUNS,1); 
beta_array = zeros(RUNS,p); 
Lipsch=(norm(x,'fro')^2)/(n*delta);


B=abs(bnd(1,2));

[beta,~] = L1_penalized_QR(y,x,tau,bnd,tuning);

for i = 1: RUNS
   
for j = 1: maxiter    
res = (y - x*beta)/delta;
what = res.*((res<=tau)&(tau-1<=res)) + (tau-1).*(res<tau-1) + tau.*(res>tau);
grd = -x'*what/n;
beta_old=beta;
t=beta-grd/Lipsch;

B_check = (B*B-2*B*t+lamda<0);
beta = B*B_check.*(t>B) + t.*((t.^2)>lamda).*(abs(t)<=B) - B*B_check.*(t<-B);

q=sum(abs(beta)>1e-6);
if (q>k0) && (k0 < p)
[~,IDS] = sort(abs(beta),'descend');
id=IDS(k0+1:end);
beta(id)=zeros(length(id),1);
end 

Q_old = Qn(beta_old,delta,y,x,tau,lamda);
Q=Qn(beta,delta,y,x,tau,lamda);
if ( (abs(Q - Q_old) < tol)&&(j>1))
%disp(['converged in:' num2str(j) ' steps for RUNS ' num2str(i)]);    
break;
end    
end    

uhat=y-x*beta;
objvals(i) = mean(uhat.*(tau-(uhat<=0)))+lamda*sum(abs(beta)>1e-6);
beta_array(i,:) = beta;


z= (abs(beta)>1e-5);

[tmp,~,~,~,~]=QR_grb(y,x(:,z),tau,0,0,bnd(z,:));

beta(z) = tmp; beta(~z) = 0;
end    

[obj_v,ind ] = min(objvals); bhat = beta_array(ind,:); rtime=toc(tStart);

end

function Q=Qn(beta,delta,y,x,tau,lamda)
res = (y - x*beta)/delta;
what = res.*((res<=tau)&(tau-1<=res)) + (tau-1).*(res<tau-1) + tau.*(res>tau);
Q = ((what'*(y-x*beta)-0.5*delta*(what'*what))/length(y))+lamda*sum(abs(beta)>1e-6);
end