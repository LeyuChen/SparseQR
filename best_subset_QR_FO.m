function  [bet_best,OBJVALS]  = best_subset_QR_FO(y,x,tau,k,bnd) 

p=size(x,2);

maxiter=1000; RUNS=50; TOL = 10^-4;
%% counter for obj-vals
objvals=zeros(maxiter,1); 

OBJVALS=zeros(RUNS,1);  %% stores the best obj-val for each run
bet_array = zeros(RUNS,p); %% array of regression coefficients

betak= zeros(p,1);

hub_rho = 5;
Lipsh= (norm(x,'fro')^2)/hub_rho;

for start = 1: RUNS

%% continuation

hub_rho = hub_rho*0.8;
Lipsh = Lipsh/0.8 ;

for ii = 1: maxiter

res = (y - x*betak)/hub_rho;
what = res.*((res<=tau)&(tau-1<=res)) + (tau-1).*(res<tau-1) + tau.*(res>tau);
grad = -x'*what;

betakold=betak;

vec= betak - grad/Lipsh;

[~,IDS] = sort(abs(vec),'descend');
ids=IDS(1:k);
betak=zeros(p,1);
betak(ids) = vec(ids);


objvals(ii) = QR_obj(betakold,hub_rho,y,x,tau);

if ( (norm(betak - betakold)/max(norm(betak),1) < TOL)&&(ii>1))
break;
end

end

%% ii
z0= (abs(betak)>1e-5);
[tmp,~,~,~,~]=QR_grb(y,x(:,z0),tau,0,0,bnd(z0,:));
betak(z0) = tmp; betak(~z0) = 0;
uhat = y-x*betak;  
OBJVALS(start)= sum(uhat.*(tau-(uhat<=0)));
bet_array(start,:) = betak;

end

[~,cc ] = min(OBJVALS);

bet_best = bet_array(cc(1),:); bet_best = bet_best(:);

end

function Q=QR_obj(beta,delta,y,x,tau)
res = (y - x*beta)/delta;
what = res.*((res<=tau)&(tau-1<=res)) + (tau-1).*(res<tau-1) + tau.*(res>tau);
Q = what'*(y-x*beta)-0.5*delta*(what'*what);
end
