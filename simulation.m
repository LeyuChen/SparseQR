
bhat=zeros(p,R);  

Linfty_bound=10;
bnd=[-Linfty_bound*ones(size(bhat,1),1) Linfty_bound*ones(size(bhat,1),1)];    

if method ~= 1  
rate=log(p)/N;    
gap=zeros(R,1); % MIO gap
rtime=zeros(R,1); % MIO running time
ncount=zeros(R,1); % MIO node count
obj_v=zeros(R,1); % MIO score
tuning_gap=zeros(R,length(tuning)); % MIO gap
tuning_ncount=zeros(R,length(tuning)); % MIO node count
tuning_obj_v=zeros(R,length(tuning)); % MIO score
warm_start_time=zeros(R,length(tuning));
end

DGP_risk=zeros(R,1); % in-sample risk at the DGP parameter vector
val_risk=zeros(R,1); % in-sample risk at the estimated parameter vector
DGP_risk_test=zeros(R,1); % out-of-sample risk at the DGP parameter vector
val_risk_test=zeros(R,1); % out-of-sample risk at the estimated parameter vector

sel = zeros(R,1); sel_all = zeros(R,1);
num_irrel=zeros(R,1); num_sel=zeros(R,1);

bias_norm = zeros(R,1); Reg_Risk = zeros(R,1); Reg_Risk_test = zeros(R,1);
bias = zeros(p,R); mse = zeros(R,1);

aux = (beta~=0); sparsity=sum(aux);
tol=1e-5; 

opt_tuning=zeros(R,1); bhat_tuning=zeros(p,R,length(tuning));
tuning_risk=zeros(R,length(tuning)); tuning_rtime=zeros(R,length(tuning));

for i=1:R

y=data(1:N,1,i);
datax=data(1:N,2:end,i);

y_test=data(N+1:2*N,1,i);
datax_test=data(N+1:2*N,2:end,i);

if method > 2
lamda = tuning*mean(abs(y))*rate;      
end    
      
        for j=1:length(tuning)
            disp(['q : ' num2str(tuning(j))]);
  
  if method == 1 % L1-PQR    
[bhat_tuning(:,i,j),tuning_rtime(i,j)] = L1_penalized_QR(y,datax,tau,bnd,tuning(j));       
  elseif method == 2 % L0-CQR
   
[bhat_tuning(:,i,j),tuning_obj_v(i,j),tuning_gap(i,j),tuning_rtime(i,j),tuning_ncount(i,j)] = best_subset_LAD(y,datax,tuning(j),T,0,Linfty_bound);

 else
  try
  switch method    
      case 3  % MIO based L0-PQR 
 [bhat_tuning(:,i,j),tuning_obj_v(i,j),tuning_gap(i,j),tuning_rtime(i,j),tuning_ncount(i,j),warm_start_time(i,j)] = L0_penalized_QR_warm_k0(y,datax,tau,lamda(j),T,0,bnd,tuning(j),k0);
      case 4 % FO based L0-PQR
 [bhat_tuning(:,i,j),tuning_obj_v(i,j)] = L0_PQR_FO_k0(y,datax,tau,lamda(j),bnd,tuning(j),k0);
      case 5 % adaptive Lasso PQR with SCAD penalty 
 [bhat_tuning(:,i,j),tuning_rtime(i,j)] = QR_Adaptive_lasso(y,datax,tau,bnd,tuning(j),1);
      case 6 % adaptive Lasso PQR with MCP penalty
 [bhat_tuning(:,i,j),tuning_rtime(i,j)] = QR_Adaptive_lasso(y,datax,tau,bnd,tuning(j),0);
      case 7 % SCAD based non-convex PQR
 [bhat_tuning(:,i,j),tuning_rtime(i,j)] = QR_SCAD_MCP(y,datax,tau,bnd,tuning(j),1);
      case 8 % MCP based non-convex PQR
 [bhat_tuning(:,i,j),tuning_rtime(i,j)] = QR_SCAD_MCP(y,datax,tau,bnd,tuning(j),0);         
  end
 
 
catch gurobiError
    fprintf('Error reported\n');
 end
        end
  
    uhat_test = y_test-datax_test*bhat_tuning(:,i,j);
    tuning_risk(i,j) = mean(uhat_test.*(tau-(uhat_test<=0)));
        end
        
        [~,opt_tuning(i)] = min(tuning_risk(i,:));
      disp(['optimal tuning : ' num2str(opt_tuning(i))]);
      bhat(:,i)=bhat_tuning(:,i,opt_tuning(i));    


u0 = y-datax*beta; uhat = y-datax*bhat(:,i);
supp_vec=(abs(bhat(:,i))>tol);
disp(['L0 norm :' num2str(sum(supp_vec))]);

DGP_risk(i) = mean(u0.*(tau-(u0<=0)));
val_risk(i) = mean(uhat.*(tau-(uhat<=0)));

bias(:,i) = bhat(:,i)-beta;
bias_norm(i) = norm(bias(:,i));
mse(i) = bias(:,i)'*bias(:,i);
Reg_Risk(i) = bias(:,i)'*(datax'*datax)*bias(:,i)/N;

if sum(supp_vec(aux))==sparsity
 sel(i) = 1;
end

  if aux == supp_vec
 sel_all(i) = 1;    
 end
 
 num_irrel(i) = sum(abs(bhat(~aux,i))>tol);
 num_sel(i)=sum(supp_vec);

  if N_val>0
  
  y_val=data(2*N+1:end,1,i);
  datax_val=data(2*N+1:end,2:end,i);
   
  u0_val = y_val-datax_val*beta;
  uhat_val = y_val-datax_val*bhat(:,i);
  DGP_risk_test(i) = mean(u0_val.*(tau-(u0_val<=0)));
  val_risk_test(i) = mean(uhat_val.*(tau-(uhat_val<=0)));
  Reg_Risk_test(i) = bias(:,i)'*(datax_val'*datax_val)*bias(:,i)/N_val;
 end

end



