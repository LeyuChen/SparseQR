
rng(1,'twister');

partition_ind = zeros(R,n); 
for j=1:R
partition_ind(j,:)=randperm(n);       
end

for i=1:R
disp(i);

ind1 = partition_ind(i,1:psize);
ind2 = partition_ind(i,psize+1:2*psize);
ind3 = partition_ind(i,2*psize+1:3*psize);
ind4 = partition_ind(i,3*psize+1:end);

y1=data(ind1,1); datax1=data(ind1,2:end);
y2=data(ind2,1); datax2=data(ind2,2:end);
y3=data(ind3,1); datax3=data(ind3,2:end);
y4=data(ind4,1); datax4=data(ind4,2:end);

if method > 2
lamda = tuning*mean(abs(y1))*rate;      
end    
      
        for j=length(tuning)-length(tuning_ext)+1:length(tuning)
            disp(['q : ' num2str(tuning(j))]);
            
         
  if method == 1 % L1-PQR    
      
[bhat_tuning(:,i,j,1),tuning_rtime(i,j,1)] = L1_penalized_QR(y1,datax1,tau_L,bnd,tuning(j));
[bhat_tuning(:,i,j,2),tuning_rtime(i,j,2)] = L1_penalized_QR(y1,datax1,tau_U,bnd,tuning(j));  

  elseif method == 2 % L0-CQR
 
[bhat_tuning(:,i,j,1),tuning_obj_v(i,j,1),tuning_gap(i,j,1),tuning_rtime(i,j,1),tuning_ncount(i,j,1)] = best_subset_QR(y1,datax1,tuning(j),T,tau_L,0,bnd);
[bhat_tuning(:,i,j,2),tuning_obj_v(i,j,2),tuning_gap(i,j,2),tuning_rtime(i,j,2),tuning_ncount(i,j,2)] = best_subset_QR(y1,datax1,tuning(j),T,tau_U,0,bnd);
  
 elseif method==3 || method==4
      
      try
          
 if method == 3 % MIO based L0-PQR 
  [bhat_tuning(:,i,j,1),tuning_obj_v(i,j,1),tuning_gap(i,j,1),tuning_rtime(i,j,1),tuning_ncount(i,j,1),warm_start_time(i,j,1)] = L0_penalized_QR_warm_k0(y1,datax1,tau_L,lamda(j),T,0,bnd,tuning(j),k0);
  [bhat_tuning(:,i,j,2),tuning_obj_v(i,j,2),tuning_gap(i,j,2),tuning_rtime(i,j,2),tuning_ncount(i,j,2),warm_start_time(i,j,2)] = L0_penalized_QR_warm_k0(y1,datax1,tau_U,lamda(j),T,0,bnd,tuning(j),k0);
  
 else % FO based L0-PQR 
  [bhat_tuning(:,i,j,1),tuning_obj_v(i,j,1),tuning_rtime(i,j,1)] = L0_PQR_FO_k0(y1,datax1,tau_L,lamda(j),bnd,tuning(j),k0);
  [bhat_tuning(:,i,j,2),tuning_obj_v(i,j,2),tuning_rtime(i,j,2)] = L0_PQR_FO_k0(y1,datax1,tau_U,lamda(j),bnd,tuning(j),k0);

 end
 
      catch gurobiError
    fprintf('Error reported\n');
      end
      
  elseif method == 5 % AL-SCAD
  [bhat_tuning(:,i,j,1),tuning_rtime(i,j,1)] = QR_Adaptive_lasso(y1,datax1,tau_L,bnd,tuning(j),1);
  [bhat_tuning(:,i,j,2),tuning_rtime(i,j,2)] = QR_Adaptive_lasso(y1,datax1,tau_U,bnd,tuning(j),1);
  
  elseif method == 6 % AL-MCP
  [bhat_tuning(:,i,j,1),tuning_rtime(i,j,1)] = QR_Adaptive_lasso(y1,datax1,tau_L,bnd,tuning(j),0);
  [bhat_tuning(:,i,j,2),tuning_rtime(i,j,2)] = QR_Adaptive_lasso(y1,datax1,tau_U,bnd,tuning(j),0);
  
  elseif method == 7 % QR-SCAD
  [bhat_tuning(:,i,j,1),tuning_rtime(i,j,1)] = QR_SCAD_MCP(y1,datax1,tau_L,bnd,tuning(j),1);
  [bhat_tuning(:,i,j,2),tuning_rtime(i,j,2)] = QR_SCAD_MCP(y1,datax1,tau_U,bnd,tuning(j),1);  
  
  else  % QR-MCP
  [bhat_tuning(:,i,j,1),tuning_rtime(i,j,1)] = QR_SCAD_MCP(y1,datax1,tau_L,bnd,tuning(j),0);
  [bhat_tuning(:,i,j,2),tuning_rtime(i,j,2)] = QR_SCAD_MCP(y1,datax1,tau_U,bnd,tuning(j),0);  
  
  end
  
    uhat_test_L = y2-datax2*bhat_tuning(:,i,j,1);
    uhat_test_U = y2-datax2*bhat_tuning(:,i,j,2);
    tuning_risk(i,j,1) = mean(uhat_test_L.*(tau_L-(uhat_test_L<=0)));
    tuning_risk(i,j,2) = mean(uhat_test_U.*(tau_U-(uhat_test_U<=0)));
      
        end
        
end
