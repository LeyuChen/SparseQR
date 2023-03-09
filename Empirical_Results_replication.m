
% Code for replicating numerical results reported in the paper:
% Chen Le-Yu and Lee Sokbae (2023), Sparse Quantile Regression, forthcoming at the Journal of
% Econometrics

% These codes can be used to replicate empirical results of the paper.

% change covariate_spec here to configurate the covariate specification   
% used in the empirical study

covariate_spec = 1; % 1 for replicating the results with p=21
                    % 2 for replicating the results with p=49
                    % 3 for replicating the results with p=609
                    % 4 for replicating the results with p=1281
                    % 5 for replicating the results with p=1617

% change method here to replicate the empirical results for a given
% quantile regression method 

method =8;  % 1 for L1-PQR; 
            % 2 for L0-CQR; 
            % 3 for MIO based L0-PQR; 
            % 4 for FO based L0-PQR;
            % 5 for adaptive Lasso PQR with SCAD penalty 
            % 6 for adaptive Lasso PQR with MCP penalty
            % 7 for SCAD based non-convex PQR
            % 8 for MCP based non-convex PQR

clearvars -except covariate_spec method;

rng(1,'twister');

switch covariate_spec  
    case 1
    load('Cattaneo_data_1.mat');     
    case 2
    load('Cattaneo_data_2.mat');
    knstep = 4;  % number of interior knots
    interaction = 0;
    case 3
    load('Cattaneo_data_2.mat');
    knstep = 4;  % number of interior knots
    interaction = 1;
    case 4
    load('Cattaneo_data_2.mat');
    knstep =12;  % number of interior knots
    interaction = 1;
    case 5
    load('Cattaneo_data_2.mat');
    knstep =16;  % number of interior knots
    interaction = 1;
end

y=data(:,1)/1000; x = [ones(size(data,1),1) data(:,2:end)];
data = [y x(:,1) zscore(x(:,2:end),1)];

if covariate_spec~=1  
B_terms = knstep + 3;
B_x = zeros(size(x,1),4*B_terms);
for k=1:4

      xmin = min(x(:,k+21));
      xmax = max(x(:,k+21));
      knq_x = (1/knstep:1/(knstep+1):1/knstep+(knstep-1)/(knstep+1))';     
      kt_x = quantile(x(:,k+21),knq_x);   % equi-spaced in terms of sample quantiles
      B_x(:,(1+B_terms*(k-1)):(B_terms*k))=B_splines(x(:,k+21),3,kt_x,xmin,xmax,0);  
end

temp=[x(:,1:21) B_x];

if interaction == 1
for j=1:20
temp=[temp x(:,j+1).*B_x];     
end    
end

data = [y temp(:,1) zscore(temp(:,2:end),1)]; clear temp;

end

R=10;

n = size(data,1); p =size(data,2)-1; T=300; tau=0.1;

tau_L=tau/2; tau_U = 1 - tau_L;

partition_ind = zeros(R,n); psize = floor(n/4);
for j=1:R
partition_ind(j,:)=randperm(n);       
end       

bhat_L=zeros(p,R); bhat_U=zeros(p,R); tol=1e-5;

bsel_L = zeros(p,1);bsel_U = zeros(p,1);

q=min([p;25]); k0=min([p;100]); thresh_size=min([psize;100]);  
tuning_ext = (0.7.^(1:8))*0.1;

if method == 2   % tuning for L0-CQR;
tuning=1:q;
else             % tuning for penalized QR;  
 if p < thresh_size  
 tuning=[0:0.1:2 tuning_ext];
 else
 tuning=[0.1:0.1:2 tuning_ext];
end
end

Linfty_bound=10;
bnd=[-Linfty_bound*ones(p,1) Linfty_bound*ones(p,1)];    

if method ~= 1   
rate=log(p)/psize;    
gap=zeros(R,2); % MIO gap
rtime=zeros(R,2); % MIO running time
ncount=zeros(R,2); % MIO node count
obj_v=zeros(R,2); % MIO score
tuning_gap=zeros(R,length(tuning),2); % MIO gap
tuning_ncount=zeros(R,length(tuning),2); % MIO node count
tuning_obj_v=zeros(R,length(tuning),2); % MIO score
warm_start_time=zeros(R,length(tuning),2);
end

in_risk=zeros(R,2); % in-sample risk at the estimated parameter vector
out_risk=zeros(R,2); % out-of-sample risk at the estimated parameter vector

num_sel=zeros(R,2);

opt_tuning=zeros(R,2); bhat_tuning=zeros(p,R,length(tuning),2);
tuning_risk=zeros(R,length(tuning),2); tuning_rtime=zeros(R,length(tuning),2);
alpha = (1-tau)*(1+1/psize); bias_corr = zeros(R,1); 
CI = zeros(R,n-3*psize,2); coverage = zeros(R,1);
CI_length = zeros(R,1); 

run('empirical_results_tuning_std.m');
run('empirical_results_tuning_ext.m');

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
     
        for j=1:length(tuning)
              
    uhat_test_L = y2-datax2*bhat_tuning(:,i,j,1);
    uhat_test_U = y2-datax2*bhat_tuning(:,i,j,2);
    tuning_risk(i,j,1) = mean(uhat_test_L.*(tau_L-(uhat_test_L<=0)));
    tuning_risk(i,j,2) = mean(uhat_test_U.*(tau_U-(uhat_test_U<=0)));
      
        end
        
        [~,opt_tuning(i,1)] = min(tuning_risk(i,:,1));
        [~,opt_tuning(i,2)] = min(tuning_risk(i,:,2));
        
      bhat_L(:,i)=bhat_tuning(:,i,opt_tuning(i,1),1); 
      bhat_U(:,i)=bhat_tuning(:,i,opt_tuning(i,2),2); 

uhat_L = y1-datax1*bhat_L(:,i); uhat_U = y1-datax1*bhat_U(:,i);
supp_vec_L=(abs(bhat_L(:,i))>tol); supp_vec_U=(abs(bhat_U(:,i))>tol);
num_sel(i,1)=sum(supp_vec_L); num_sel(i,2)=sum(supp_vec_U);

in_risk(i,1) = mean(uhat_L.*(tau_L-(uhat_L<=0)));
in_risk(i,2) = mean(uhat_U.*(tau_U-(uhat_U<=0)));

  u3_L = y3-datax3*bhat_L(:,i);  u3_U = y3-datax3*bhat_U(:,i); 
  out_risk(i,1) = mean(u3_L.*(tau_L-(u3_L<=0)));
  out_risk(i,2) = mean(u3_U.*(tau_U-(u3_U<=0)));
bias_corr(i) = quantile(max([-u3_L u3_U],[],2),alpha); 
CI(i,:,:) = [datax4*bhat_L(:,i)-bias_corr(i) datax4*bhat_U(:,i)+bias_corr(i)]; 
coverage(i) = mean((y4>=(CI(i,:,1)'))&(y4<=(CI(i,:,2)')));
CI_length(i) = mean(CI(i,:,2)-CI(i,:,1),2);

end

switch method   
    case 1
        QR_method = 'L1-PQR';
    case 2
        QR_method = 'L0-CQR';
    case 3
        QR_method = 'MIO based L0-PQR';
    case 4 
        QR_method = 'FO based L0-PQR';
    case 5
        QR_method = 'AL-SCAD';
    case 6
        QR_method = 'AL-MCP';
    case 7
        QR_method = 'QR-SCAD';
    case 8
        QR_method = 'QR-MCP';
end 

avg_out_risk=mean(out_risk);
avg_num_sel=mean(num_sel);

disp(['Empirical Results for ' QR_method ' with p = ' num2str(p)]);
disp(['Out-of-sample quantile prediction risk at quantile level 0.05: ' num2str(avg_out_risk(1))]);
disp(['Out-of-sample quantile prediction risk at quantile level 0.95: ' num2str(avg_out_risk(2))]);
disp(['Estimated sparsity at quantile level 0.05: ' num2str(avg_num_sel(1))]);
disp(['Estimated sparsity at quantile level 0.95: ' num2str(avg_num_sel(2))]);
disp(['Confidence interval length: ' num2str(mean(CI_length))]);
disp(['Confidence interval coverage: ' num2str(mean(coverage))]);

run('variable_selection_results.m');
