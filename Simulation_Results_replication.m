
% Code for replicating numerical results reported in the paper:
% Chen Le-Yu and Lee Sokbae (2023), Sparse Quantile Regression, forthcoming at the Journal of
% Econometrics

% These codes can be used to replicate simulation results of the paper.

% change Results here to configurate the simulation setting for 
% replication of the results reported in Tables 5.1 ~ 5.3.

Results = 3; % 1 for replicating the results reported in Table 5.1
             % 2 for replicating the results reported in Table 5.2
             % 3 for replicating the results reported in Table 5.3
           

% change method here to replicate the simulation results for a given
% quantile regression method 

method =8;  % 1 for L1-PQR; 
            % 2 for L0-CQR; 
            % 3 for MIO based L0-PQR; 
            % 4 for FO based L0-PQR;
            % 5 for adaptive Lasso PQR with SCAD penalty 
            % 6 for adaptive Lasso PQR with MCP penalty
            % 7 for SCAD based non-convex PQR
            % 8 for MCP based non-convex PQR

clearvars -except Results method;

switch Results   
    case 1
     p=10; s=5;   
     beta = zeros(p,1); beta(1:floor(p/s):p)=1; 
     configuration=1; conf='i';
    case 2
     p=500; s=5;
     beta = zeros(p,1); beta(1:floor(p/s):p)=1; 
     configuration=1; conf='i';
    case 3
     p=500; s=20;
     beta_star = [ones(5,1);0.5.^(1:s-5)'];
     beta = zeros(p,1); beta(1:floor(p/s):p)=beta_star; 
     configuration=2; conf='ii';
end

N = 100; % size of the training sample
N_val = 5000; % size of the validation sample
R = 100; % simulation repetitions

T=600; tau=0.5; rho=0.5; v=0.25; 

data=generate_data(N,N_val,R,beta,rho,v);

rng(1,'twister');

q=min([p;25]); k0=min([p;100]);

if method == 2   % tuning for L0-CQR;
tuning=1:q;

else             % tuning for penalized QR;  
 if p < N    
 tuning=0:0.1:2;
 else
 if configuration == 1
 tuning=0.1:0.1:2;   
 else
 tuning=0.01:0.01:2; 
 end
 end
 end

run("simulation.m");

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

output_str1 = ["Corr_sel" "Orac_sel" "Num_irrel" "Avg_sparsity" "in_RR" "out_RR"];
disp(['Table 5.' num2str(Results) ': Simulation results for ' QR_method  ' under parameter configuration (' conf ') with p = ' num2str(p)]);
disp(output_str1);
disp(num2str(mean([sel sel_all num_irrel num_sel val_risk./DGP_risk val_risk_test./DGP_risk_test])));
output_str2 = ["average parameter estimation error" " average regression function estimation error"];
disp(output_str2);
disp(mean([bias_norm Reg_Risk_test]));

