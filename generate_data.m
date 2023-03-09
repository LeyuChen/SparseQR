
function [data]=generate_data(N,N_val,R,beta,rho,v)

rng(9999,'twister');

p=length(beta);

sigma=ones(p-1,1);
for i=1:p-2
sigma(i+1)=rho^i;
end
sigma=toeplitz(sigma);

data = zeros(2*N+N_val,p+1,R); 

for i=1:R
[y,datax] = simulation_data(N,beta,sigma,v); % training sample
data(1:N,:,i)=[y datax];
[y_test,datax_test] = simulation_data(N,beta,sigma,v); % testing sample   
data(N+1:2*N,:,i)=[y_test datax_test];  
[y_val,datax_val] = simulation_data(N_val,beta,sigma,v); % validation sample
data(2*N+1:end,:,i)=[y_val datax_val];  
end

end