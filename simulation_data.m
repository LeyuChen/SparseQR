
function [y,datax] = simulation_data(n,beta,sigma,v)

p=length(beta);

e=v*randn(n,1);

datax = [ones(n,1) mvnrnd(zeros(1,p-1),sigma,n)];
h=datax(:,2);
y = datax*beta + h.*e;
end