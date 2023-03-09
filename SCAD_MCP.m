
% Code for replicating numerical results reported in the paper:
% Chen Le-Yu and Lee Sokbae (2023), Sparse Quantile Regression, forthcoming at the Journal of
% Econometrics

% The SCAD and MCP penalty function derivatives
% For SCAD, a > 2.
% For MCP, a > 1.

function [weight] = SCAD_MCP(beta,a,lamda,type)

if lamda>0
if type ==1 % SCAD derivative
weight = lamda*(beta<=lamda) + (beta>lamda).*max(a*lamda-beta,0)/(a-1);
else % MCP derivative
temp=a*lamda;
weight = lamda*(1-beta/temp).*(beta<=temp);
end
else
weight = zeros(length(beta),1);
end
end
