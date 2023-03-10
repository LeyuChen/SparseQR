
% Code for replicating numerical results reported in the paper:
% Chen Le-Yu and Lee Sokbae (2023), Sparse Quantile Regression, forthcoming at the Journal of
% Econometrics

%**********************************************************%
%******           B-splines                         *******%
%**********************************************************%

%{
Input:
    x = (n*1) vector of data
  deg = degree (order = df + 1)
knots = interior knots
 xmin = lower bound
 xmax = upper bound
intercept = 1 (intercept included); = 0 (intercept excluded)

Output:
(n*kappa) matrix of B-spline basis, where kappa = #(knots)+order.
%}

function m1tmp=B_splines(x,deg,knots,xmin,xmax,intercept)

       k = deg + 1;                               % order 
       t = [xmin;knots;xmax];                     % an initial knot sequence 
     tr1 = trimr(t,0,1);                          % a vector of t[i]         
     tr2 =  trimr(t,1,0);                         % a vector of t[i+1]       
   m1tmp = (x >= tr1').*(x < tr2');               % B-splines of degree 0   

if k > 1

j = 2;

while j <= k

      t = [xmin*ones(j,1);knots;xmax*ones(j,1)];               % an updated knot sequence 
  m1tmp = [zeros(size(m1tmp,1),1) m1tmp zeros(size(m1tmp,1),j)];
      
    tr1 = trimr(t,0,j);
    tr2 = trimr(t,j-1,1);
    m11 = trimr(m1tmp',0,j)';
   den1 = (tr2 - tr1)';
   den1 = (1-(den1==0))./(den1 + (den1==0));               
  term1 = (x-tr1').*m11.*den1;

    tr3 = trimr(t,j,0);
    tr4 = trimr(t,1,j-1);
    m12 = trimr(m1tmp',1,j-1)';
   den2 = (tr3 - tr4)';
   den2 = (1-(den2==0))./(den2 + (den2==0));
  term2 = (tr3'-x).*m12.*den2;

  m1tmp = term1 + term2;                                     

j = j + 1;
end

end

if intercept == 0
m1tmp = m1tmp(:,2:size(m1tmp,2));
end

end

function t=trimr(x,t,b)
t=x(t+1:size(x,1)-b,:);   
end
