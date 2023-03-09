function  [bet_array, OBJVALS, bet_best]  = FUN_bs_LP_FO( X,Y, k, bb0, maxiter, RUNS , TOL) 

%% try best subset via GD:
%% solve 0.5\|Y - Xbeta\|_1 \sbt \| beta\|_0 \leq k


[n,pp]=size(X);


if (isempty(maxiter))
maxiter=1000; 
end

if (isempty(RUNS))
RUNS=50; 
end

if (isempty(TOL))
TOL = 10^-4;
end




%% counter for obj-vals
objvals=zeros(maxiter,1); 
LBvals=objvals; sparsity=objvals;


%% number of RUNS [ every run is a random initialization ]

OBJVALS=zeros(RUNS,1);  %% stores the best obj-val for each run
bet_array = zeros(RUNS,pp); %% array of regression coefficients


%% start from LAD solution [ thresholded ]:
%{
bbLAD = FUN_LAD_grb(X, Y) ; 
[aa,bb]=sort(abs(bbLAD),'descend');
bbLAD(bb((k+1):end)) = 0;
%}

if (isempty(bb0))
betak= zeros(pp,1);
else
betak = bb0;
end


hub_rho = 5;
Lipsh= normest(X)^2/hub_rho;
sub_grad_norm = norm(sum(abs(X),1));


for start = 1: RUNS

%% continuation

hub_rho = hub_rho*0.8;
Lipsh = Lipsh/0.8 ;

zk = (abs(betak)>1e-7);

for ii = 1: maxiter

resk = Y - X*betak;
what= min(abs(resk)/hub_rho,1).*sign(resk) ;
grad = -X'*what;


betakold=betak;

vec= betak - grad/Lipsh;

[vals,IDS] = sort(abs(vec),'descend');
ids=IDS(1:k);
betak=zeros(pp,1);
betak(ids) = vec(ids);


objvals(ii) = resk'*what - hub_rho*norm(what,2)^2/2   ;


if ( (norm(betak - betakold)/max(norm(betak),1) < TOL)&&(ii>1))
break;
end




end

%% ii
z0= (abs(betak)>1e-5);
tmp = FUN_LAD_grb(X(:,z0), Y) ;
betak(z0) = tmp; betak(~z0) = 0;
OBJVALS(start)= norm(Y - X*betak,1);
bet_array(start,:) = betak;

end

[ss,cc ] = min(OBJVALS);

bet_best = bet_array(cc(1),:); bet_best = bet_best(:);


