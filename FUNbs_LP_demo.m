function  [results2, betavec, zvec ]  = FUNbs_LP_demo( X,Y, kk, sideInfo , params) 

[nn,pp ] = size(X);


if(isfield(sideInfo,'beta_start')) 
beta_start = sideInfo.beta_start;
else
beta_start=[];
end


if(isfield(sideInfo,'Linfty_bound')) 
Linfty_bound = sideInfo.Linfty_bound ;
end


if isempty(beta_start)
start_warm=(1<0);
else
start_warm=(1>0);

%% supply model warm-start:

betavec=beta_start; 
zvec = 0*beta_start;
zvec(beta_start==0)=1;
tvec = abs(  Y - X*betavec(:) );
muvec= abs(betavec);

end


%% form implied upper and lower bounds:
zetaInfbound= zeros(nn,1);
for ii = 1:nn
aa=sort(abs(X(ii,:)),'descend');
zetaInfbound(ii) = Linfty_bound*sum(aa(1:kk)); %% zeta = Xbeta ; 
clear aa;
end
t_UB = abs(Y) + zetaInfbound;  t_LB = 0*abs(abs(Y) - zetaInfbound); 
%% mu_UB = L1betabound*ones(pp,1) ; mu_LB = zeros(pp,1);

BigM = Linfty_bound;

%% sum (abs(  Y - Xbeta ))  sbt nnz(beta) <= k; sum(abs(beta)) <= Linfty_bound*k
%% t_{i} >= Y_{i} - x_{i}'beta >= -t_{i}
%% beta<= Mz_{i} ; -beta_{i} <= Mz_{i} 
%% vars p + p + n = (beta, z, t)

numvars= 2*pp + nn;
 
AA= [X , zeros(nn,pp) , -eye(nn)];
AA1=[-X, zeros(nn,pp), -eye(nn)];
lastrow=[zeros(1, pp), ones(1,pp), zeros(1,nn) ];
AA = sparse([AA;AA1; lastrow]); clear AA1

AA1 = [ eye(pp) , -BigM*eye(pp) , zeros(pp,nn)]; %% beta - BigM*z <=0 ; 
AA2 = [ -eye(pp) , -BigM*eye(pp) , zeros(pp,nn)]; %% -beta  - BigM*z <=0

A= [ AA ; AA1 ; AA2 ];

clear model

model.modelsense = 'min';
model.A=A;

%% vars p + p + n  = (beta, z, t)

model.sense= [repmat('<',size(A ,1),1) ];

model.rhs=[Y; -Y; kk; zeros(2*pp,1)];

%% vars p + p + n  = (beta, z, t)

model.ub = [repmat(Linfty_bound,pp,1); repmat(1,pp,1); t_UB ] ;

model.lb = [repmat(-Linfty_bound,pp,1); repmat(0,pp,1); t_LB] ;


model.obj= [zeros( 2*pp ,1); ones(nn,1)];

model.vtype=[ repmat('C',pp,1) ; repmat('B', pp ,1); repmat('C',nn,1) ] ;


if (start_warm == (1>0))
model.start=[betavec;zvec; tvec];
end



if (isempty(params))
results2 = gurobi(model);
else
results2=gurobi(model,params);
end



%% vars p + p + n  + p = (beta, z, t , mu)

betavec = results2.x(1:pp);
zvec=results2.x((1+pp) : (2*pp) );


