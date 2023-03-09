function  bet_LAD = FUN_LAD_grb( X, Y) 
%% computes the LAD soln via grb


[n,p] =size(X);  

AA= [ X , eye(n,n)  ]; %% taui  + Xbet >= Y 
AA1= [ -X, eye(n,n) ]; %% taui - Xbet >= -Y 
% vars (bet, taui): p . n 

A=sparse([AA; AA1]);

clear model

model.A=A;
model.lb = [ repmat(-Inf,p,1)  ; repmat(0,n,1)];

model.rhs= [ Y ; -Y ] ;
model.sense= repmat('>', 2*n  , 1); 
params.outputflag = 0;

model.obj = [ zeros(p, 1) ; ones(n , 1) ];

results=gurobi(model,params);

bet_LAD= results.x(1:p);

