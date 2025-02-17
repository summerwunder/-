function [u,y_sp,ufull,y_sp_full,ti] = CtrlSegment(uini,yini,Fcasts,umax,umin,u0,Ctrlparams)
global catcher

T = Ctrlparams.T;
Tini = Ctrlparams.Tini;
Tf = Ctrlparams.Tf;
lamg = Ctrlparams.lamg;
lamU = 1;
lamcomf = 1;
lamini = 1;

ysp = Fcasts.ysp_lo;
ysp_up = Fcasts.ysp_hi;
cost = Fcasts.cost;
ti=0;

num_g = (T-Tini-Tf+1);
folds = ceil(length(cost)/Tf);

Up=Ctrlparams.Up;
Uf=Ctrlparams.Uf;
Yp=Ctrlparams.Yp;
Yf=Ctrlparams.Yf;

costlong = cost;

costvect = zeros(1,folds*num_g);
for ii=1:folds
    costvect((ii-1)*num_g+1:ii*num_g) = costlong(1:Tf)'*Uf;
end
spvect = zeros(numel(ysp),1);
spvect_up = zeros(numel(ysp),1);
for ii=1:folds
    spvect((ii-1)*Tf+1:ii*Tf) = reshape(ysp((ii-1)*Tf+1:ii*Tf,:),Tf,1);
    spvect_up((ii-1)*Tf+1:ii*Tf) = reshape(ysp_up((ii-1)*Tf+1:ii*Tf,:),Tf,1);
end

%% Formulate constraints
% Up*g = uini:
nrows = Tini*folds;
Aeq1 = [kron(diag(ones(folds-1,1),-1),-Uf)+kron(eye(folds),Up),zeros(nrows,Tf*folds),zeros(nrows,Tini*folds)];
beq1 = [uini;zeros(Tini*(folds-1),1)];

% Uf*g<=umax; -Uf*g<=-umin
nrows = Tf*folds;
Aineq1 = [kron(eye(folds),Uf),zeros(nrows,Tf*folds),zeros(nrows,Tini*folds);kron(eye(folds),-Uf),zeros(nrows,Tf*folds),zeros(nrows,Tini*folds)];
bineq1 = [repmat(umax*ones(Tf,1),folds,1);repmat(-umin*ones(Tf,1),folds,1)];

% Yf*g-eps<=Tsp_hi; -Yf*g-eps2<=-Tsp_lo
nrows = Tf*folds;
Aineq2 = [kron(eye(folds),Yf),-eye(nrows),zeros(nrows,Tini*folds);-kron(eye(folds),Yf),-eye(nrows),zeros(nrows,Tini*folds)];
bineq2 = [spvect_up;-spvect];

% Yp*g-eps2<=yini; -Yp*g-eps2<=-yini
nrows = Tini*folds;
Aineq3 = [kron(diag(ones(folds-1,1),-1),-Yf)+kron(eye(folds),Yp),zeros(nrows,Tf*folds),-eye(nrows);kron(diag(ones(folds-1,1),-1),Yf)+kron(eye(folds),-Yp),zeros(nrows,Tf*folds),-eye(nrows)];
bineq3 = [yini;zeros(Tini*(folds-1),1);-yini;zeros(Tini*(folds-1),1)];

% -eps<=0
nrows = Tf*folds;
Aineq4 = [zeros(nrows,num_g*folds),-eye(nrows),zeros(nrows,Tini*folds)];
bineq4 = zeros(nrows,1);

% -eps2<=0
nrows = Tini*folds;
Aineq5 = [zeros(nrows,num_g*folds),zeros(nrows,Tf*folds),-eye(nrows)];
bineq5 = zeros(nrows,1);

Aineq = sparse([Aineq1;Aineq2;Aineq3;Aineq4;Aineq5]);
bineq = sparse([bineq1;bineq2;bineq3;bineq4;bineq5]);
Aeq = sparse(Aeq1);
beq = sparse(beq1);

dropto=1;

f0 = [0*lamU*costvect,kron(linspace(1,dropto,folds),0*lamcomf*ones(1,Tf)),kron(linspace(1,dropto,folds),lamini*ones(1,Tini))];
f1 = [0*lamU*costvect,kron(linspace(1,dropto,folds),1*lamcomf*ones(1,Tf)),kron(linspace(1,dropto,folds),0*lamini*ones(1,Tini))];
f2 = [1*lamU*costvect,kron(linspace(1,dropto,folds),0*lamcomf*ones(1,Tf)),kron(linspace(1,dropto,folds),0*lamini*ones(1,Tini))];
H = sparse(diag([1*lamg*ones(1,num_g*folds),0*kron(linspace(1,dropto,folds),0*lamcomf*ones(1,Tf)),kron(linspace(1,dropto,folds),0*lamini*ones(1,Tini))]));

optionsL = optimoptions('linprog','Display','none');
options = optimoptions('quadprog','LinearSolver','sparse','Display','off');
% Objective
try
    [gopt0,acc] = linprog(1*f0+0*f1+0*f2,Aineq,bineq,Aeq,beq,[],[],optionsL);
catch
    catcher=catcher+1;
    acc = 0.01;
end
try
    tic
    [gopt,~] = quadprog(1*H,0*f0+1*f1+0*f2,[Aineq;f0],[bineq;max(0,acc)],Aeq,beq,[],[],[],options);
    ti=toc;
catch
    catcher=catcher+1;
    [gopt,~] = quadprog(1*H,0*f0+1*f1+0*f2,[Aineq;f0],[bineq;0.01],Aeq,beq);
end

    

y_sp_full = kron(eye(folds),Yf(1:Tf,:))*gopt(1:num_g*folds);
ufull = kron(eye(folds),Uf(1:Tf,:))*gopt(1:num_g*folds);



y_sp = y_sp_full(1);
u = ufull(1);


end