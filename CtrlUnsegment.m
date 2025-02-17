function [u,y_sp,ufull,y_sp_full,ti] = CtrlUnsegment(uini,yini,Fcasts,umax,umin,Ctrlparams)
global catcher
T = Ctrlparams.T;
Tini = Ctrlparams.Tini;
Tf = Ctrlparams.Tf;
lamg = Ctrlparams.lamg;
lamU = 1;
lamY = 1;
lamini = 1;


ti=0;
ysp = Fcasts.ysp_lo;
ysp_up = Fcasts.ysp_hi;

num_g = (T-Tini-Tf+1);

Up=Ctrlparams.Up;
Uf=Ctrlparams.Uf;
Yp=Ctrlparams.Yp;
Yf=Ctrlparams.Yf;

cost = ones(Tf,1);
    

%% Formulate constraints
% Variables: [g; eps_y; eps_ini]

% Up*g = uini:
Aeq1 = [Up,zeros(Tini,Tf),zeros(Tini,Tini)];
beq1 = uini;

% Uf*g<=umax; -Uf*g<=-umin
Aineq1 = [Uf,zeros(Tf,Tf),zeros(Tf,Tini);-Uf,zeros(Tf,Tf),zeros(Tf,Tini)];
bineq1 = [umax*ones(Tf,1);-umin*ones(Tf,1)];


% Yf*g-eps_y<=Ysp_hi; -Yf-eps_y<=-Ysp_lo
Aineq2 = [Yf,-eye(Tf),zeros(Tf,Tini);-Yf,-eye(Tf),zeros(Tf,Tini)];
bineq2 = [ysp_up;-ysp];


% Yp*g-eps_ini<=yini; -Yp*g-eps_ini<=-yini
Aineq3 = [Yp,zeros(Tini,Tf),-eye(Tini);-Yp,zeros(Tini,Tf),-eye(Tini)];
bineq3 = [yini;-yini];

% -eps_y<=0
Aineq4 = [zeros(Tf,num_g),-eye(Tf),zeros(Tf,Tini)];
bineq4 = zeros(Tf,1);

% -eps_ini<=0
Aineq5 = [zeros(Tini,num_g),zeros(Tini,Tf),-eye(Tini)];
bineq5 = zeros(Tini,1);

Aineq = sparse([Aineq1;Aineq2;Aineq3;Aineq4;Aineq5]);
bineq = sparse([bineq1;bineq2;bineq3;bineq4;bineq5]);
Aeq = sparse(Aeq1);
beq = sparse(beq1);

%% Formulate objective
f0 = [0*cost'*Uf,zeros(1,Tf),lamini*ones(1,Tini)];
f1 = [0*cost'*Uf,lamY*ones(1,Tf),zeros(1,Tini)];
f2 = [lamU*cost'*Uf,zeros(1,Tf),zeros(1,Tini)];
H = sparse(diag([lamg*ones(1,num_g),zeros(1,Tf),0*lamini*ones(1,Tini)]));

optionsL = optimoptions('linprog','Display','none');
options=optimoptions('quadprog','LinearSolver','sparse','Display','off');
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

    
% Solution
y_sp_full = Yf(1:Tf,:)*gopt(1:num_g);
y_sp = y_sp_full(1);

ufull = Uf(1:Tf,:)*gopt(1:num_g);
u = ufull(1);


end