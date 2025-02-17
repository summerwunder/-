disp('Run one set of params:')

% Prediction horizon:
N = 30;
% Initialisation length:
Tini = 5;
% Regularisation weight:
LamG = 0.5;

% Sim setup
simlength_tr = 300; % number of training samples
simlength_run = 100; % number of run samples

% Model physical params
m1 = 0.5;
m2 = 1.5;
k1 = 2;
k2 = 2;
d1 = 1;
d2 = 1;
y0 = [0;0;0;0]; % initial condition

% Max/min input limits
umax = 1;
umin = -1;


dtc = 1; % control sample time
trlen = simlength_tr*dtc; % training length in seconds
runlen = simlength_run*dtc; % run length in seconds


% Build a setpoint to follow:
spch1 = 25;
spch2 = 20;
yspfull = 0.8*[0.5*ones(spch1,1);-0.3*ones(spch2,1);-0.1*ones(spch1,1);-0.5*ones(spch2,1);0.5*ones(spch1,1);-0.1*ones(spch2,1);0.5*ones(spch1,1);-0.5*ones(runlen-spch1*4-spch2*3+N*dtc,1)];
sp_band = 0; % Set point band

% Build a random input for training
rng(1)
trgap = 10;
urand=max(umin*1,min(umax*1,4*(-0.5+rand(ceil(trlen/(trgap)),1))));

% Build a random disturbance component
wrandtr = 1*(0.2*(1+sin((1/15)*(1:trlen)))'+0.3*(-0.5+rand(trlen,1)));
wrand = 1*(0.2*(1+sin((1/15)*(1:runlen/dtc)))'+0.3*(-0.5+rand(runlen/dtc,1)));

% Build an input cost
cost = ones(runlen/dtc+N,1);

% Preallocation
timey = zeros(runlen/dtc,1); % for collecting computational times
errfad = zeros(4,1); % for collecting set-point errors
timeyfad = zeros(4,1); % for collecting average computational times
yfad = zeros(4,simlength_tr+simlength_run+1); % for collecting outputs
ufad = zeros(4,simlength_tr+simlength_run+1); % for collecting inputs
wfad = zeros(4,simlength_tr+simlength_run+1); % for collecting disturbances



for includedist = [0,1]

    scenariocount = 1;

    for segmented = [0,1] % First run is unsegmeneted, second run is segmented
           
        % Controller params
        Ctrlparams.Tini = Tini;
        Ctrlparams.Tf = N;

        if segmented == 1
            Horiz = floor(Ctrlparams.Tf/Ctrlparams.Tini)*Ctrlparams.Tini;
            Ctrlparams.Tf = Ctrlparams.Tini;
        else
            Horiz = Ctrlparams.Tf;
        end

        n_ord = 10;
        Ctrlparams.T = (1+1)*(Ctrlparams.Tf+Ctrlparams.Tini+1*n_ord)-1;
        


        % Objective penalty weight (\lambda_g)
        Ctrlparams.lamg = LamG;    
        
        
        %% 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Train  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Run simulation model for training period
        disturb = includedist*wrandtr;
        [ttr,ytr,u_tr,ws] = another_mass(linspace(1,trlen+1,ceil(trlen/dtc)+1),y0,1,m1,m2,k1,k2,d1,d2,[],urand,disturb,trlen,trgap);
        
        % Gather data
        Meas.uin = u_tr(1:end);
        Meas.y = ytr(1:end,2);
        
        % Build data matrix structures
        [Data] = BuildDataMatrices(trlen/dtc,Meas,Ctrlparams);
        [Ctrlparams.Up,Ctrlparams.Uf,Ctrlparams.Yp,Ctrlparams.Yf] = getHankels(Data.utr,Data.ytr,Ctrlparams);
        
        % Get set-point forecasts
        Fcasts.ysp_lo = yspfull(1:Horiz);
        Fcasts.ysp_hi = yspfull(1:Horiz)+sp_band;
        Fcasts.cost = cost(1:Horiz);
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Run  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Call data predictive control algorithm for 1st time to get first input
        if segmented==1
            usp = CtrlSegment(Data.uini,Data.yini,Fcasts,umax,umin,u_tr(end),Ctrlparams);
        else
            usp = CtrlUnsegment(Data.uini,Data.yini,Fcasts,umax,umin,Ctrlparams);
        end
        
        t = ttr;
        yr = ytr(end,:);
        
        % Pre allocation
        trun = zeros(runlen/dtc,1); % Time vector
        yrun = zeros(runlen/dtc,4); % Output vector
        wrun = zeros(runlen/dtc,1); % Disturbance vector
        urun = zeros(runlen/dtc,1); % Input vector

        
        for k = 1:runlen/dtc

            % Update setpoints
            Fcasts.ysp_lo = yspfull(k+1:k+Horiz);
            Fcasts.ysp_hi = yspfull(k+1:k+Horiz)+sp_band;
            Fcasts.cost = cost(k+1:k+Horiz);
            w = includedist*wrand(k);

            % Call simulation model for 1 timestep
            [t,yr,u_r,wsr] = another_mass([t(end),t(end)+dtc],yr(end,:),0,m1,m2,k1,k2,d1,d2,usp,urand,w,trlen,trgap);
            
            yrun(k,:) = yr(end,:);
            trun(k) = t(end);
            urun(k) = u_r(end);
            wrun(k) = wsr(end);
            
            % Gather data
            Meas.uin = [u_tr(1:end);urun(1:k)];
            Meas.y = [ytr(1:end,2);yrun(1:k,2)];
            
            % Update initialisation matrices
            [Data] = BuildDataMatrices(round(t(end)/dtc),Meas,Ctrlparams);
                    

            % Call data predictive control algorithm
            if segmented==1
                [usp,~,~,~,ti] = CtrlSegment(Data.uini,Data.yini,Fcasts,umax,umin,u_r(end),Ctrlparams);
            else
                [usp,~,~,~,ti] = CtrlUnsegment(Data.uini,Data.yini,Fcasts,umax,umin,Ctrlparams);
            end

            % Collect computational time measurement
            timey(k,1) = ti;
            
        end
        
        % Measure and store set-point error and computational time
        err = abs(yrun(:,2)-yspfull(1:runlen/dtc));
        errfad(scenariocount+includedist*2,1)=sum(err);

        % Store output info
        yfad(scenariocount+includedist*2,:) = [ytr(:,2);yrun(:,2)]';
        ufad(scenariocount+includedist*2,:) = [u_tr;urun]';
        wfad(scenariocount+includedist*2,:) = [ws;wrun]';
        
        timeyfad(scenariocount+includedist*2,1)=mean(timey(:,1));

        

        if scenariocount==1
            disp('Unsegmented scenario complete, beginning segmented scenario')
        else
            disp('Segmented scenario complete')
        end

        scenariocount = scenariocount+1;

    end


    if includedist==0
        disp('*************************************************************')
        disp('Undisturbed scenarios complete, beginning disturbed scenarios')
        disp('*************************************************************')
    else
        disp('Disturbed scenarios complete')
    end


end

%% Plot and save
save('results/SetPointError_OneScenario.mat','errfad')
save('results/OutputVars_OneScenario.mat','yfad','ufad','wfad')

x=trun-trlen;
y1=yfad(2,simlength_tr+2:end)';
y2=yfad(4,simlength_tr+2:end)';
y3=yfad(3,simlength_tr+2:end)';
y4=yspfull(1:runlen/dtc);
xlab='Time (sec)';
ylab='Displacement (m)';
leg1='y_{2} Disturbance-free';
leg2='y_{2} Segmented formulation';
leg3='y_{2} Unsegmented formulation';
leg4='y_{sp}Set-point';

figure
figure_size = [800, 600, 600, 300];
plot(x,y1,'g-','Linewidth',2)
set(gca,'FontSize',12, 'FontName','Times New Roman')
set(gcf, 'Position', figure_size )
hold all
plot(x,y2,'-','Linewidth',2)
plot(x,y3,'r-','Linewidth',2)
plot(x,y4,'k--','Linewidth',2)
grid on
xlabel(xlab)
ylabel(ylab)
legend(leg1,leg2,leg3,leg4)
saveas(gcf, 'results/OutputsOneScenario.png')


x=[ttr;trun];
y1=wfad(3,:)';
y2=ufad(4,:)';
y3=ufad(3,:)';
xlab='Time (sec)';
ylab='Input force (N)';
leg1='Unmeasured disturbance';
leg2='Segmented formulation input';
leg3='Unsegmented formulation input';

figure
figure_size = [800, 600, 600, 300];
hold off
stairs(x,y1,'k-','Linewidth',1)
set(gca,'FontSize',12, 'FontName','Times New Roman')
set(gcf, 'Position', figure_size )
hold all
stairs(x,y2,'-','Linewidth',1)
stairs(x,y3,'r--','Linewidth',1)
grid on
xlabel(xlab)
ylabel(ylab)
legend(leg1,leg2,leg3)
saveas(gcf, 'results/InputsOneScenario.png')
