% This function arranges "training" data in appropriate Hankel form and updates "initialisation" data at each sample

function [ZoneControlData] = BuildDataMatrices(k,Meas,Ctrlparams)

uin = Meas.uin;
y = Meas.y;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Zone training data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
utr = uin(2:min(k,1+Ctrlparams.T),:);
ytr = y(2:min(k,1+Ctrlparams.T),:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Zone initialise data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if k-Ctrlparams.Tini>0    
    uini = uin(k-Ctrlparams.Tini+1:k,:);
    yini = y(k-Ctrlparams.Tini+1:k,:);
else
    uini = zeros(1,numZ);
    yini = zeros(1,numZ);
end


ZoneControlData.utr = utr;
ZoneControlData.ytr = ytr;
ZoneControlData.uini = uini;
ZoneControlData.yini = yini;