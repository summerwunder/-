function [Up,Uf,Yp,Yf] = getHankels(utr,ytr,Ctrlparams)

function [hp,hf] = buildhank(tr,T,Tf,Tini,n)
    hp = zeros(Tini*n,T-Tf-Tini+1);
    hf = zeros((Tini+Tf-Tini)*n,T-Tf-Tini+1);

    for jj = 1:n
        pcol = tr(1:Tini,jj);
        prow = tr(Tini:T-Tf,jj);
        fcol = tr(Tini+1:Tini+Tf,jj);
        frow = tr(Tini+Tf:T,jj);
        pj = hankel(pcol,prow);
        fj = hankel(fcol,frow);
        for kk = 1:length(pj(:,1))
            hp((kk-1)*n+jj,:) = pj(kk,:);
        end
        for kk = 1:length(fj(:,1))
            hf((kk-1)*n+jj,:) = fj(kk,:);
        end
    end

end

Tini = Ctrlparams.Tini;
T = Ctrlparams.T;
Tf = Ctrlparams.Tf;

num_g = T-Tini-Tf+1;

Up = zeros(Tini,num_g);
Uf = zeros(Tf,num_g);
Yp = zeros(Tini,num_g);
Yf = zeros(Tf,num_g);

[Up(1:Tini,1:(T-Tini-Tf+1)),Uf(1:Tf,1:(T-Tini-Tf+1))] = buildhank(utr,T,Tf,Tini,1);
[Yp(1:Tini,1:(T-Tini-Tf+1)),Yf(1:Tf,1:(T-Tini-Tf+1))] = buildhank(ytr,T,Tf,Tini,1);


end