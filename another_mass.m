%% Double mass spring simulator
function [t,y,up,ws] = another_mass(tspan,y0,Tr,m1,m2,k1,k2,d1,d2,usp,urand,w0,trlen,trgap)

    function [dydt,u_p,w] = odefun(t,y)
        
        function fval = f(t,gap,ran)
            fval = ran(min(trlen/gap,(ceil(t/gap))));
        end

        A = [0 0 1 0;0 0 0 1;-(k1+k2)/m1 k2/m1 -(d1+d2)/m1 d2/m1;k2/m2 -k2/m2 d2/m2 -d2/m2];
        B = [0;0;1/m1;0];

        if Tr==1
            u_p = f(t,trgap,urand);
            w = f(t,1,w0);

        else
            u_p = usp;
            w = w0;
        end
        
        dydt = A*y+B*u_p+B*w;

    end

[t,y]=ode45(@odefun,tspan,y0);
up = 0*t;
ws = 0*t;
for ii = 1:length(t)
    [~,up(ii),ws(ii)]=odefun(t(ii),y(ii,:)');
end

end