function [t1]=f_transfer11(v)
% vreme/duzina za totalni transfer na osnovu koeficijenta sprege
gama=0;
kubna=0;
saturaciona=0; % moguce nelinearnosti

t_pocetak=0;
t_kraj=50;
dt=0.001;
n_t=uint32((t_kraj-t_pocetak)/dt+1);
t=linspace(t_pocetak,t_kraj,n_t);   % vremenski/duzinski opseg

H=[0 v;v 0];
poc_uslov=[1 0];
maxx=0;
i_t=0;
b=1;
while b
    i_t=i_t+1;
    tt=[t_pocetak+(i_t-1)*dt i_t*dt];
    options = odeset('RelTol',1e-12,'AbsTol',1e-12);
    [ttt,vek_t]=ode45(@f_zavisno_z, tt, poc_uslov, options, H, kubna, saturaciona, gama);  
    poc_uslov=vek_t(max(size(ttt)),:);
    if i_t~=1
        vek_pravi(i_t,:)=vek_t(max(size(ttt)),:);
    else
        vek_pravi(i_t,:)=poc_uslov;
    end
    if abs(vek_pravi(i_t,2)^2)>maxx
        maxx=abs(vek_pravi(i_t,2)^2);
    else
        b=0;
    end
end
t1=t(i_t);
end