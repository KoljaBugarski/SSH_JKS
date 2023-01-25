function [L1,vek]=f_transfer11_z(coupling_coeff,L)
% vreme/duzina za totalni transfer na osnovu koeficijenta sprege
gama=0;
kubna=0;
saturaciona=0; 
on_site=0;
ran=zeros(1,2);

% moguce nelinearnosti
H=[0 coupling_coeff;coupling_coeff 0];
poc_uslov=[1 0];
maxx=0;
i_t=0;
dl=L(3)-L(2);
b=1;
while b
    i_t=i_t+1;
    tt=[L(1)+(i_t-1)*dl i_t*dl];
    options = odeset('RelTol',1e-12,'AbsTol',1e-12);
    [ttt,vek_t]=ode45(@f_evolucija, tt, poc_uslov, options, H, kubna, saturaciona, gama, ran, on_site);  
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
L1=L(i_t);
vek=abs(vek_pravi.^2);
end