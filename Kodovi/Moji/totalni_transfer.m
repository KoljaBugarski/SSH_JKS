clear all
close all
clc

%%
z_pocetak=1e-5;
z_kraj=4.5e-5;
dz=1e-7;
n_z=((z_kraj-z_pocetak)/dz)+1;
z=linspace(z_pocetak,z_kraj,n_z);

A = 25.911e2;
B = 0.1921e6;
coupling_coeff=A*exp(-B*z);

L_pocetak=0;
L_kraj=10;
dL=1e-3;
n_L=((L_kraj-L_pocetak)/dL)+1;
L=linspace(L_pocetak,L_kraj,n_L);
%%
for i=1:n_z

    [L_simulirano(i),vek]=f_transfer11_z(coupling_coeff(i),L);
    maxx(i)=abs(vek(max(size(vek)),2))^2;

end
%%

L_teorija=(pi/2)./coupling_coeff;
figure
plot(z,L_simulirano)

