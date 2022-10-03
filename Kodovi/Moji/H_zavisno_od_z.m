clear all;
close all;
clc;

%%

n=5;
t_pocetak=0;
t_kraj=1000;
dt=0.01;
t_br_tacaka=(t_kraj-t_pocetak)/dt;
t=linspace(t_pocetak,t_kraj,t_br_tacaka);
w=0.4*(2+cos(0.3*t));
v=0.4*(2+sin(0.3*t-(pi/2)));

figure
plot(t,w);
hold on;
plot(t,v);
hold on;

%%
H=zeros(2*n,2*n,t_br_tacaka);
for i_t=1:t_br_tacaka
    for i=1:2*n
        for j=1:2*n
            if (abs(i-j)==1)
                if (mod(i+j,4)==3)
                    H(i,j,i_t)=v(i_t);
                    H(j,i,i_t)=v(i_t);
                end
                if (mod(i+j,4)==1)
                    H(i,j,i_t)=w(i_t);
                    H(j,i,i_t)=w(i_t);
                end
            end
        end
    end
end

[vek_E,d_E]=eig(H(:,:,1)); 
E=zeros((2*n),1);
i=0;
for j=1:1:2*n
    E(j)=d_E(j,j);             
end


%%

poc_uslov=zeros(2*n,1);
poc_uslov(5)=1;
%poc_uslov=vek_E(:,5);
vek_pravi=zeros(t_br_tacaka,2*n);
vek_pravi(1,:)=poc_uslov;

for i_t=1:t_br_tacaka
    tt=[t_pocetak+(i_t-1)*dt i_t*dt];
    options = odeset('RelTol',1e-9,'AbsTol',1e-9);
    [tt,vek_t]=ode45(@zavisno_z, tt, poc_uslov, options, H(:,:,i_t));  
    poc_uslov=vek_t(2,:);
    if i_t~=1
        vek_pravi(i_t,:)=vek_t(2,:);
    end
end

X=linspace(1,2*n,2*n);
figure;
mesh(X,t,abs(vek_pravi).^2);
colorbar;
title('Propagacija u vremnu');
xlabel('Cvorovi');
ylabel('Vreme');
xlim([1 (2*n)]);
ylim([0 t_kraj]);
