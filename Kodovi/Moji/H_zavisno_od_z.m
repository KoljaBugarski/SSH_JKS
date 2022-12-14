clear all;
close all;
clc;

%%

n=11; % broj talasovoda
t_pocetak=0;
t_kraj=500;
dt=0.01;
t_br_tacaka=(t_kraj-t_pocetak)/dt;
t=linspace(t_pocetak,t_kraj,t_br_tacaka);

kubna_nl=0;
saturaciona_nl=0;
gama=1;

w=0.1+0.8*cos(0.00314*t);
v=0.1+0.8*sin(0.00314*t);
figure
plot(t,w,'b');
hold on;
plot(t,v,'r');
legend('w','v');

%%
H=zeros(n,n,t_br_tacaka);
for i_t=1:t_br_tacaka
    for i=1:n
        for j=1:n
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
E=zeros((n),1);
i=0;
for j=1:1:n
    E(j)=d_E(j,j);             
end

E_sve=zeros(t_br_tacaka,n);
for i_t=1:t_br_tacaka
    E_sve(i_t,:)=eig(H(:,:,i_t));
end

figure
for i=1:n
    plot(t,E_sve(:,i));
    hold on;
end

%%

% poc_uslov=zeros(n,1);  % Riplovi se dobijaju jer ovo nije svojsveni
% poc_uslov(1)=1;        % vektor pocetnog H, ali ovaj ispod jeste
poc_uslov=vek_E(:,6);
vek_pravi=zeros(t_br_tacaka,n);
vek_pravi(1,:)=poc_uslov;

power=zeros(t_br_tacaka,1);
fidelity=zeros(t_br_tacaka,1);
IP=zeros(t_br_tacaka,1);

for i_t=1:t_br_tacaka
    tt=[t_pocetak+(i_t-1)*dt i_t*dt];
    options = odeset('RelTol',1e-12,'AbsTol',1e-12);
    [ttt,vek_t]=ode45(@f_zavisno_z, tt, poc_uslov, options, H(:,:,i_t), kubna_nl, saturaciona_nl, gama);  
    poc_uslov=vek_t(max(size(ttt)),:);
    if i_t~=1
        vek_pravi(i_t,:)=vek_t(max(size(ttt)),:);
    end
    power(i_t)=sum(abs(vek_pravi(i_t,:).^2));
    fidelity(i_t)=abs(conj(vek_pravi(i_t,:))*transpose(vek_pravi(1,:)));
    IP(i_t)=power(i_t)/(sum(abs(vek_pravi(i_t,:).^4)));

end



X=linspace(1,n,n);
figure;
mesh(X,t,abs(vek_pravi).^2);
colorbar;
title('Propagacija u vremnu');
xlabel('Cvorovi');
ylabel('Vreme');
xlim([1 (n)]);
ylim([0 t_kraj]);

vek_tt=zeros(t_br_tacaka,2*(n)+1);
for tt=1:1:t_br_tacaka
    for i=1:1:n
        vek_tt(tt,i*2-1)=vek_pravi(tt,i);
        vek_tt(tt,i*2)=vek_pravi(tt,i);
    end
end
vek_tt(:,2*(n)+1)=vek_pravi(:,n);
X=linspace(1,2*(n)+1,2*(n)+1);
figure;
mesh(X,t,abs(vek_tt).^2);
colorbar;
title('Propagacija u vremnu');
xlabel('Cvorovi');
ylabel('Vreme');
xlim([1 2*(n)+1]);
ylim([0 t_kraj]);