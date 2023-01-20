clear all;
close all;
clc;

%%

n=11; % broj talasovoda
t_pocetak=0;
t_kraj=500;
dt=0.01;
n_t=((t_kraj-t_pocetak)/dt)+1;
t=linspace(t_pocetak,t_kraj,n_t);


w=0.1+0.8*cos(0.00314*t);
v=0.1+0.8*sin(0.00314*t);
figure
plot(t,w,'b');
hold on;
plot(t,v,'r');
legend('w','v');

%%
H=zeros(n,n,n_t);
for i_t=1:n_t
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

E_sve=zeros(n_t,n);
for i_t=1:n_t
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
kubna_nl=0;
saturaciona_nl=0;
gama=1;
ran=zeros(1,max(size(vek_E)));
on_site=0;
vek_pravi=f_evolucija_zavisno_z(t,poc_uslov,H,kubna_nl, saturaciona_nl, gama,ran,on_site);



X=linspace(1,n,n);
figure;
mesh(X,t,abs(vek_pravi).^2);
colorbar;
title('Propagacija u vremnu');
xlabel('Cvorovi');
ylabel('Vreme');
xlim([1 (n)]);
ylim([0 t_kraj]);

vek_tt=zeros(n_t,2*(n)+1);
for tt=1:1:n_t
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