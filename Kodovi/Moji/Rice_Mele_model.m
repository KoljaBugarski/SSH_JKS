clear all
close all
clc

x=[1.2 1.1 1 0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2 0.1]; % v i w
y=[1.3 1.42 1.56 1.73 1.96 2.24 2.61 3.12 3.9 5.21 7.85 15.72]; % t 
%%

a=1.56;
b=-1.003; 
sprega=0.8;        % jedno od ova dva mora biti =0. Drugi ce biti
t1=0;              % izracunat tako da dodje do max transfera

if sprega~=0
    t1=round(a*sprega^b,2);
else
    sprega=round(((t1/a)^(1/b)),2);
end

n=5;

dt=0.01;
t_pocetak=0;
t_kraj=120;

k=round(t1);
if t1>k
    t2=k+1-t1;
else
    t2=k-t1;
end

n1=t1/dt;
n2=t2/dt;
n_t=(t_kraj-t_pocetak)/dt+1;
kk=round((n_t-1)/(2*n1+2*n2));
t=linspace(t_pocetak,t_kraj,n_t);


kubna_nl=0;
saturaciona_nl=0;
gama=0;

% f_t=zeros(1,3*n_t);
% for i=1:n_t
%     if t(i)<t_kraj/8 || t(i)==t_kraj/8
%         f_t(i)=(8/t_kraj)*t(i);
%         j=i;
%     end
%     if t(i)<3*t_kraj/8 && t(i)>t_kraj/8
%         f_t(i)=f_t(j);
%         k=i;
%     end
%     if t(i)==3*t_kraj/8
%         f_t(i)=f_t(k);
%     end
%     if t(i)<t_kraj/2 && t(i)>3*t_kraj/8
%         f_t(i)=1-(8/t_kraj)*(t(i)-t_kraj*3/8);
%     end
%     if t(i)>t_kraj/2
%         f_t(i)=0;
%     end
% end
% f_t(n_t+1:2*n_t)=f_t(1:n_t);
% f_t(2*n_t+1:3*n_t)=f_t(1:n_t);
% % figure
% % plot(t,f_t(1:n_t));
% 
% 
% b=0;
% for j=1:n_t
%     if abs(t_kraj/2-t(j))<0.00001
%         b=b+1;
%         bb=j;
%     end
% end
% u_t=zeros(1,n_t);
% for i=1:n_t
%     j=n_t+i;
%     u_t(i)=0;%f_t(j)-f_t(j+bb);
% end
% 
% 
% b=0;
% for j=1:n_t
%     if abs(t_kraj/4-t(j))<0.00001
%         b=b+1;
%         bb=j;
%     end
% end
% w_t=zeros(1,n_t);
% for i=1:n_t
%     j=n_t+i;
%     w_t(i)=2*f_t(j+bb);
% end
% 
% 
% v_t=zeros(1,n_t);
% for i=1:n_t
%     j=n_t+i;
%     v_t(i)=f_t(j-bb);
% end

%%

u_t=zeros(1,n_t);
v_t=zeros(1,n_t);
w_t=zeros(1,n_t);

for i=1:1:uint32(n2)
    u_t(i)=u_t(i)+sprega;
end
for i=uint32(n2+1):1:uint32(n1+n2)
    v_t(i)=v_t(i)+sprega;
end
for i=uint32(n1+n2+1):1:uint32(n1+2*n2)
    u_t(i)=u_t(i)+sprega;
end
for i=uint32(n1+2*n2+1):1:uint32(2*n1+2*n2)
    w_t(i)=w_t(i)+sprega;
end
for j=2:kk
    v_t((j-1)*(2*n1+2*n2)+1:j*(2*n1+2*n2))=v_t(1:2*n1+2*n2);
    w_t((j-1)*(2*n1+2*n2)+1:j*(2*n1+2*n2))=w_t(1:2*n1+2*n2);
    u_t((j-1)*(2*n1+2*n2)+1:j*(2*n1+2*n2))=u_t(1:2*n1+2*n2);
    u_t((j-1)*(2*n1+2*n2)+1:j*(2*n1+2*n2))=u_t(1:2*n1+2*n2);
end

figure 
plot(t,u_t, 'g')
hold on;
plot(t,v_t, 'b')
hold on;
plot(t,w_t, 'r')
hold on;


%%
H=zeros(2*n,2*n,n_t);
for i=1:2*n
    for j=1:i
        if abs(i-j) ==1
           if (mod(i+j,4)==3)
               H(i,j,:)=v_t;
               H(j,i,:)=v_t;
           end
           if (mod(i+j,4)==1)
               H(i,j,:)=w_t;
               H(j,i,:)=w_t;
           end
        end
        if i==j
            if mod(i,2)==0
                H(i,i,:)=u_t;
            else
                H(i,i,:)=-u_t;
            end
        end
    end
end
%% Snapshots eigenvalues

E=zeros(2*n,n_t);
vek_E=zeros(2*n,2*n,n_t);
for i=1:n_t
    [vek_E(:,:,i),d]=eig(H(:,:,i));
    E(:,i)=diag(d);
end
% figure
% plot(t,E(n-1,:));
% figure
% plot(t,E(n+1+1,:));


%%
poc_uslov=zeros(2*n,1);
poc_uslov(1)=1;
% poc_uslov=vek_E(:,n,1);
vek_pravi=zeros(n_t,2*n);
vek_pravi(1,:)=poc_uslov;

power=zeros(n_t,1);
fidelity=zeros(n_t,1);
IP=zeros(n_t,1);
b=1;
bb=1;
maxx=0;
for i_t=1:n_t
    tt=[t_pocetak+(i_t-1)*dt i_t*dt];
    options = odeset('RelTol',1e-9,'AbsTol',1e-9);
    [ttt,vek_t]=ode45(@zavisno_z, tt, poc_uslov, options, H(:,:,i_t), kubna_nl, saturaciona_nl, gama);  
    poc_uslov=vek_t(max(size(ttt)),:);
    if i_t~=1
        vek_pravi(i_t,:)=vek_t(max(size(ttt)),:);
    end
    power(i_t)=sum(abs(vek_pravi(i_t,:).^2));
    fidelity(i_t)=abs(conj(vek_pravi(i_t,:))*transpose(vek_pravi(1,:)));
    IP(i_t)=power(i_t)/(sum(abs(vek_pravi(i_t,:).^4)));
    if abs(vek_pravi(i_t,2))^2>0.9
       vreme(b)=t(i_t);
       snaga(b)=abs(vek_pravi(i_t,2))^2;
       b=b+1;
    end
    if abs(vek_pravi(i_t,2*n))^2>0.99
        trenutak(bb)=t(i_t);
        snaga_u_pos(bb)=abs(vek_pravi(i_t,2*n))^2;
        bb=bb+1;
        if abs(vek_pravi(i_t,2*n))^2>maxx
            maxx=abs(vek_pravi(i_t,2*n))^2;
            trenutak_max=t(i_t);
        end
    end
end

vek_tt=zeros(n_t,2*(2*n)+1);
for tt=1:1:n_t
    for i=1:1:2*n
        vek_tt(tt,i*2-1)=vek_pravi(tt,i);
        vek_tt(tt,i*2)=vek_pravi(tt,i);
    end
end
vek_tt(:,2*(n*2)+1)=vek_pravi(:,n*2);
X=linspace(1,2*(n*2)+1,2*(n*2)+1);
figure;
mesh(X,t,abs(vek_tt).^2);
colorbar;
title('Propagacija u vremnu');
xlabel('Cvorovi');
ylabel('Vreme');
xlim([1 2*(n*2)+1]);
ylim([0 t_kraj]);