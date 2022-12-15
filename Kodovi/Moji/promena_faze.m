clear all;
close all;
clc;

%% parametri

n=20; % koliko dimera bez umetnutih, racunajuci dva kranja
v=0.5;
w =1.2;
na_koja_mesta=[0]; % mora 0 na kraju da stoji, defekat ce u novoj strukturi biti talasovod sa tim rednim brojem
koliko_defekata=max(size(na_koja_mesta))-1;
srednji_defekat=int8((koliko_defekata+1)/2);
                                                                  
t_pocetak=0;
t_kraj=50;
dt=0.01;
t_br_tacaka=(t_kraj-t_pocetak)/dt;
t=linspace(t_pocetak,t_kraj,t_br_tacaka);

%nelinearnost
on_site=0; %onsite disorder (mora i gh=1)
gama=0; % za onsite disorder
kubna_nl=0; % za kubnu
saturaciona_nl=0; %saturacionu

%1-ima hopping disorder ili 0-nema
hd=0;
jacina_hd=1e-2;

% perturbacija na ulazni vekotr
per=0;
per_na_poc_vek=0.05; % koliko jaka perturbacija

tri_d=1; % da li da crrta mesh ili ne
tir_d_cvor_predstavljen_sa_dve_tacke_na_x_osi=1;% lepse se vidi ravan xy

pocetni_uslov=1; % 0-u srednji defekat sva snaga, 1-neki od svojstevih vekota,2-lin superpozicija edge moda, 3-sta god 
svo_polje_u_koji_talasovod=n+1;
koji_sv_vektor=18;
koliko_vek_u_sp_poz=3;




%% formiranje Hamiltonijana H

Hn=zeros(2*n+koliko_defekata,2*n+koliko_defekata); 
i=1;
ii=1;
bul=0;
while i<2*n+koliko_defekata                  
    for j=1:1:i                             % Gledamo samo donju trougaonu matricu jer je simetricna u odnosnu na gl dijagonalu
        if abs(i-j)==1
                if na_koja_mesta(ii) ~= i   % Ako nema defekta unosi normalno s tim sto unosenjem defekta menja uslov da li na tom mestu treba da bude
                    if bul                  % v ili w. Kada se ubaci defekat, uslovi se zamene sto je resno preko promenljive bul i operacije xor
                        if mod(i+j,4)==3
                            Hn(i,j)=w;
                            Hn(j,i)=w;
                        end
                        if mod(i+j,4)==1
                            Hn(i,j)=v;
                            Hn(j,i)=v;
                        end
                    else
                        if mod(i+j,4)==3
                            Hn(i,j)=v;
                            Hn(j,i)=v;
                        end
                        if mod(i+j,4)==1
                            Hn(i,j)=w;
                            Hn(j,i)=w;
                        end
                    end
                else
                    Hn(i,j)=v;              % Kada postoji defekat moram na to mesto da ubacim v i na sledece mesto u matrici opet zbog cega se i povecava ovde
                    Hn(j,i)=v;
                    i=i+1;
                    Hn(i,i-1)=v;
                    Hn(i-1,i)=v;
                    ii=ii+1;
                    bul=bitxor(bul,1);
                end
        end
    end
    i=i+1;
end 
Hn(i,i-1)=v;
Hn(i-1,i)=v;


%% Svojstvene vrednosti i vekori Hamiltonijana H

[vek_E,d_E]=eig(Hn); 
E=zeros((2*n+koliko_defekata),1);
i=0;
for j=1:1:2*n+koliko_defekata
    E(j)=d_E(j,j);              % koje su zero mode, uslov se menja u zavinsoti od
    if abs(E(j))<0.2            % v i w jer se povecava odnosno smanjuje procep
        i=i+1;
        zero_mode(i)=j;
    end
end


%% unintarni operator

for j=1:2*n+koliko_defekata
    faza(j)=angle(exp((-1i)*E(j)*t_kraj));
end



U=zeros(2*n+koliko_defekata);
for j=1:2*n+koliko_defekata
    U=U+exp((-1i)*E(j)*t_kraj)*vek_E(:,j)*(vek_E(:,j))';%
end

H1=zeros(2*n+koliko_defekata);
for i=1:2*n+koliko_defekata
    H1(i,i)=pi/t_kraj;
end

[vek_E1,d_E1]=eig(H1); 
E1=zeros((2*n+koliko_defekata),1);
for j=1:1:2*n+koliko_defekata
    E1(j)=d_E1(j,j);              
end

U1=zeros(2*n+koliko_defekata);
for j=1:2*n+koliko_defekata
    U1=U1+exp((-1i)*E1(j)*t_kraj)*vek_E1(:,j)*(vek_E1(:,j))';%
end

UU=expm(-1i*Hn*t_kraj);

%% Evolucija

% pocetni uslovi 
poc_uslov=zeros(2*n+koliko_defekata,1);
switch pocetni_uslov
    case 0
        poc_uslov(svo_polje_u_koji_talasovod)=1;
    case 1
        poc_uslov=vek_E(:,koji_sv_vektor);
    case 2
        fi=0;
        poc_uslov(1)=1;
        poc_uslov(max(size(poc_uslov)))=1*exp((1i)*fi);
        poc_uslov=1/(sqrt(2))*poc_uslov;
    case 3
        poc_uslov=vek_E(:,n);
    case 4
        poc_uslov=(1/sqrt(2*n+koliko_defekata))*ones(2*n+koliko_defekata,1);
end


% perturbacija na ulazni vektor
if per
    for j=1:1:2*n+koliko_defekata
        poc_uslov(j)=poc_uslov(j)+(rand*2-1)*per_na_poc_vek;
    end
    norm=sqrt(poc_uslov'*poc_uslov);
    poc_uslov=(1/norm)*poc_uslov;
end

% hoping disorder
if hd
    for j=1:1:2*n+koliko_defekata
        for l=1:1:2*n+koliko_defekata
          if(l==j+1) || (j==l+1)
              switch(o)
                  case 0
                      Hn(j,l)=Hn(j,l)+jacina_hd*(rand*2-1);
                      o=o+1;
                  case 1
                      Hn(j,l)=Hn(j,l)+jacina_hd*(rand*2-1);
                      o=o+1;
                  case 2
                      Hn(j,l)=Hn(j,l)+jacina_hd*(rand*2-1);
                      o=o+1;
                  case 3 
                      Hn(j,l)=Hn(j,l)+jacina_hd*(rand*2-1);
                      o=0;
              end
          end
        end
    end
end

% onsite disorder
ran=ones(2*n+koliko_defekata,1);
for j=1:2*n+koliko_defekata
    ran(j)=gama*ran(j)*(rand*2-1);
end

options = odeset('RelTol',1e-9,'AbsTol',1e-9);
[ttt,vek_t]=ode45(@f_nelinearni, t, poc_uslov, options, Hn, on_site, kubna_nl, saturaciona_nl, gama, ran); 
vek_pravi=vek_t;

power=zeros(t_br_tacaka,1);
fidelity=zeros(t_br_tacaka,1);
IP=zeros(t_br_tacaka,1);
Pd=zeros(t_br_tacaka,n);
for i_t=1:t_br_tacaka
    power(i_t)=sum(abs(vek_pravi(i_t,:).^2));
    fidelity(i_t)=abs(conj(vek_pravi(i_t,:))*transpose(vek_pravi(1,:)));
    IP(i_t)=power(i_t)/(sum(abs(vek_pravi(i_t,:).^4)));
    if koliko_defekata==0
        for i=1:2:2*n-1
            Pd(i_t,i)=((i+1)/2)*(abs(vek_pravi(i_t,i))^2+abs(vek_pravi(i_t,i+1))^2);
        end
    end
end
Pdd=zeros(t_br_tacaka,1);
for i=1:t_br_tacaka
    Pdd(i)=sum(Pd(i,:));
end

% Sta je uslo a sta je izaslo
figure;                
bar(1:2*n+koliko_defekata,vek_pravi(1,:));    
title('Vektor na ulazu');
xlabel('Cvor');
ylim([-1 1])
figure;                
bar(1:2*n+koliko_defekata,vek_pravi(t_br_tacaka,:));    
title('Vektor na izlazu preko H');
xlabel('Cvor');
ylim([-1 1])
figure;                
bar(1:2*n+koliko_defekata,U1*poc_uslov);    
title('Vektor na izlazu preko U');
xlabel('Cvor');
ylim([-1 1])

figure;                
bar(1:2*n+koliko_defekata,angle(vek_pravi(1,:)));    
title('pocetna Faza');
xlabel('Cvor');
figure;                
bar(1:2*n+koliko_defekata,angle(vek_pravi(t_br_tacaka,:)));    
title('krajnja Faza preko H');
xlabel('Cvor');
figure;
bar(1:2*n+koliko_defekata,angle(U*poc_uslov));    
title('kranja Faza preko U');
xlabel('Cvor');
figure;
bar(1:2*n+koliko_defekata,angle(UU*poc_uslov));    
title('kranja Faza preko UU');
xlabel('Cvor');

%%

%U dodaje globalnu fazu na svojstveni vektor ali faza uzima vrednosti
%iintervala (-pi,pi) tako da kada predje pi onda se na -pi dodajde ostatak
%faze odnosno (fi-n*pi). Konacna forula je fi-(n+1)*fi, Jedinu ako je
%vektor edge moda ondnosno E blisko nuli, onda poludi

%

for j=1:2*n+koliko_defekata
    faza(j)=angle(exp((-1i)*E(j)*t_kraj));
end
faza_poc_vek=angle(transpose(poc_uslov));
glob_faza=-E(koji_sv_vektor)*t_kraj;
glob_faza=mod(glob_faza,2*pi);
faza_izlaz=faza_poc_vek+faza(koji_sv_vektor);
% for j=1:2*n+koliko_defekata
%     if faza_izlaz(j)>pi
%         faza_izlaz(j)=mod(faza_izlaz(j),pi);
%     end
% end

faza1=angle(U*poc_uslov)-angle(poc_uslov);
% 
% 
% for i=1:2*n
%     pom(:,i)=vek_E(:,i)*(vek_E(:,i))'*vek_E(:,n);
% end
% vek_tt=zeros(t_br_tacaka,2*(2*n)+1);
% for tt=1:1:t_br_tacaka
%     for i=1:1:2*n
%         vek_tt(tt,i*2-1)=vek_pravi(tt,i);
%         vek_tt(tt,i*2)=vek_pravi(tt,i);
%     end
% end
% vek_tt(:,2*(n*2)+1)=vek_pravi(:,n*2);
% X=linspace(1,2*(n*2)+1,2*(n*2)+1);
% figure;
% mesh(X,t,abs(vek_tt).^2);
% colorbar;
% title('Propagacija u vremnu');
% xlabel('Cvorovi');
% ylabel('Vreme');
% xlim([1 2*(n*2)+1]);
% ylim([0 t_kraj]);
% figure

