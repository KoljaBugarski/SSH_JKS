clear all
close all
clc
%%
v1=1; % koeficijent sprege za 3 talasovoda
n=11; % broj talasovoda 
v2=1; % koeficijent sprege za n talasovoda

t_pocetak=0;
t_kraj=50;
dt=0.001;
n_t=uint32((t_kraj-t_pocetak)/dt+1);
t=linspace(t_pocetak,t_kraj,n_t);   % vremenski/duzinski opseg

t1=f_transfer12(v1,t); % vreme/duzina potrebna da se desi 50/50 za 3 talasovoda
[t2,vek_1x3]=f_transfer12(v2,t); % vreme/duzina potrebna da se desi 50/50 za 3 talasovoda ali ovde posle nastavljamo propagaciju
[t3,vek_1x1]=f_transfer11(v2,t); % vreme/duzina potrebna da se u potpunosti transferuje snaga iz jednog talasovada u drugi
[t4,vek_1x2]=f_transfer1ujenda_polovina(v2,t); % vreme/duzina potrebna da se u pola snage iz jendog transferuje u drugi talasovod

% provera za 3 talasovoda
psi_inneisp=[0 1 0]';
H_neisp=[0 v1 0; v1 0 v1; 0 v1 0];
psi_outneisp=expm(-1i*H_neisp*t1)*psi_inneisp;
figure 
bar(1:3,abs(psi_outneisp).^2);

% provera za n talasovoda
tt=zeros(1,(((n-3)/2)+1));
sredina=fix(n/2)+1;
tt(1)=t2;
for i=2:(((n-3)/2)+1)
    tt(i)=t3;
end
psi_inisp=zeros(1,n);
for i=1:5
    psi_inisp(sredina)=1;
end
psi_inisp=psi_inisp';
H=zeros(n,n,uint32((n-1)/2));
H(sredina,sredina-1,1)=v2;
H(sredina,sredina+1,1)=v2;
H(sredina-1,sredina,1)=v2;
H(sredina+1,sredina,1)=v2;
for ii=2:fix(n/2)
    dole_desno=sredina+(ii-1);
    gore_levo=sredina-(ii-1);
    H(gore_levo,gore_levo-1,ii)=v2;
    H(gore_levo-1,gore_levo,ii)=v2;
    H(dole_desno,dole_desno+1,ii)=v2;
    H(dole_desno+1,dole_desno,ii)=v2;
end
izlaz=psi_inisp;
for i=1:(((n-3)/2)+1)
    izlaz=expm(-1i*H(:,:,i)*tt(i))*izlaz;
end
figure 
bar(1:n,abs(izlaz).^2);

% ukupna vremena/duzine 
T_n_talasovoda=t2+((n-3)/2)*t3;
T_tri_talasovoda=t1;

% rabi oscilacije =1
rabi=t3*v2/(pi/2);

x_1x3=t(1:max(size(vek_1x3)));
y_1x3=vek_1x3(:,2);
x_1x1=t(1:max(size(vek_1x1)));
y_1x1=vek_1x1(:,1);
x_1x2=t(1:max(size(vek_1x2)));
y_1x2=vek_1x2(:,1);               % iz fitovanja se vidi da se snaga u talasovodima u koje smo ubacili svetlost menja kao cos^2
                            % stim sto je u slucaju 1x2 omega=v(koeficijent sprege), a 1x3 w=sqrt(2)*v 


