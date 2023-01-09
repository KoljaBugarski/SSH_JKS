clear all
close all
clc
%%
v1=1; % koeficijent sprege za 3 talasovoda
n=11; % broj talasovoda 
v2=2; % koeficijent sprege za n talasovoda

t1=f_transfer12(v1); % vreme/duzina potrebna da se desi 50/50 za 3 talasovoda
t2=f_transfer12(v2); % vreme/duzina potrebna da se desi 50/50 za 3 talasovoda ali ovde posle nastavljamo propagaciju
t3=f_transfer11(v2); % vreme/duzina potrebna da se u potpunosti transferuje snaga iz jednog talasovada u drugi

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