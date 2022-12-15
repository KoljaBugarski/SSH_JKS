clear all;
close all;
clc;
%% ISCRTAVA IZLAZ IZ STUKTURE U ZAVISNOSTI OD FAZNE RAZLIKE U KRAJNJIM TALASOVODIMA
%      n=5             n=7
%       |              |
%      | |            | |
%     |   |          |   |
%      | |          |     | <- u ovim talasovodima se menja faza
%       |            |   |
%                     | |
%                      |
%                                    

n=9; % broj talasovoda
v=1.2; % koeficijent sprege

L12=iz_keca_u_dvojku(v);
Ltotal=totalni_transfer(v);
k=(1/2)*n-(3/2);
LL=zeros(1,k);
if n==3
    L=[L12 1 L12];
else
    for i =1:k
        LL(i)=Ltotal;
    end
    L=[L12 LL 1 LL L12];
end

nn=21;
promena_faze=linspace(-pi,pi,nn);
izlaz=zeros(n,nn);
for iii=1:nn
[psi_outINF,H_INF,faza]=f_interferometar(n,L,promena_faze(iii),v);
izlaz(fix(n/2)+1,iii)=abs(psi_outINF(fix(n/2)+1))^2;
if nn<22
    figure;
    bar(1:n,abs(izlaz(:,iii)).^2);
    ylim([0 1]);
end
end
figure
plot(promena_faze/pi,izlaz(fix(n/2)+1,:));
xlabel('Faza/\pi')
ylabel('Intenzitet u srednjem talasovodu')