clear all;
close all; 
clc;

%% CRTA IZLAZ IZ DVE RAZLICITE STUKTURE ZA KAPLER ZA ZADATI IZLAZ samo za n=3,5,7,9
% "Kontinualni" n=7   "Isprekidani"(samo 50/50) n=7     "Svi osim sred" n=7 
%         |                            |                      |
%        |||                          | |               | | |   | | | 
%       |||||                        |   |
%      |||||||                      |     |

n=7; % broj talasovoda
psi_out=[1 0 0 0 0 0 1]'; % zeljeni izlaz
psi_out=(1/sqrt(sum(abs(psi_out).^2)))*psi_out;
[L_KON,psi_outKON,H_KON,SP_KON,psi_in]=f_WGA_kontinualan(n,psi_out);
figure;
bar(1:n,abs(psi_outKON).^2);
ylim([0 1])
title('Bez korekcije, kontinulani')

funKON =@(x) sum(abs(abs(expm(-1i*H_KON(:,:,3)*x(3))*expm(-1i*H_KON(:,:,2)*x(2))*expm(-1i*H_KON(:,:,1)*x(1))*psi_in).^2-abs(psi_out).^2)); % expm(-1i*H_KON(:,:,5)*x(5))*expm(-1i*H_KON(:,:,4)*x(4))*
x = fminsearch(funKON,L_KON);
izlazKON_korekcija=expm(-1i*H_KON(:,:,3)*x(3))*expm(-1i*H_KON(:,:,2)*x(2))*expm(-1i*H_KON(:,:,1)*x(1))*psi_in; % expm(-1i*H_KON(:,:,5)*x(5))*expm(-1i*H_KON(:,:,4)*x(4))*
figure;
bar(1:n,abs(izlazKON_korekcija).^2);
ylim([0 1])
title('Sa korekcijom, kontinulani')
P_target=abs(psi_out).^2;
P_simulated=abs(izlazKON_korekcija).^2;
SP_KON_korekcija=sum((P_target-P_simulated).^2)/n;

if psi_out(1)==1/sqrt(2)
    if n~=3
        [L_ISP,psi_outISP,H_ISP,SP_ISP]=f_WGA_isprekidan(n,psi_out);
        figure;
        bar(1:n,abs(psi_outISP).^2);
        ylim([0 1]);
        title('Bez korekcije, isprekidani')

        funISP =@(x) sum(abs(abs(expm(-1i*H_ISP(:,:,3)*x(3))*expm(-1i*H_ISP(:,:,2)*x(2))*expm(-1i*H_ISP(:,:,1)*x(1))*psi_in).^2-abs(psi_out).^2)); % expm(-1i*H_ISP(:,:,5)*x(5))*expm(-1i*H_ISP(:,:,4)*x(4))*
        x = fminsearch(funISP,L_ISP);
        izlazISP_korekcija=expm(-1i*H_ISP(:,:,3)*x(3))*expm(-1i*H_ISP(:,:,2)*x(2))*expm(-1i*H_ISP(:,:,1)*x(1))*psi_in; % expm(-1i*H_ISP(:,:,5)*x(5))*expm(-1i*H_ISP(:,:,4)*x(4))*
        figure;
        bar(1:n,abs(izlazISP_korekcija).^2);
        ylim([0 1])
        title('Sa korekcijom, isprekidani')
        P_target=abs(psi_out).^2;
        P_simulated=abs(izlazISP_korekcija).^2;
        SP_ISP_korekcija=sum((P_target-P_simulated).^2)/n;
    end
end

if (psi_out(fix(n/2)+1))==0
    [L_SVI_OSIM_SREDNJEG,psi_outSVI_OSIM_SREDNJEG,H_SVI_OSIM_SREDNJEG,kk,SP_SVI_OSIM_SREDJNEG,psi_in]=f_WGA_SVI_OSIM_SREDNJEG(n,psi_out);
    figure;
    bar(1:n,abs(psi_outSVI_OSIM_SREDNJEG).^2);
    ylim([0 1]);

end

% n=5 za razlicite korake, kontinualna struktura
% 5.8815 5.7125 korak 0.1
% 0.7829 1.5427 korak 0.01
% 3.660032597135474   2.908965863430914 korak 0.001 [1 1 1 1 1]
