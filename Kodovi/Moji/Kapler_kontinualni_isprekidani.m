clear all;
close all; 
clc;

%% CRTA IZLAZ IZ DVE RAZLICITE STUKTURE ZA KAPLER ZA ZADATI IZLAZ samo za n=3,5,7,9
% "Kontinualni" n=7   "Isprekidani"(samo 50/50) n=7     "Svi osim sred" n=7 
%         |                            |                      |
%        |||                          | |               | | |   | | | 
%       |||||                        |   |
%      |||||||                      |     |

n=9; % broj talasovoda
psi_out=[1 0 0 0 0 0 0 0 1]'; % zeljeni izlaz
psi_out=(1/sqrt(sum(abs(psi_out).^2)))*psi_out;
[L_KON,psi_outKON,funKON,H_KON,SP_KON]=f_WGA_kontinualan(n,psi_out);
figure;
bar(1:n,abs(psi_outKON).^2);
ylim([0 1])
if psi_out(1)==1/sqrt(2)
    if n~=3
        [L_ISP,psi_outISP,H_ISP,SP_ISP]=f_WGA_isprekidan(n,psi_out);
        figure;
        bar(1:n,abs(psi_outISP).^2);
        ylim([0 1]);
    end
end

if (psi_out(fix(n/2)+1))==0
    [L_SVI_OSIM_SREDNJEG,psi_outSVI_OSIM_SREDNJEG,funSVI_OSIM_SREDNJEG,H_SVI_OSIM_SREDNJEG,kk,SP_SVI_OSIM_SREDJNEG]=f_WGA_SVI_OSIM_SREDNJEG(n,psi_out);
    figure;
    bar(1:n,abs(psi_outSVI_OSIM_SREDNJEG).^2);
    ylim([0 1]);
end

% n=5 za razlicite korake, kontinualna struktura
% 5.8815 5.7125 korak 0.1
% 0.7829 1.5427 korak 0.01
% 3.660032597135474   2.908965863430914 korak 0.001 [1 1 1 1 1]
