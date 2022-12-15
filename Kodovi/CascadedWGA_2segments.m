clear all
close all

M = 5; 
L=2*pi;
%% Wanted output state
psic=[1 1 1 1 1]';
psic=psic/sqrt(sum(abs(psic).^2));
%% Input state
psi0=[zeros(1,M)]';
if mod(M,2)==0 psi0(M/2)=1; else psi0(ceil(M/2))=1; end
%% Search L1-L2 space to find which combination gives the output closest to the wanted output state
L1v=(0.2:0.1:2)*pi;
L2v=(0.2:0.1:2)*pi;
for brL1=1:length(L1v)
    for brL2=1:length(L2v)
% Construct coupling matrix of the first WGA segment
C1 = diag(ones(1,2),1) + diag(ones(1,2),-1);
C1M=padarray(C1,[(5-size(C1,1))/2 (5-size(C1,2))/2],0,'both');
T1=expm(-i*C1M*L1v(brL1));
psi1=T1*psi0;

% Construct coupling matrix of the second WGA segment
C2 = diag(ones(1,M-1),1) + diag(ones(1,M-1),-1);
C2M=padarray(C2,[(M-size(C2,1))/2 (M-size(C2,2))/2],0,'both');
T2=expm(-i*C2M*L2v(brL2));
psi2=T2*psi1;

% Fit function, here error function. Minimum fit function gives the best solution.
fitf(brL1, brL2)=sum(abs(abs(psi2).^2-abs(psic).^2));
    end
end

%% Plot fit function
contourf(L1v/pi,L2v/pi,fitf'), colorbar

%% Find fit func minimum
[fitf_min,br2Dmin]=min(fitf(:));
[brL1min,brL2min]=ind2sub(size(fitf),br2Dmin);
fitf_min, L1v(brL1min)/pi, L2v(brL2min)/pi