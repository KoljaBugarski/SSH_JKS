clear all
close all
%% LSA ANALYSIS for SSH
N=32;

%v=0.6;
w=1;

Kx0=0;

deltamax=0.1*2;

NLsteps=51;


deltalist = linspace(deltamax,-deltamax,NLsteps);
deltaplist = linspace(0,1,NLsteps);
stability=zeros(NLsteps);
largest=stability;
type=stability;

ksize = 2*pi/N;
klist = ksize * [(0:N-1)-N/2];
Kx =klist;

%Delta=deltalist;
[Delta,Deltap]=meshgrid(deltalist,deltaplist);

for p=1:NLsteps
  for q=1:NLsteps
     delta=Delta(p,q);
     v=Deltap(p,q);

[vecs,betas]= eig(ssh(Kx0,v,w));
betas=diag(betas);
%  if D > 2*sqrt(2) %|| D >= 0
%   beta0=betas(2);
%   vec0=vecs(:,2);    
%   else
beta0=betas(1);
vec0=vecs(:,1);
%   end

     
HLSA = zeros(4);

energy = beta0 + delta*(abs(vec0(1)).^2+abs(vec0(2)).^2); %/(1+(abs(vec0(1)).^2+abs(vec0(2)).^2));
f1 = delta*abs(vec0(1)).^2;  %/(1+abs(vec0(1)).^2);    % \delta f(I0)   for first vector
fp1 = f1;  %f1*deltap ;  %\delta f(I0) *deltap, where deltap=df/dI*I/f
f2 = delta*abs(vec0(2)).^2;  %/(1+abs(vec0(2)).^2);
fp2 = f2;   %*deltap;

betaLSA=zeros(N,4);
maxgrowth=zeros(N,1);

for n=1:N
  
        kx=Kx(n);
      
    
        HLSA(1:2,1:2) = ssh(Kx0+kx,v,w);
        HLSA(3:4,3:4) = -conj(ssh(Kx0-kx,v,w));
        
        HLSA(1,1) = HLSA(1,1) - energy + f1 + fp1;
        HLSA(2,2) = HLSA(2,2) - energy + f2 + fp2;
        
        HLSA(3,3) = HLSA(3,3) + energy - f1 - fp1;
        HLSA(4,4) = HLSA(4,4) + energy - f2 - fp2;       
        
        HLSA(1,3) = fp1;
        HLSA(2,4) = fp2;
        HLSA(3,1) = -fp1;
        HLSA(4,2) = -fp2;                     
        
        betaLSA(n,:) = eig(HLSA);
        maxgrowth(n) = max(imag(betaLSA(n,:)));
        
        if kx==0 
            maxgrowth(n)=nan; %exclude goldstone mode
        end

  end


    
stability(p,q) = min(maxgrowth(:));

[fastest,index]=max(maxgrowth(:));


largest(p,q) = fastest;

   end

   end
figure
pcolor(Deltap,Delta,stability./(1e-6+stability))
shading flat; axis square;
title('Does instability span whole BZ?')
% measured by the smallest nonzero imaginary part of LSA eigenvalue, excluding
% Goldstone mode
% 1 = yes, 0 = no
colorbar
xlabel('v/w ')
ylabel('Nonlinearity')

figure
pcolor(Deltap,Delta,largest./(1e-6+abs(largest)))
shading flat; axis square;
title('Instability type')
% -1 = complex instability, 1 = real instability, 0 = stable
% this classification assumes the particle-hole symmetry, i.e. LSA
% eigenvalues must occur in +- and complex conjugate pairs, so complex
% instability means all four eigenvalues are complex
caxis([-1 1]); colorbar
ylabel('Nonlinear frequency shift')
xlabel('v/w')

figure
pcolor(Deltap,Delta,abs(largest))
shading flat; axis square;
title('Largest instability growth rate')
% largest imaginary part of LSA eigenvalue
%caxis([0 2]); 
colorbar
ylabel('Nonlinear frequency shift')
xlabel('v/w')

% figure
% pcolor(Deltap,Delta,abs(klargestx))
% shading flat; axis square;
% title('k - largest growth rate')
% % largest imaginary part of LSA eigenvalue
% %caxis([0 1]); 
% colorbar
% ylabel('Nonlinear frequency shift')
% xlabel('Detunning parameter \Delta')
% 
% figure
% pcolor(Deltap,Delta,klargesty)
% shading flat; axis square;
% title('ky - largest growth rate')
% % largest imaginary part of LSA eigenvalue
% %caxis([0 1]); 
% colorbar
% ylabel('Nonlinear frequency shift')
% xlabel('Detunning parameter \Delta')
% %% 
% %% 



figure
subplot(1,2,1)
pcolor(Deltap,Delta,abs(largest))
shading flat; axis square;
title('Instability rate')
% largest imaginary part of LSA eigenvalue
%caxis([0 2]); 
colorbar
ylabel('\Gamma I_0/J_1')
xlabel('v/w')
subplot(1,2,2)
pcolor(Deltap,Delta,largest./(1e-6+abs(largest)))
shading flat; axis square;
title('Instability type')
caxis([-1 1]); colorbar
ylabel('\Gamma I_0/J_1')
xlabel(' v/w')




