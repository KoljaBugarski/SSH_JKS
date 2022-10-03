clear all
close all
% parameters
% 

N=40; % grid size
%N=32*2;

% t1=1;
% t2=2;

v=0.3;
   
w=1;
D=0;
 % nonlinear coefficient
delta=1;

W=0; % disorder strength

intensity = 1 ; %plane wave intensity

Lz = 40; % propagation distance
L=N;

numbermax=1;

N_steps=2*5000; %number of propagation steps
dz = Lz/N_steps;
dx = L/N;
recording_steps=N_steps/200;
saved_steps = N_steps/recording_steps+1;
ksize = 2*pi/N;
coord = (0:N-1) * L/N - L/2;
X=transpose(coord);

klist = ksize * [(0:N-1)-N/2];
Kx=klist;

%Kx = fftshift(Kx);

evolution_operator = zeros(N,2,2);
h_operator=evolution_operator;
beta=zeros(N,2);
vec_A = zeros(N,2);
vec_B = vec_A;

for p=1:N
  
        
        h_operator(p,:,:) = ssh(Kx(p),v,w);        
        evolution_operator(p,:,:) = expm(-1i*dz*squeeze(h_operator(p,:,:)));        
        beta(p,:) = eig(squeeze(h_operator(p,:,:)));
        
        [Vec,betat]=eig(squeeze(h_operator(p,:,:)));
        
        vec_A(p,:)=Vec(:,1);
        vec_B(p,:)=Vec(:,2);

end


figure
plot(fftshift(Kx),fftshift((beta(:,1))),'o')
hold on 
plot(fftshift(Kx),fftshift((beta(:,2))),'o')
title('band diagram- linear case')
hold off

projector = zeros(N,2,2,saved_steps);
 
pnumber_R = zeros(numbermax,saved_steps);
pnumber_K = zeros(numbermax,saved_steps);

pnumber_Kx = zeros(numbermax,saved_steps);
pnumber_Ky = zeros(numbermax,saved_steps);

avg_pop1 = zeros(saved_steps,N);
avg_pop2 = avg_pop1;


load 'ulaz.txt';
xinput=ulaz;

for number=1:numbermax
    
number
psi_A = zeros(N,1);
psi_B = zeros(N,1);

%psi_A = sqrt(intensity)*exp(-1i*(X*pi));
psi_A=xinput(1:2:2*N-1); %vec_A(:,1);
psi_B=xinput(2:2:2*N); %vec_A(:,2);

input_A=psi_A;
input_B=psi_B;

psi_A = psi_A + 0e-2.*rand(N,1).*exp(1i*2*pi*rand(N,1));
psi_B = psi_B + 0e-2.*rand(N,1).*exp(1i*2*pi*rand(N,1));

total_power = sum(abs(psi_A(:)).^2+abs(psi_B(:)).^2).*dx.^2;

psi_A_ft = fft(psi_A);
psi_B_ft = fft(psi_B);
fourier_power = sum((abs(psi_A_ft(:)).^2 + abs(psi_B_ft(:)).^2));
count=0;

disorderA = W*(rand(N,1)-0.5);
disorderB = W*(rand(N,1)-0.5);


  
        count=count+1;
    
                              
        %% HAMILTONIANS FOR PRINTING
        
%        H_LA = L^2*(1/N).^4.*sum(sum( conj(psi_A_ft).*( psi_A_ft.*h_operator(:,:,1,1) + psi_B_ft.*h_operator(:,:,1,2))));
%        H_LB = L^2*(1/N).^4.*sum(sum( conj(psi_B_ft).*( psi_B_ft.*h_operator(:,:,2,2) + psi_A_ft.*h_operator(:,:,2,1))));
           
%        HLn(number,count)= real(H_LA + H_LB)/(total_power);
%        HNLn(number,count)=-delta*(dx.^2 .* sum( log(1 + abs(psi_A(:)).^2) + log(1+abs(psi_B(:)).^2)) - power_A(count)-power_B(count) )/(total_power);
            
        %-------------

        %% FORIER COMPONENT OF FIELDS
         
        %% THESE QUANTITIES SHOWS THE 'CORRELATUIONS'? PROJECTORS

        FT_int = abs(psi_A_ft).^2 + abs(psi_B_ft).^2;
        
        projector(:,1,1,count) = projector(:,1,1,count) + abs(psi_A_ft).^2./FT_int/numbermax;
        projector(:,2,2,count) = projector(:,2,2,count) + abs(psi_B_ft).^2./FT_int/numbermax;        
        projector(:,1,2,count) = projector(:,1,2,count) + psi_A_ft.*conj(psi_B_ft)./FT_int/numbermax;
        projector(:,2,1,count) =  projector(:,2,1,count) + psi_B_ft.*conj(psi_A_ft)./FT_int/numbermax;    
               
        projbloch1=abs(psi_A_ft(:).*squeeze(conj(vec_A(:,1))) +psi_B_ft(:).*squeeze(conj(vec_A(:,2)))).^2./fourier_power;
        projbloch2=abs(psi_A_ft(:).*squeeze(conj(vec_B(:,1))) +psi_B_ft(:).*squeeze(conj(vec_B(:,2)))).^2./fourier_power;
        
        for in=1:N
        avg_pop1(count,in) = (avg_pop1(count,in)) + projbloch1(in,1)/numbermax;
        avg_pop2(count,in) = (avg_pop2(count,in)) + projbloch2(in,1)/numbermax;
        end
        
        %% PARTICIPATION NUMBER (REAL SPACE AND FOURIER SPACE) 
        
        pnumber_R(number,count) = 1./(sum( abs(psi_A(:)).^4 + abs(psi_B(:)).^4)./total_power.^2)/(2*N^2);     
        pnumber_K(number,count) = 1./(sum(abs(psi_A_ft(:)).^4 + abs(psi_B_ft(:)).^4)./fourier_power.^2)/(2*N^2);   
        
        
   %     marginal_Kx = sum( abs(psi_A_ft).^2 + abs(psi_B_ft).^2,1)./fourier_power;
   %     
        
   %     pnumber_Kx(number,count) = 1./sum(marginal_Kx.^2);
   %    






for j=1:N_steps
        %%%FIRST STEP OF THE SPLIT-STEP PROCEDURE -NL OPERATOR ACTS AN THE
        %%%HALF OF STEP
    phase_A = dz*(disorderA+delta*nonlinearity(abs(psi_A).^2)) / 2;
    phase_B = dz*(disorderB+delta*nonlinearity(abs(psi_B).^2)) / 2;    
    
    psi_A = psi_A.*exp(-1i*phase_A);
    psi_B = psi_B.*exp(-1i*phase_B);
    
        %% SECOND STEP- LINEAR OPERATOR ON THE WHOLE TIME ISNTANT
        
    psi_A_ft = fft(psi_A);
    psi_B_ft = fft(psi_B);
    
    psi_A = ifft( psi_A_ft.*evolution_operator(:,1,1) + psi_B_ft.*evolution_operator(:,1,2));
    psi_B = ifft( psi_A_ft.*evolution_operator(:,2,1) + psi_B_ft.*evolution_operator(:,2,2));   

    %% FINAL STEP : THE NL OPERATOR ACTS ON THE REST OF THE HALF OF THE TIME ISNATNT
    
    phase_A = dz*(disorderA+delta*nonlinearity(abs(psi_A).^2)) / 2;
    phase_B = dz*(disorderB+delta*nonlinearity(abs(psi_B).^2)) / 2;      
    
    psi_A = psi_A.*exp(-1i*phase_A);
    psi_B = psi_B.*exp(-1i*phase_B);
    
    
    %% CONDITION FOR PLOTTING 
    if mod(j,recording_steps)==0
        
        count=count+1;
    
        psi_A_ft = fft(psi_A);
        psi_B_ft = fft(psi_B);
                                
        %% HAMILTONIANS FOR PRINTING
        
%        H_LA = L^2*(1/N).^4.*sum(sum( conj(psi_A_ft).*( psi_A_ft.*h_operator(:,1,1) + psi_B_ft.*h_operator(:,1,2))));
%        H_LB = L^2*(1/N).^4.*sum(sum( conj(psi_B_ft).*( psi_B_ft.*h_operator(:,2,2) + psi_A_ft.*h_operator(:,2,1))));
           
%        HLn(number,count)= real(H_LA + H_LB)/(total_power);
%        HNLn(number,count)=-delta*(dx.^2 .* sum( log(1 + abs(psi_A(:)).^2) + log(1+abs(psi_B(:)).^2)) - power_A(count)-power_B(count) )/(total_power);
            
        %-------------

        %% FORIER COMPONENT OF FIELDS
         
        %% THESE QUANTITIES SHOWS THE 'CORRELATUIONS'? PROJECTORS

        FT_int = abs(psi_A_ft).^2 + abs(psi_B_ft).^2;
        
        projector(:,1,1,count) = projector(:,1,1,count) + abs(psi_A_ft).^2./FT_int/numbermax;
        projector(:,2,2,count) = projector(:,2,2,count) + abs(psi_B_ft).^2./FT_int/numbermax;        
        projector(:,1,2,count) = projector(:,1,2,count) + psi_A_ft.*conj(psi_B_ft)./FT_int/numbermax;
        projector(:,2,1,count) =  projector(:,2,1,count) + psi_B_ft.*conj(psi_A_ft)./FT_int/numbermax;    
               
        projbloch1=abs(psi_A_ft(:).*squeeze(conj(vec_A(:,1))) +psi_B_ft(:).*squeeze(conj(vec_A(:,2)))).^2./fourier_power;
        projbloch2=abs(psi_A_ft(:).*squeeze(conj(vec_B(:,1))) +psi_B_ft(:).*squeeze(conj(vec_B(:,2)))).^2./fourier_power;
        
        for in=1:N
        avg_pop1(count,in) = (avg_pop1(count,in)) + projbloch1(in,1)/numbermax;
        avg_pop2(count,in) = (avg_pop2(count,in)) + projbloch2(in,1)/numbermax;
        end
        
        
        %% PARTICIPATION NUMBER (REAL SPACE AND FOURIER SPACE) 
        
        pnumber_R(number,count) = 1./(sum( abs(psi_A(:)).^4 + abs(psi_B(:)).^4)./total_power.^2)/(2*N);     
        pnumber_K(number,count) = 1./(sum(abs(psi_A_ft(:)).^4 + abs(psi_B_ft(:)).^4)./fourier_power.^2)/(2*N);   
        
        marginal_Kx = sum( abs(psi_A_ft).^2 + abs(psi_B_ft).^2,1)./fourier_power;
        
        pnumber_Kx(number,count) = 1./sum(marginal_Kx.^2);
      
        
        
        
    end   
   
end   % time

end    % ensembles

%%%% compute observables

popratio_mean=zeros(count,1);
popratio_std=popratio_mean;
beta_mean = zeros(count,1);
beta_std = zeros(count,1);
mu_mean=beta_mean;
mu_std=mu_mean;
purity=zeros(N,1);
purity_gap = beta_mean;
chern_number=purity_gap;

ueff=zeros(N,2);

PR_mean = mean(pnumber_R,1);
PR_std = std(pnumber_R,0,1);

PK_mean = mean(pnumber_K,1);
PK_std = std(pnumber_K,0,1);

azimuth_field = zeros(N,count);

for j=1:count
    p1(:,1) = ( avg_pop1(j,:)); 
    p2(:,1) = ( avg_pop2(j,:));
    
    beta_eff= (p1-p2)./(2*p1.*p2.*(beta(:,2)));
    mu_eff = -((p1+p2).*(beta(:,2)))./(p1-p2);


    beta_mean(j,1) = mean(beta_eff(:));
    beta_std(j,1)=std(beta_eff(:));

    mu_mean(j) = mean(mu_eff(:));
    mu_std(j)=std(mu_eff(:));
        
    for n=1:N
       
        
    proj = squeeze(projector(n,:,:,j));
      
    neff_matrix = 2*proj - eye(2);
    
    theta = acos(neff_matrix(1,1));
    phi = angle(neff_matrix(1,2));
    
    ueff(n,1) = cos(theta/2);
    ueff(n,2) = sin(theta/2).*exp(1i*phi);  
    
    purity(n) = 2*trace(proj*proj)-1;
        
  
    end
    
   
    
    purity_gap(j) = min(purity(:));
    
  

end 




tlist = squeeze(linspace(0,Lz,saved_steps));



figure
plot(tlist, PR_mean)
xlabel('Time t')
ylabel('P_r/(2N^2)')
hold on



figure
plot(tlist, PK_mean)
xlabel('Time t')
ylabel('P_k /(2N^2)')
hold on


figure
plot(tlist,purity_gap)
xlabel('Time t')
ylabel('Purity gap')



figure
plot(tlist,chern_number)
xlabel('Time t')
ylabel('Chern number')




gap=10;
figure
plot(tlist, beta_mean,'b')
hold on
errorbar(tlist(1:gap:end), beta_mean(1:gap:end),beta_std(1:gap:end),'b','LineStyle','none')
hold off
xlabel('Time t')
ylabel('Effective inverse temperature \beta')


figure
plot(tlist, mu_mean,'b')
hold on
errorbar(tlist(1:gap:end), mu_mean(1:gap:end),mu_std(1:gap:end),'b','LineStyle','none')
hold off
xlabel('Time t')
ylabel('Effective chemical potential \mu')



figure
imagesc(FT_int)
xlabel('Time t')
ylabel('Intensity in Fourier space')

figure
imagesc(abs(psi_A).^2 + abs(psi_B).^2)
xlabel('Time t')
ylabel('Intensity in lattice space')

figure
subplot(2,1,1)
plot([1:N],(abs(input_A(:,1))),'o')
hold on 
plot([1:N],(abs(input_B(:,1))),'o')
xlabel('n')
ylabel('Mode amplitudes  (t=0)')
subplot(2,1,2)
plot([1:N],(abs(psi_A(:,1))),'o')
hold on 
plot([1:N],(abs(psi_B(:,1))),'o')
xlabel('n')
ylabel('Mode amplitudes(t=tend)')


figure
plot(real(input_A))
hold on
plot(imag(input_A))