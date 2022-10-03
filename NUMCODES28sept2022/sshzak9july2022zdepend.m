clear all
close all


%Probing bulk topology by leaky modes and modulation instability (modification of the Daniel's code)
% Leaky modes --> extension to bulk+termal bath (i.e. main WGs + auxilliary waveguides
% MI - only main WGs are considered


N_array = 64*4; % number of main waveguides
N_env = 150; % number of auxilliary waveguides per radiation channel

N_env=200;

N_ensemble=1; % number of disorder realizations
N_waveguides = N_array+2*N_env;

% coupling constants (in units of /um)
lda0 = 0.8; % [um]
t0=2*pi/lda0*0.025; % *0.025
t2=0.35*t0; % t2>t1 trivial % t1=0.5*t0;
t1=1*t0; % t2=1*t0;

t1=0.202855997;
t2=0.065495291; 
t1=0.6;
% t1=0.065495291; 
% t2=0.6;
omega_0 = 0;


%% Environment effective params 
% 
% t_env=0.35*t0; 
% t_eps=0.45*t0; % between SSH and cladding 
% m_env=-0.8*t0;
% 
gamma=t0*0.0025*0; % loss 
t_env = 0.416 * t0; % calculated from dimer splitting in the auxiliary array: 
% t_env = 0.0104 *(2*pi/lda0);
m_env= -0.7578 * t0 ; % calculated from detuning:
% m_env = (1.7444 - data1wg(7,3))*(2*pi/lda0); 
t_eps = 0.44 * t0; % estimated from dimer splitting g_s: 
% t_eps =  0.011*(2*pi/lda0); 



m_env = -0.1770;
t_env = 0.108465956;
t_eps=0.1394551;

t_env = 0.076848494;


%%
% disorder 
W=0.1*t0*0;

% propagation distance (in um)
L = 2000;

Npoints=1000;

zlist = linspace(0,L,Npoints);
tsteps=Npoints;

gamma=0.0;
gg=0.05;
com0=zeros(N_ensemble,tsteps);
com1=com0;
Zak_phase=com0;
P_array0=zeros(N_ensemble,tsteps);
P_array1=P_array0;
P_array2=P_array0;
H0 = zeros(N_waveguides);

%% SSH Hamiltonian 
H_NL = zeros(N_array,N_array)+ omega_0*eye(N_array,N_array);
for j=1:N_array-1
    
    if mod(j,2) == 0
        H_NL(j,j+1) = t1;
        H_NL(j+1,j) = t1;
    else
        H_NL(j,j+1) = t2;
        H_NL(j+1,j) = t2;
    end
    
end

%% auxiliary array Hamiltonian 

H_env = zeros(N_env,N_env);

H_env  = kron(eye(N_env),(m_env-1i*gamma)) + kron(diag(ones(1,(N_env-1)),1), t_env) + kron(diag(ones(1,(N_env-1)),-1), transpose(conj(t_env)));

H0=[H_env,zeros(N_env,N_array), zeros(N_env,N_env);
    zeros(N_array,N_env), H_NL, zeros(N_array,N_env);
	zeros(N_env,N_env),zeros(N_env,N_array), H_env];

H0(N_env,N_env+1) = t_eps; H0(N_env+1,N_env) = t_eps;
H0(N_env+N_array,N_env+N_array+1) = t_eps; H0(N_env+N_array+1,N_env+N_array) = t_eps;

%% add disorder (in SSH main array)
Hdisorder=zeros(N_waveguides);

for nd=1:N_ensemble

disorder=W*(rand(N_array,1)-0.5);
    for j=N_env+1:N_env+N_array-1
        Hdisorder(j,j+1)=disorder(j-N_env);
        Hdisorder(j+1,j)=disorder(j-N_env);
    end

% Full Hamiltonian = main WGs+ disorder, thermal bath
H = H0 + Hdisorder;
%% solve evolution 

% initial condition - index0 - main WGs+disorder; 1 - +thermal bath; 2-
% main- NL
psi0 = zeros(N_waveguides,1);
psi1=zeros(N_array,1);
psi2=zeros(N_array,1);

input_pos= N_array/2;
psi0(N_env + N_array/2) = 1;
psi1(input_pos) = 1;

% dynamics with auxiliary waveguides
[~,y0] = ode45( @(t,y) -1i*H*y, zlist, psi0);

% dynamics without auxiliary waveguides + linear lattice
[~,y1] = ode45( @(t,y) -1i*H_NL*y, zlist, psi1);

% initial condition for NL lattice
psi2(input_pos) = 1;
x_ul=psi2;
x_ul = x_ul + 1e-3.*rand(N_array,1).*exp(1i*2*pi*rand(N_array,1));
disorder_strength=0;
       
 options = odeset('RelTol',1e-8,'AbsTol',1e-8);
  % dynamics without auxiliary waveguides + nonlinearity                      
 [~,y2]=ode45(@(t,y) -1i.*H_NL*y-1i.*gg.*(abs(y).^2).*y,zlist,x_ul,options); % Evolution Runge-Kutta

%% take psi from the main array only 
psi0=y0(:,N_env+1:N_env+N_array);  % leaky effect
psi1=y1(:,1:N_array);   % diosrder
psi2=y2(:,1:N_array);   % NL

% norm in all cases
P_array0(nd,:) = sum(abs(psi0).^2,2);
P_array1(nd,:) = sum(abs(psi1).^2,2);
P_array2(nd,:) = sum(abs(psi2).^2,2);

X=transpose([1:N_array]);

%projector=zeros(tsteps,N_array/2,2,2);

% compute observables

for j=1:tsteps
    com0(nd,j) = sum((X(:)-input_pos).*squeeze(transpose(abs(psi0(j,:)).^2)))./P_array0(nd,j);
    com1(nd,j) = sum((X(:)-input_pos).*squeeze(transpose(abs(psi1(j,:)).^2)))./P_array1(nd,j);    
    com2(nd,j) = sum((X(:)-input_pos).*squeeze(transpose(abs(psi2(j,:)).^2)))./P_array2(nd,j); 
    
    psiA0 = psi0(j,1:2:end-1);
    psiB0 = psi0(j,2:2:end);

    psiA0ft=fft(psiA0);
    psiB0ft=fft(psiB0);
    
    
    psiA1 = psi1(j,1:2:end-1);
    psiB1 = psi1(j,2:2:end);

    psiA1ft=fft(psiA1);
    psiB1ft=fft(psiB1);

    psiA2 = psi2(j,1:2:end-1);
    psiB2 = psi2(j,2:2:end);

    psiA2ft=fft(psiA2);
    psiB2ft=fft(psiB2);
    
    
 % don't need to calculate these observables as they are not easily measurable in the experiment    
   projectorl(j,:,1,1) = abs(psiA0ft).^2;
   projectorl(j,:,2,2) = abs(psiB0ft).^2;
   projectorl(j,:,1,2) = conj(psiB0ft).*psiA0ft;
   projectorl(j,:,2,1) = conj(psiA0ft).*psiB0ft;   

   connection=eye(2);
   
   for n=1:N_array/2
   connection = squeeze(projectorl(j,n,:,:))*connection;    
   end
% leaky effect
   Zak_phasel(nd,j) = angle(trace(connection))/pi; 
   Zak_phase_2l(nd,j) = angle(trace(connection));


   projectorln(j,:,1,1) = abs(psiA1ft).^2;
   projectorln(j,:,2,2) = abs(psiB1ft).^2;
   projectorln(j,:,1,2) = conj(psiB1ft).*psiA1ft;
   projectorln(j,:,2,1) = conj(psiA1ft).*psiB1ft;   
    
   connection=eye(2);
   
   for n=1:N_array/2
   connection = squeeze(projectorln(j,n,:,:))*connection;    
   end
% disordered lattice - linear

   Zak_phaseln(nd,j) = angle(trace(connection))/pi; 
   Zak_phase_2ln(nd,j) = angle(trace(connection));

   projector(j,:,1,1) = abs(psiA2ft).^2;
   projector(j,:,2,2) = abs(psiB2ft).^2;
   projector(j,:,1,2) = conj(psiB2ft).*psiA2ft;
   projector(j,:,2,1) = conj(psiA2ft).*psiB2ft;   
    
   connection=eye(2);
   
   for n=1:N_array/2
   connection = squeeze(projector(j,n,:,:))*connection;    
   end
% NL lattice
   Zak_phase(nd,j) = angle(trace(connection))/pi; 
   Zak_phase_2(nd,j) = angle(trace(connection));
    
end

end

% imagesc(abs(H)); colorbar
%%
figure
plot(zlist,mean(P_array0,1),'b')
hold on
plot(zlist,mean(P_array1,1),'r')
hold on
plot(zlist,mean(P_array2,1),'g')
plot(zlist,0.5*ones(tsteps,1),'--k')
% interval=99;
% errorbar(zlist(1:interval:end),mean(P_array0(:,1:interval:end),1),std(P_array0(:,1:interval:end),1),'b','LineStyle','none')
xlim([0 L])
%ylim([0 1.05])
axis square
box on
plot([100 100],[0 1],'--k')
xlabel('Propagation distance (\mum)')
title('Total wavepacket norm N')


figure
% pcolor(abs(y0))
 imagesc(abs(y0))


figure
if (t1<t2)
scatter(zlist,Zak_phase_2) % correct 
else 
scatter(zlist,Zak_phase_2) % correct 
%plot(zlist,wrapTo2Pi(Zak_phase_2))
end
hold on 
xlabel('Propagation distance (\mum)')
ylabel('Zak phase \nu (MI induced)')
%ylim([-2 2])
axis square
box on


figure
if (t1<t2)
scatter(zlist,Zak_phase_2l) % correct 
else 
scatter(zlist,Zak_phase_2l) % correct 
%plot(zlist,wrapTo2Pi(Zak_phase_2))
end
hold on 
xlabel('Propagation distance (\mum)')
ylabel('Zak phase L\nu (leaky effect)')
%ylim([-2 2])
axis square
box on


figure
if (t1<t2)
scatter(zlist,Zak_phase_2ln) % correct 
else 
scatter(zlist,Zak_phase_2ln) % correct 
%plot(zlist,wrapTo2Pi(Zak_phase_2))
end
hold on 
xlabel('Propagation distance (\mum)')
ylabel('Zak phaseLN \nu (disordered linear lattice')
%ylim([-2 2])
axis square
box on




figure
plot(zlist,mean(com2,1),'b') 
hold on
% plot(zlist,mean(com1,1),'r') % SSH
plot(zlist,0.5*ones(tsteps,1),'--k')
% plot(zlist,mean(com0-com1,1),'k') 
%errorbar(zlist(1:interval:end),mean(com0(:,1:interval:end),1),std(com0(:,1:interval:end),1),'b','LineStyle','none')
hold off
xlim([0 L])
xlabel('Propagation distance (\mum)')
title('Mean displacement <\Deltax>')
axis square
box on

