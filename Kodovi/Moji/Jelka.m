_novaclear all
close all
clc
%%

tir_d_cvor_predstavljen_sa_dve_tacke_na_x_osi=1;
n=9;
x=[2.24919613457334	5.29974567278782	5.67592669104971	1.55404682865558]; % Vrednosti dobijene iz koda za racunanje duzina kaplera
n_t=10001;
t=linspace(0,sum(x),10001);

psi_in=zeros(1,n)';
psi_in(fix(n/2)+1)=1;

lpom1=zeros(1,max(size(x)));
lpom1(1)=x(1);
for i=2:max(size(x))
    lpom1(i)=sum(x(1:i));
end



options = odeset('RelTol',1e-8,'AbsTol',1e-8);
[Ti,vek_pravi]=ode45(@(tt,y)f_evolucija_nova(tt,y,n,lpom1),t,psi_in,options); % Evolution Runge-Kutta

figure
bar(n,abs(vek_pravi(n_t,:)).^2)
ylim([0 1])

% Jedan cvor 2 tacke na x osi
if tir_d_cvor_predstavljen_sa_dve_tacke_na_x_osi
    vek_tt=zeros(n_t,n);
    for tt=1:1:n_t
        for i=1:1:n
            vek_tt(tt,i*2-1)=vek_pravi(tt,i);
            vek_tt(tt,i*2)=vek_pravi(tt,i);
        end
    end
    vek_tt(:,2*n+1)=vek_pravi(:,n);
    X=linspace(1,2*n+1,2*n+1);
    figure;
    mesh(X,t,abs(vek_tt).^2);
    colorbar;
    title(['Broj talasova=',num2str(n)]);
    xlabel('Sites');
    ylabel('Time');
    xlim([1 2*n+1]);
    ylim([0 sum(x)]);
    view(2)
end


