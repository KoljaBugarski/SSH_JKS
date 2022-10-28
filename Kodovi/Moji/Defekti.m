clear all;
%close all;
clc;

%% parametri

n=50; % koliko dimera bez umetnutih, racunajuci dva kranja
w=0.5;
v =1.2;
for i =1:10
    defekti(i)=uint32(rand*100);
end
defekti=sort(defekti);
na_koja_mesta=[0]; % mora 0 na kraju da stoji, defekat ce u novoj strukturi biti talasovod sa tim rednim brojem
koliko_defekata=max(size(na_koja_mesta))-1;
srednji_defekat=int8((koliko_defekata+1)/2);
                                                     %pr n=98 [20 39 52 65 84]
                                                                  
% Pr. sa jednim defektom, n=10 defekt na 6. mestu:
% Ako je w>v:
% O   OO   OO   O   OO   OO   O
% 1   23   45   6   78   91   1
%                         0   1
%
% Ako je w<v:
% OO   OO   OOO   OO   OO
% 12   34   567   89   11
%                      01
%
% Pr. sa vise defekata, n=14 defekt na 6. i 11. mestu:
% Ako je w>v:
% O   OO   OO   O   OO   OO   O   OO   OO   O
% 1   23   45   6   78   91   1   11   11   1 
%                         0   1   23   45   6
%
% Ako je w<v:
% OO   OO   OOO   OO   OOO   OO   OO
% 12   34   567   89   111   11   11
%                      012   34   56

t_pocetak=0;
t_kraj=500;
dt=0.01;
t_br_tacaka=(t_kraj-t_pocetak)/dt;
t=linspace(t_pocetak,t_kraj,t_br_tacaka);

%nelinearnost
on_site=0; %onsite disorder (mora i gh=1)
gama=5; % za onsite disorder
kubna_nl=1; % za kubnu
saturaciona_nl=0; %saturacionu

%1-ima hopping disorder ili 0-nema
hd=0;
jacina_hd=1e-2;

% perturbacija na ulazni vekotr
per=0;
per_na_poc_vek=0.05; % koliko jaka perturbacija

tri_d=1; % da li da crrta mesh ili ne
tir_d_cvor_predstavljen_sa_dve_tacke_na_x_osi=1;% lepse se vidi ravan xy

pocetni_uslov=4; % 0-u srednji defekat sva snaga, 1-neki od svojstevih vekota,2-lin superpozicija edge moda, 3-sta god 
svo_polje_u_koji_talasovod=n;
koji_sv_vektor=n;
koliko_vek_u_sp_poz=3;

%thermal bath
n_okoline_levo=5;
n_okoline_desno=5;
veza_okoline_i_strukture=0.1;
v_env=0.8;
w_env=0.2;
N=2*(n_okoline_levo+n_okoline_desno+n)+koliko_defekata;
uporedi_kada_nema_thermal_bath=1;


%% formiranje Hamiltonijana H

Hn=zeros(2*n+koliko_defekata,2*n+koliko_defekata); 
i=1;
ii=1;
bul=0;
while i<2*n+koliko_defekata                  
    for j=1:1:i                             % Gledamo samo donju trougaonu matricu jer je simetricna u odnosnu na gl dijagonalu
        if abs(i-j)==1
                if na_koja_mesta(ii) ~= i   % Ako nema defekta unosi normalno s tim sto unosenjem defekta menja uslov da li na tom mestu treba da bude
                    if bul                  % v ili w. Kada se ubaci defekat, uslovi se zamene sto je resno preko promenljive bul i operacije xor
                        if mod(i+j,4)==3
                            Hn(i,j)=w;
                            Hn(j,i)=w;
                        end
                        if mod(i+j,4)==1
                            Hn(i,j)=v;
                            Hn(j,i)=v;
                        end
                    else
                        if mod(i+j,4)==3
                            Hn(i,j)=v;
                            Hn(j,i)=v;
                        end
                        if mod(i+j,4)==1
                            Hn(i,j)=w;
                            Hn(j,i)=w;
                        end
                    end
                else
                    Hn(i,j)=v;              % Kada postoji defekat moram na to mesto da ubacim v i na sledece mesto u matrici opet zbog cega se i povecava ovde
                    Hn(j,i)=v;
                    i=i+1;
                    Hn(i,i-1)=v;
                    Hn(i-1,i)=v;
                    ii=ii+1;
                    bul=bitxor(bul,1);
                end
        end
    end
    i=i+1;
end 
Hn(i,i-1)=v;
Hn(i-1,i)=v;

%% Izgled strukture

ii=1;
figure
if w>v                                     % ako je w>v  dimeri su "pocepani" pa zapravo ono sto na slici izgleda kao dimer, to su cvorovi razlicitih 
    for i=1:1:2*n+koliko_defekata          % jedinicnih celija. Na slici su blizi cvorovi cija je medjusobna veza jaca
        for j=1:1:i
            if i-j==1
                if Hn(i,j)==w
                    y1(ii)=0;
                    jj=ii;
                    ii=ii+1;
                    line([jj,ii],[0,0],'color','blue')
                    hold on
                end
                if Hn(i,j)==v
                    if i==1
                        y1(ii)=0;
                        jj=ii;
                        ii=ii+1;
                        y1(ii)=1;
                        ii=ii+1;
                        figure
                        line([jj,ii],[0,0],'color','red')
                        hold on;
                    else
                        y1(ii)=0;
                        jj=ii;
                        ii=ii+1;
                        y1(ii)=1;
                        ii=ii+1;
                        line([jj,ii],[0,0],'color','red')
                        hold on;
                    end
                end
            end
        end
    end
else
    for i=1:1:2*n+koliko_defekata          
        for j=1:1:i
            if i-j==1
                if Hn(i,j)==v
                    y1(ii)=0;
                    jj=ii;
                    ii=ii+1;
                    line([jj,ii],[0,0],'color','red')
                    hold on
                end
                if Hn(i,j)==w
                    if i==1
                        y1(ii)=0;
                        jj=ii;
                        ii=ii+1;
                        y1(ii)=1;
                        ii=ii+1;
                        figure
                        line([jj,ii],[0,0],'color','blue')
                        hold on;
                    else
                        y1(ii)=0;
                        jj=ii;
                        ii=ii+1;
                        y1(ii)=1;
                        ii=ii+1;
                        line([jj,ii],[0,0],'color','blue')
                        hold on;
                    end
                end
            end
        end
    end
end
y1(ii)=0;
scatter(1:ii,y1);
hold on;
ylim([-0.5 0.5]);
xlim([1 ii]);
title(['Izgled strukture: w=' num2str(w), ' v=' num2str(v), ' .w-plavo, v-crveno']);

% %Izgled strukture
% y1=zeros(1,2*n+1);
% y1(1)=0;
% bit=0;
% for i=1:1:2*n+1
%     for j=1:1:2*n+1
%         if j-i==1
%             if Hn(i,j) == v
%                 bit=bitxor(bit,1);
%                 if bit
%                     y1(i+1)=1;
%                 else
%                     y1(i+1)=0;
%                 end
%             end
%             if Hn(i,j) == w
%                 y1(i+1)=y1(i);
%             end
%         end
%     end
% end
% figure
% scatter(1:2*n+1,y1);
% ylim ([-0.5 1.5]);

%% Svojstvene vrednosti i vekori Hamiltonijana H

[vek_E,d_E]=eig(Hn); 
E=zeros((2*n+koliko_defekata),1);
i=0;
for j=1:1:2*n+koliko_defekata
    E(j)=d_E(j,j);              % koje su zero mode, uslov se menja u zavinsoti od
    if abs(E(j))<0.2            % v i w jer se povecava odnosno smanjuje procep
        i=i+1;
        zero_mode(i)=j;
    end
end

% Crtanje svojstvenih energija
figure;
xx=1:(2*n+koliko_defekata);
scatter(1:(2*n+koliko_defekata),E);
if w > v
    title('Energije svojstvenih stanja, w>v')
else
    title('Energije svojstvenih stanja w<v')
end
% 
% axes('position',[.65 .175 .25 .25])
% box on % put box around new pair of axes
% indexOfInterest = (xx < zero_mode(max(size(zero_mode)))+2) & (xx > zero_mode(1)-2);
% scatter(xx(indexOfInterest),E(indexOfInterest));
% ylim([-0.5 0.5])
% xlim([zero_mode(1)-1 zero_mode(max(size(zero_mode)))+1])
%  

% Crtanje zero moda
if w>v
    vek_zero=zeros(2*n+koliko_defekata,koliko_defekata+2);
    for i=1:koliko_defekata+2
        vek_zero(:,i)=vek_E(:,n+i-1);
        figure;
        bar(1:(2*n+koliko_defekata),abs((vek_E(:,n+i-1))).^2);
        ylim([0 1]);
        a=koliko_defekata-1;
        if (a)
            if w>v
                title('Zero moda. Defekati na vise mesta. w>v');
            else
                title('Zero moda. Defekati na vise mesta. w<v');
            end
        else
            if w>v
                title(['Zero moda. Defekat na ' num2str(na_koja_mesta(1)), '. mestu. w>v']);
            else
                title(['Zero moda. Defekat na ' num2str(na_koja_mesta(1)), '. mestu. w<v']);
            end
        end
    end
else
    vek_zero=zeros(2*n+koliko_defekata,koliko_defekata);
    for i=1:koliko_defekata
        vek_zero(:,i)=vek_E(:,n+i-1);
        figure;
        bar(1:(2*n+koliko_defekata),abs((vek_E(:,n+i-1))).^2);
        ylim([0 1]);
        a=koliko_defekata-1;
        if (a)
            if w>v
                title('Zero moda. Defekati na vise mesta. w>v');
            else
                title('Zero moda. Defekati na vise mesta. w<v');
            end
        else
            if w>v
                title(['Zero moda. Defekat na ' num2str(na_koja_mesta(1)), '. mestu. w>v']);
            else
                title(['Zero moda. Defekat na ' num2str(na_koja_mesta(1)), '. mestu. w<v']);
            end
        end
    end
end

%% Evolucija

% pocetni uslovi 
poc_uslov=zeros(2*n+koliko_defekata,1);
switch pocetni_uslov
    case 0
        poc_uslov(svo_polje_u_koji_talasovod)=1;
    case 1
        poc_uslov=vek_E(:,koji_sv_vektor);
    case 2
        fi=0;
        poc_uslov(1)=1;
        poc_uslov(max(size(poc_uslov)))=1*exp((1i)*fi);
        poc_uslov=1/(sqrt(2))*poc_uslov;
    case 3
        poc_uslov=vek_E(:,n);
    case 4
        poc_uslov=(1/sqrt(2*n+koliko_defekata))*ones(2*n+koliko_defekata,1);
end

% perturbacija na ulazni vektor
if per
    for j=1:1:2*n+koliko_defekata
        poc_uslov(j)=poc_uslov(j)+(rand*2-1)*per_na_poc_vek;
    end
    norm=sqrt(poc_uslov'*poc_uslov);
    poc_uslov=(1/norm)*poc_uslov;
end

% hoping disorder
if hd
    for j=1:1:2*n+koliko_defekata
        for l=1:1:2*n+koliko_defekata
          if(l==j+1) || (j==l+1)
              switch(o)
                  case 0
                      Hn(j,l)=Hn(j,l)+jacina_hd*(rand*2-1);
                      o=o+1;
                  case 1
                      Hn(j,l)=Hn(j,l)+jacina_hd*(rand*2-1);
                      o=o+1;
                  case 2
                      Hn(j,l)=Hn(j,l)+jacina_hd*(rand*2-1);
                      o=o+1;
                  case 3 
                      Hn(j,l)=Hn(j,l)+jacina_hd*(rand*2-1);
                      o=0;
              end
          end
        end
    end
end

% onsite disorder
ran=ones(2*n+koliko_defekata,1);
for j=1:2*n+koliko_defekata
    ran(j)=gama*ran(j)*(rand*2-1);
end

options = odeset('RelTol',1e-9,'AbsTol',1e-9);
[ttt,vek_t]=ode45(@nelinerani, t, poc_uslov, options, Hn, on_site, kubna_nl, saturaciona_nl, gama, ran); 
vek_pravi=vek_t;

power=zeros(t_br_tacaka,1);
fidelity=zeros(t_br_tacaka,1);
IP=zeros(t_br_tacaka,1);
Pd=zeros(t_br_tacaka,n);
for i_t=1:t_br_tacaka
    power(i_t)=sum(abs(vek_pravi(i_t,:).^2));
    fidelity(i_t)=abs(conj(vek_pravi(i_t,:))*transpose(vek_pravi(1,:)));
    IP(i_t)=power(i_t)/(sum(abs(vek_pravi(i_t,:).^4)));
    if koliko_defekata==0
        for i=1:2:2*n-1
            Pd(i_t,i)=((i+1)/2)*(abs(vek_pravi(i_t,i))^2+abs(vek_pravi(i_t,i+1))^2);
        end
    end
end
Pdd=zeros(t_br_tacaka,1);
for i=1:t_br_tacaka
    Pdd(i)=sum(Pd(i,:));
end
W=2*(1/t(t_br_tacaka))*trapz(t,Pdd)


% Sta je uslo a sta je izaslo
figure;                
bar(1:2*n+koliko_defekata,abs(vek_pravi(1,:)).^2);    
title('Vektor na ulazu');
xlabel('Cvor');
ylim([0 1])
figure;                
bar(1:2*n+koliko_defekata,abs(vek_pravi(t_br_tacaka,:)).^2);    
title('Vektor na izlazu');
xlabel('Cvor');
ylim([0 1])

%propagacija u vremenu 3d
if tri_d
    X=linspace(1,2*n+koliko_defekata,2*n+koliko_defekata);
    figure;
    mesh(X,t,abs(vek_pravi).^2);
    colorbar;
    title('Propagacija u vremnu');
    xlabel('Cvorovi');
    ylabel('Vreme');
    xlim([1 (2*n+koliko_defekata)]);
    ylim([0 t_kraj]);
end

% Jedan cvor 2 tacke na x osi
if tir_d_cvor_predstavljen_sa_dve_tacke_na_x_osi
    if n<6
        vek_tt=zeros(t_br_tacaka,2*(2*n+koliko_defekata)+1);
        for tt=1:1:t_br_tacaka
            for i=1:1:2*n+koliko_defekata
                vek_tt(tt,i*2-1)=vek_pravi(tt,i);
                vek_tt(tt,i*2)=vek_pravi(tt,i);
            end
        end
        vek_tt(:,2*(2*n+koliko_defekata)+1)=vek_pravi(:,2*n+koliko_defekata);
        X=linspace(1,2*(2*n+koliko_defekata)+1,2*(2*n+koliko_defekata)+1);
        figure;
        mesh(X,t,abs(vek_tt).^2);
        colorbar;
        title('Propagacija u vremnu');
        xlabel('Cvorovi');
        ylabel('Vreme');
        xlim([1 2*(2*n+koliko_defekata)+1]);
        ylim([0 t_kraj]);
    end
end

% Racunanje snage na ivicama
maxx=0;
P=zeros(t_br_tacaka,2*n+koliko_defekata);
for tt=1:t_br_tacaka
    for j=1:2*n+koliko_defekata
        P(tt,j)=conj(vek_pravi(tt,j))*transpose(vek_pravi(tt,j));
    end
    pp=P(tt,1) + P(tt,2*n+koliko_defekata);
    if (pp>maxx)
        maxx=pp;
        vremeski_trenutak_tacka_snaga_max_na_ivicama=tt;
        vremeski_trenutak_snaga_max_na_ivicama=t(tt);
    end
end

% crtanje vekotra sa najvecom snagaom na ivicama
figure;                
bar(1:2*n+koliko_defekata,abs(vek_pravi(vremeski_trenutak_tacka_snaga_max_na_ivicama,:)).^2);    
title(['Najveca snaga na ivicama. Vreme: ', num2str(vremeski_trenutak_snaga_max_na_ivicama)]);
xlabel('Cvor');
ylim([0 1])

%% thermal bath

% H_env=zeros(N,N);
% for i=2:2*n_okoline_levo
%     if mod(i,2)==0
%         H_env(i,i-1)=v_env;
%         H_env(i-1,i)=v_env;
%     else
%         H_env(i,i-1)=w_env;
%         H_env(i-1,i)=w_env;
%     end
% end
% H_env(2*n_okoline_levo+1,2*n_okoline_levo)=veza_okoline_i_strukture;
% H_env(2*n_okoline_levo,2*n_okoline_levo+1)=veza_okoline_i_strukture;
% 
% for i=2:2*n+koliko_defekata
%     H_env(2*n_okoline_levo+i,2*n_okoline_levo+i-1)=Hn(i,i-1);
%     H_env(2*n_okoline_levo+i-1,2*n_okoline_levo+i)=Hn(i-1,i);
% end
% 
% H_env(2*n_okoline_levo+2*n+koliko_defekata+1,2*n_okoline_levo+2*n+koliko_defekata)=veza_okoline_i_strukture;
% H_env(2*n_okoline_levo+2*n+koliko_defekata,2*n_okoline_levo+2*n+koliko_defekata+1)=veza_okoline_i_strukture;
% for i=2*n_okoline_levo+2*n+koliko_defekata+2:N
%     if mod(n,2)==0
%         if mod(i,2)==0
%             H_env(i,i-1)=v_env;
%             H_env(i-1,i)=v_env;
%         else
%             H_env(i,i-1)=w_env;
%             H_env(i-1,i)=w_env;
%         end
%     else
%         if mod(i,2)==0
%             H_env(i,i-1)=w_env;
%             H_env(i-1,i)=w_env;
%         else
%             H_env(i,i-1)=v_env;
%             H_env(i-1,i)=v_env;
%         end
%     end
% end
% 
% [vek_E_env,d_E_env]=eig(H_env); 
% E_env=zeros((N),1);
% i=0;
% for j=1:1:N
%     E_env(j)=d_E_env(j,j);           
%     if abs(E_env(j))<0.2           
%         i=i+1;
%         zero_mode_env(i)=j;
%     end
% end
% 
% % Crtanje svojstvenih energija sa thermal bathom
% figure;
% scatter(1:N,E_env);
% title('Energije svojstvenih stanja sa thermal bathom')
% xlim([1 N]);
% 
% % Zero mode sa thermal bathom
% for i=1:max(size(zero_mode))
%     if uporedi_kada_nema_thermal_bath==1
%         vek_pom=zeros(1,N);
%         vek_pom(2*n_okoline_levo+1:2*n_okoline_levo+2*n+koliko_defekata)=vek_E(:,zero_mode(i));
%         figure;
%         bar(1:N,abs(vek_pom).^2);
%         ylim([0 1]);
%         title(['Bez thermal bata ' num2str(i), '. zero moda']);
%         figure;
%         bar(1:N,abs((vek_E_env(:,zero_mode_env(i)))).^2);
%         ylim([0 1]);
%         title([num2str(i), '. zero moda sa thermal bathom. Br. talasovoda levo:' num2str(2*n_okoline_levo), '; Br. talasovoda desno:' num2str(2*n_okoline_desno)]);
%     else
%         figure;
%         bar(1:N,abs((vek_E_env(:,zero_mode_env(i)))).^2);
%         ylim([0 1]);
%         a=koliko_defekata-1;
%         title([num2str(i), '. zero moda sa thermal bathom. Br. talasovoda levo:' num2str(2*n_okoline_levo), '; Br. talasovoda desno:' num2str(2*n_okoline_desno)]);
%     end
% end
% 
% %% Evolucija sa thermal bathom
% 
% % pocetni uslovi 
% poc_uslov=zeros(N,1);
% switch pocetni_uslov
%     case 0
%         poc_uslov(2*n_okoline_levo+na_koja_mesta(srednji_defekat))=1;
%     case 1
%         poc_uslov=vek_E_env(:,koji_sv_vektor);
%     case 2
%         for i=1:koliko_vek_u_sp_poz
%             poc_uslov=(1/(sqrt(koliko_vek_u_sp_poz)))*vek_E_env(:,n+i-1)+poc_uslov;
%         end
%     case 3
%         poc_uslov=vek_E_env(:,n);
% end
% 
% ran=ones(N,1);
% options = odeset('RelTol',1e-9,'AbsTol',1e-9);
% [t,vek_t_env]=ode45(@nelinerani, t, poc_uslov, options, H_env, on_site, snaga, kubna_nl, gh, ran);
% 
% % Sta je uslo a sta je izaslo
% figure;                
% bar(1:N,abs(vek_t_env(1,:)).^2);    
% title('Vektor na ulazu, thermal bath');
% xlabel('Cvor');
% ylim([0 1])
% figure;                
% bar(1:N,abs(vek_t_env(t_br_tacaka,:)).^2);    
% title('Vektor na izlazu, thermal bath');
% xlabel('Cvor');
% ylim([0 1])
% 
% %propagacija u vremenu 3d
% if tri_d
%     X=linspace(1,N,N);
%     figure;
%     mesh(X,t,abs(vek_t_env).^2);
%     colorbar;
%     title('Propagacija u vremnu, thermal bath');
%     xlabel('Cvorovi');
%     ylabel('Vreme');
%     xlim([1 (N)]);
%     ylim([0 t_kraj]);
% end