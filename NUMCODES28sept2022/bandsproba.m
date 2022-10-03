
% Fig. 1.4 - Toplogical insulators - SSH - finite lattice

clear all
close all

N=40; % number of unit cells; each cell consists of two sites: a and b 
ksize = 2*pi/N;
klist = ksize * [(0:N-1)-N/2];
Kx=klist;

% two parameters of SSH are v and w
vend=3;
vinit=0;
vstep=0.3;
vnumb=(vend-vinit)/vstep;


    
evolution_operator = zeros(N,N);
h_operator=evolution_operator;
beta=zeros(vnumb,2*N);
vec_A = zeros(vnumb,N);
vec_B = vec_A;


w=1;
vbr=0;

    for v=vinit:vstep:vend
    vbr=vbr+1;
     vv(vbr,1)=v;      
       
   
        h_operator= sshfinite(v,w,N);     % Lattice Hamiltonian    
       
        beta(vbr,:) = eig(squeeze(h_operator(:,:)));
        
        [Vec,betat]=eig(squeeze(h_operator(:,:)));
        
        vector(vbr,:,:)=Vec(:,:);
        
        

  for k=1:N % in np.arange(-pi, pi, delta_2):
      
        H = ssh(Kx(k),v,w);
       [vecc, ev]=eig(squeeze(H(:,:)));
       evv1(:,vbr)=ev(1);
       evv2(:,vbr)=ev(2);
        evv3(:,vbr)=ev(3);
       evv4(:,vbr)=ev(4);
  end
    end


figure
plot(vv(:),(real(beta(:,1:2*N))))
title('band diagram- linear case')
 
cell=[1:2*N];
cella=transpose(cell);
veca(:,1)=vector(2,:,10);
figure
bar(cella,real(veca))

cell=[1:2*N];
cella=transpose(cell);
veca(:,1)=vector(3,:,40);
figure
bar(cella,real(veca))

