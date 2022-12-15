function [izlaz,H,faza]=f_interferometar(n,L,promena_faze,v)
% Pravi odgovarajuci H i U za interferometar i racuna izlaz
psi_in=zeros(1,n)';
psi_in(fix(n/2)+1)=1;
sredina=fix(n/2)+1;
H=zeros(n,n,n);
H(sredina,sredina-1,1)=v;
H(sredina,sredina+1,1)=v;
H(sredina-1,sredina,1)=v;
H(sredina+1,sredina,1)=v;

for ii=2:fix(n/2)
    dole_desno=sredina+(ii-1);
    gore_levo=sredina-(ii-1);
    H(gore_levo,gore_levo-1,ii)=v;
    H(gore_levo-1,gore_levo,ii)=v;
    H(dole_desno,dole_desno+1,ii)=v;
    H(dole_desno+1,dole_desno,ii)=v;
end
H(1,1,fix(n/2)+1)=promena_faze;
H(n,n,fix(n/2)+1)=-promena_faze;

b=2;
bb=n-1;
for ii=(fix(n/2)+2):(2*fix(n/2)+1)
    H(b,b-1,ii)=v;
    H(b-1,b,ii)=v;
    H(bb,bb+1,ii)=v;
    H(bb+1,bb,ii)=v;
    b=b+1;
    bb=bb-1;
end
izlaz=psi_in;
faza=zeros(n,n);
for i=1:n
    izlaz=expm(-1i*H(:,:,i)*L(i))*izlaz;
    faza(:,i)=angle(izlaz);
end



end



