function [La,izlaz,H,sp,psi_in]=f_WGA_isprekidan(n,psi_out)
% Racuna duzine odgovarajucih delova strukture kako bi se dobio zeljeni
% kapler
od_koliko=0.2;
do_koliko=2;
korak=0.1;

pom=do_koliko+korak;
broj=(pom-od_koliko)/korak+1;
L=zeros(fix(n/2),broj);
for i=1:fix(n/2)
    L(i,:)=(od_koliko:korak:pom)*pi;
end
aa=size(L);
vrsta=aa(1);
kolona=aa(2);

psi_in=zeros(1,n)';
psi_in(fix(n/2)+1)=1;

sredina=fix(n/2)+1;
H=zeros(n,n,fix(n/2));
H(sredina,sredina-1,1)=1;
H(sredina,sredina+1,1)=1;
H(sredina-1,sredina,1)=1;
H(sredina+1,sredina,1)=1;
for ii=2:fix(n/2)
    dole_desno=sredina+(ii-1);
    gore_levo=sredina-(ii-1);
    H(gore_levo,gore_levo-1,ii)=1;
    H(gore_levo-1,gore_levo,ii)=1;
    H(dole_desno,dole_desno+1,ii)=1;
    H(dole_desno+1,dole_desno,ii)=1;
end

min=10;
i1=ones(1,vrsta);
while i1(vrsta)~=kolona
    psi_outA=psi_in;
    for i2=1:vrsta
        
        psi_outA=expm(-1i*H(:,:,i2)*L(i2,i1(i2)))*psi_outA;
        if min > sum(abs(abs(psi_out).^2-abs(psi_outA).^2))
            min=sum(abs(abs(psi_out).^2-abs(psi_outA).^2));
            for i=1:vrsta
                LL(i)=L(i,i1(i));
            end
        end
    end
    for jj=1:vrsta
        if i1(jj)==kolona && jj~=vrsta
            i1(jj+1)=i1(jj+1)+1;
            i1(jj)=1;
        end
    end
    i1(1)=i1(1)+1;
end
La=LL; 
izlaz=psi_in;
for i=1:vrsta
    izlaz=expm(-1i*H(:,:,i)*La(i))*izlaz;
end
P_target=abs(psi_out).^2;
P_simulated=abs(izlaz).^2;
sp=sum((P_target-P_simulated).^2)/n;
    
end