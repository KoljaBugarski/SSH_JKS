function [La,izlaz,fun,H,kk,sp]=f_WGA_SVI_OSIM_SREDNJEG(n,psi_out)
% Racuna duzine odgovarajucih delova strukture kako bi se dobio zeljeni
% kapler
od_koliko=0.2;
do_koliko=15;
korak=0.1;
broj=uint32((do_koliko-od_koliko)/korak+1);
psi_in=zeros(1,n)';
psi_in(fix(n/2)+1)=1;
H=zeros(n,n,2);
L=zeros(fix(n/2),broj);
for i=1:fix(n/2)
    L(i,:)=(od_koliko:korak:do_koliko)*pi;
end
for i=1:n
    for j=1:n
        if abs(i-j)==1
            H(i,j,1)=1;
            H(j,i,1)=1;
        end
    end
end
a=zeros(n,n);
a((fix(n/2)+1),fix(n/2))=1;
a(fix(n/2),(fix(n/2)+1))=1;
a((fix(n/2)+1),(fix(n/2)+2))=1;
a((fix(n/2)+2),(fix(n/2)+1))=1;
H(:,:,2)=H(:,:,1)-a;
% for i =1:2
%     E(i)=eig(H(:,:,i));
% end

min=10;
for j=1:broj
    for k=1:broj
        psi_outA=expm(-1i*H(:,:,2)*L(2,j))*expm(-1i*H(:,:,1)*L(1,k))*psi_in;
        if min > sum(abs(abs(psi_out).^2-abs(psi_outA).^2)) && abs(psi_outA((fix(n/2)+1)))^2<0.001
            min=sum(abs(abs(psi_out).^2-abs(psi_outA).^2));
            L1=L(1,k);
            L2=L(2,j);
            kk=abs(psi_outA((fix(n/2)+1)))^2;
        end
    end
end
mat2=expm(-1i*H(:,:,2));
mat1=expm(-1i*H(:,:,1));
E2=mat2(fix(n/2)+1,:);
E1=mat1(fix(n/2)+1,:);
%sred=abs(expm(-1i*H(fix(n/2)+1,:,2)*x(2))*expm(-1i*H(fix(n/2)+1,:,1)*x(1))*psi_in).^2;
fun =@(x) sum(abs(abs(expm(-1i*H(:,:,2)*x(2))*expm(-1i*H(:,:,1)*x(1))*psi_in).^2-abs(psi_out).^2));%+abs(E1*exp(x(1))*psi_in)^2;
Lovi=[L1 L2];
x = fminsearch(fun,Lovi);
La=x;
izlaz=expm(-1i*H(:,:,2)*La(2))*expm(-1i*H(:,:,1)*La(1))*psi_in;
P_target=abs(psi_out).^2;
P_simulated=abs(izlaz).^2;
sp=(sum(P_target-P_simulated).^2)/n;
end