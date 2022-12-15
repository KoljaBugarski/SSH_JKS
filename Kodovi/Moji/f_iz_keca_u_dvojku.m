function [L]=f_iz_keca_u_dvojku(v)
% Racuna duzinu talasovoda za koje dodje do prelaska iz sredisnjeg
% talasovoda u dva bocna pri cemu je odnos snaga u bocnim talasovodima
% 50/50, za odgovarajucu vrednost koeficijenta sprezanja v
od_koliko=0.2;
do_koliko=2;
korak=0.1;
broj=(do_koliko-od_koliko)/korak+1;
H=zeros(3,3);
H(1,2)=v;
H(2,1)=v;
H(2,3)=v;
H(3,2)=v;
psi_in=zeros(1,3)';
psi_in(2)=1;
psi_out=(1/sqrt(2))*[1 0 1]';
L=(od_koliko:korak:do_koliko)*pi;
min=10;
for k=1:broj
    psi_outA=expm(-1i*H*L(k))*psi_in;
    if min > sum(abs(abs(psi_out).^2-abs(psi_outA).^2))
        min=sum(abs(abs(psi_out).^2-abs(psi_outA).^2));
        L1=L(k);
    end
end
fun =@(x) sum(abs(abs(expm(-1i*H*x(1))*psi_in).^2-abs(psi_out).^2));
L=[L1];
x = fminsearch(fun,L);
L=x;
end


