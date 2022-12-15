function [L]=f_totalni_transfer(v)
% Racuna duzinu talasovoda za koje dodje do maksimalnog transfera iz jednog
% talasovoda u drugi kada postoji sprega samo izmedju njih,
% za odgovarajucu vrednost koeficijenta sprezanja v
od_koliko=0.2;
do_koliko=2;
korak=0.1;
broj=(do_koliko-od_koliko)/korak+1;
H=zeros(2,2);
H(1,2)=v;
H(2,1)=v;
psi_in=zeros(1,2)';
psi_in(1)=1;
psi_out=[0 1]';
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


