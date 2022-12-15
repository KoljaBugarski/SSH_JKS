function dadt = f_zavisno_z(t,a,H,kubna,saturaciona,gama)
dadt=(-1i*H+diag((-1i)*gama*(abs(a)).^2)*kubna+diag((-1i)*gama*((abs(a)).^2)./(1+(abs(a)).^2))*saturaciona)*a;
end