function dadt = nelinerani(t,a,H,on_site,kubna,saturaciona,gama,ran)
dadt=(-1i*H+diag((-1i)*gama*(abs(a)).^2)*kubna+diag((-1i)*gama*((abs(a)).^2)./(1+(abs(a)).^2))*saturaciona+diag(ran)*(-1i)*on_site)*a;
end