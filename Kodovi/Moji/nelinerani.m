function dadt = nelinerani(t,a,H,on_site,sp_mod,ran)
dadt=(-1i*H+diag((-1i)*(abs(a)).^2)*sp_mod+diag(ran)*(-1i)*on_site)*a;
end