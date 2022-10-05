function dadt = zavisno_z(t,a,H,sp_mod)
dadt=(-1i*H+diag((-1i)*(abs(a)).^2)*sp_mod)*a;
end