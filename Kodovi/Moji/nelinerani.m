function dadt = nelinerani(t,a,H,on_site,snaga,sp_mod,gh,ran)
P=0.8;
dadt=-1i*H*a+a.*ran*(-1i)*on_site*gh+(-1i)*((abs(a)).^2).*a*sp_mod+snaga*P;
end