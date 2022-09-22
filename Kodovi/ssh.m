function blochH = ssh(kx,v,w)
blochH=zeros(2);
blochH(1,1) =0;
blochH(2,2) = blochH(1,1);
blochH(1,2) =v+w*exp(-1i*kx);    
blochH(2,1) = v+w*exp(1i*kx);
end

