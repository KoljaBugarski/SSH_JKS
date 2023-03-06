function vek_pravi=f_evolucija_nova(tt,y,n,x)

H=zeros(n,n);

for ii=1:max(size(x))-1
    if tt<=x(1)
        i=0;
        break
    end
    if tt > x(ii) && tt <=x(ii+1)
        i=ii;
        break
    end
end

sredina=fix(n/2)+1;
if i==0
    if tt<x(1)
        H(sredina,sredina-1)=1;
        H(sredina,sredina+1)=1;
        H(sredina-1,sredina)=1;
        H(sredina+1,sredina)=1;
    end
end
if i~=0
    if tt > x(i) && tt <= x(i+1)
        H(sredina,sredina-1)=1;
        H(sredina,sredina+1)=1;
        H(sredina-1,sredina)=1;
        H(sredina+1,sredina)=1;
        iii=i;
        while iii~=0
        dole_desno=sredina+iii;
        gore_levo=sredina-iii;
        H(gore_levo,gore_levo-1)=1;
        H(gore_levo-1,gore_levo)=1;
        H(dole_desno,dole_desno+1)=1;
        H(dole_desno+1,dole_desno)=1;
        iii=iii-1;
        end
    end
end

vek_pravi=-1i.*H*y;      

end
    

% 
% function dy = lattice_chile_gain_loss_NL(t,y,A,gamma,gain,loss,Ncell)
% zeromat = zeros(6*Ncell,6*Ncell);
% GLN = zeromat;
% 
%  for i =1:6*Ncell
%      GLN(i,i) = (gain./(1+abs(y(i)).^2)-loss)+(1i).*gamma.*(abs(y(i)).^2);
%  end
% 
%  dy=1i.*A*y+GLN*y;
% 
% end