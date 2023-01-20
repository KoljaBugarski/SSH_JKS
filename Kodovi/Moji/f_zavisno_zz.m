function vek_pravi=f_zavisno_zz(t,poc_uslov,H,kubna_nl, saturaciona_nl, gama)
n_t=max(size(t));
t_pocetak=t(1);
dt=t(3)-t(2);
for i_t=1:n_t
    tt=[t_pocetak+(i_t-1)*dt i_t*dt];
    options = odeset('RelTol',1e-9,'AbsTol',1e-9);
    [ttt,vek_t]=ode45(@f_zavisno_z, tt, poc_uslov, options, H(:,:,i_t), kubna_nl, saturaciona_nl, gama);  
    poc_uslov=vek_t(max(size(ttt)),:);
    if i_t~=1
        vek_pravi(i_t,:)=vek_t(max(size(ttt)),:);
    else
        vek_pravi(i_t,:)=poc_uslov;
    end

end
end