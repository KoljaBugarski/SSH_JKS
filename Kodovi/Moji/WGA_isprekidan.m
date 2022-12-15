function [La,izlaz,H]=WGA_isprekidan(n,psi_out)
% Racuna duzine odgovarajucih delova strukture kako bi se dobio zeljeni
% kapler
od_koliko=0.2;
do_koliko=2;
korak=0.1;
broj=(do_koliko-od_koliko)/korak+1;
od_koliko=0.2;
do_koliko=2;
korak=0.1;
broj=(do_koliko-od_koliko)/korak+1;
switch (n)
    case 3
        L=zeros(fix(n/2),broj);
        for i=1:fix(n/2)
            L(i,:)=(od_koliko:korak:do_koliko)*pi;
        end
        psi_in=zeros(1,n)';
        psi_in(fix(n/2)+1)=1;
        sredina=fix(n/2)+1;
        H=zeros(n,n,fix(n/2));
        H(sredina,sredina-1,1)=1;
        H(sredina,sredina+1,1)=1;
        H(sredina-1,sredina,1)=1;
        H(sredina+1,sredina,1)=1;
        for ii=2:fix(n/2)
            dole_desno=sredina+(ii-1);
            gore_levo=sredina-(ii-1);
            H(gore_levo,gore_levo-1,ii)=1;
            H(gore_levo-1,gore_levo,ii)=1;
            H(dole_desno,dole_desno+1,ii)=1;
            H(dole_desno+1,dole_desno,ii)=1;
        end
        min=10;
        for k=1:broj
            psi_outA=expm(-1i*H(:,:,1)*L(1,k))*psi_in;
            if min > sum(abs(abs(psi_out).^2-abs(psi_outA).^2))
                min=sum(abs(abs(psi_out).^2-abs(psi_outA).^2));
                L1=L(1,k);
            end
        end
        fun =@(x) sum(abs(abs(expm(-1i*H(:,:,1)*x(1))*psi_in).^2-abs(psi_out).^2));
        Lovi=[L1];
        x = fminsearch(fun,Lovi);
        La=x;
        izlaz=expm(-1i*H(:,:,1)*La(1))*psi_in;
    case 5
        L=zeros(fix(n/2),broj);
        for i=1:fix(n/2)
            L(i,:)=(od_koliko:korak:do_koliko)*pi;
        end
        psi_in=zeros(1,n)';
        psi_in(fix(n/2)+1)=1;
        sredina=fix(n/2)+1;
        H=zeros(n,n,fix(n/2));
        H(sredina,sredina-1,1)=1;
        H(sredina,sredina+1,1)=1;
        H(sredina-1,sredina,1)=1;
        H(sredina+1,sredina,1)=1;
        for ii=2:fix(n/2)
            dole_desno=sredina+(ii-1);
            gore_levo=sredina-(ii-1);
            H(gore_levo,gore_levo-1,ii)=1;
            H(gore_levo-1,gore_levo,ii)=1;
            H(dole_desno,dole_desno+1,ii)=1;
            H(dole_desno+1,dole_desno,ii)=1;
        end
        min=10;
        for j=1:broj
            for k=1:broj
                psi_outA=expm(-1i*H(:,:,2)*L(2,j))*expm(-1i*H(:,:,1)*L(1,k))*psi_in;
                if min > sum(abs(abs(psi_out).^2-abs(psi_outA).^2))
                    min=sum(abs(abs(psi_out).^2-abs(psi_outA).^2));
                    L1=L(1,k);
                    L2=L(2,j);
                end
            end
        end
        fun =@(x) sum(abs(abs(expm(-1i*H(:,:,2)*x(2))*expm(-1i*H(:,:,1)*x(1))*psi_in).^2-abs(psi_out).^2));
        Lovi=[L1 L2];
        x = fminsearch(fun,Lovi);
        La=x;
        izlaz=expm(-1i*H(:,:,2)*La(2))*expm(-1i*H(:,:,1)*La(1))*psi_in;
    case 7
        L=zeros(fix(n/2),broj);
        for i=1:fix(n/2)
            L(i,:)=(od_koliko:korak:do_koliko)*pi;
        end
        psi_in=zeros(1,n)';
        psi_in(fix(n/2)+1)=1;
        sredina=fix(n/2)+1;
        H=zeros(n,n,fix(n/2));
        H(sredina,sredina-1,1)=1;
        H(sredina,sredina+1,1)=1;
        H(sredina-1,sredina,1)=1;
        H(sredina+1,sredina,1)=1;
        for ii=2:fix(n/2)
            dole_desno=sredina+(ii-1);
            gore_levo=sredina-(ii-1);
            H(gore_levo,gore_levo-1,ii)=1;
            H(gore_levo-1,gore_levo,ii)=1;
            H(dole_desno,dole_desno+1,ii)=1;
            H(dole_desno+1,dole_desno,ii)=1;
        end
        min=10;
        for i=1:broj
            for j=1:broj
                for k=1:broj
                    psi_outA=expm(-1i*H(:,:,3)*L(3,i))*expm(-1i*H(:,:,2)*L(2,j))*expm(-1i*H(:,:,1)*L(1,k))*psi_in;
                    if min > sum(abs(abs(psi_out).^2-abs(psi_outA).^2))
                        min=sum(abs(abs(psi_out).^2-abs(psi_outA).^2));
                        L1=L(1,k);
                        L2=L(2,j);
                        L3=L(3,i);
                    end
                end
            end
        end
        fun =@(x) sum(abs(abs(expm(-1i*H(:,:,3)*x(3))*expm(-1i*H(:,:,2)*x(2))*expm(-1i*H(:,:,1)*x(1))*psi_in).^2-abs(psi_out).^2));
        Lovi=[L1 L2 L3];
        x = fminsearch(fun,Lovi);
        La=x;
        izlaz=expm(-1i*H(:,:,3)*La(3))*expm(-1i*H(:,:,2)*La(2))*expm(-1i*H(:,:,1)*La(1))*psi_in;
    case 9
        L=zeros(fix(n/2),broj);
        for i=1:fix(n/2)
            L(i,:)=(od_koliko:korak:do_koliko)*pi;
        end
        psi_in=zeros(1,n)';
        psi_in(fix(n/2)+1)=1;
        sredina=fix(n/2)+1;
        H=zeros(n,n,fix(n/2));
        H(sredina,sredina-1,1)=1;
        H(sredina,sredina+1,1)=1;
        H(sredina-1,sredina,1)=1;
        H(sredina+1,sredina,1)=1;
        for ii=2:fix(n/2)
            dole_desno=sredina+(ii-1);
            gore_levo=sredina-(ii-1);
            H(gore_levo,gore_levo-1,ii)=1;
            H(gore_levo-1,gore_levo,ii)=1;
            H(dole_desno,dole_desno+1,ii)=1;
            H(dole_desno+1,dole_desno,ii)=1;
        end
        min=10;
        for i=1:broj
            for j=1:broj
                for k=1:broj
                    for m=1:broj
                        psi_outA=expm(-1i*H(:,:,4)*L(4,i))*expm(-1i*H(:,:,3)*L(3,j))*expm(-1i*H(:,:,2)*L(2,k))*expm(-1i*H(:,:,1)*L(1,m))*psi_in;
                        if min > sum(abs(abs(psi_out).^2-abs(psi_outA).^2))
                            min=sum(abs(abs(psi_out).^2-abs(psi_outA).^2));
                            L1=L(1,m);
                            L2=L(2,k);
                            L3=L(3,j);
                            L4=L(4,i);
                        end
                    end
                end
            end
        end
        fun =@(x) sum(abs(abs(expm(-1i*H(:,:,4)*x(4))*expm(-1i*H(:,:,3)*x(3))*expm(-1i*H(:,:,2)*x(2))*expm(-1i*H(:,:,1)*x(1))*psi_in).^2-abs(psi_out).^2));
        Lovi=[L1 L2 L3 L4];
        x = fminsearch(fun,Lovi);
        La=x;
        izlaz=expm(-1i*H(:,:,4)*La(4))*expm(-1i*H(:,:,3)*La(3))*expm(-1i*H(:,:,2)*La(2))*expm(-1i*H(:,:,1)*La(1))*psi_in;
end
end