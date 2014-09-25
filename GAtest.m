    useeds=zeros(len,100);
    for n=1:100
        useeds(:,n)=u;
    end
    for n=1:20
        useeds(:,10:100)=useeds(:,10:100).*(rand(len,91)*0.08+0.96);
        useeds=0.5*( useeds+abs(useeds));
        for m=1:100
            uenergy(m)=gamma/2*norm(H*useeds(:,m)-g)^2+sum(abs(Dx*useeds(:,m)));
        end
        [AAA,uorder]=sort(uenergy);
        ubuff=useeds;
        for m=1:100
            useeds(:,m)=ubuff(:,uorder(1+mod(m,8)));
        end
        n
    end