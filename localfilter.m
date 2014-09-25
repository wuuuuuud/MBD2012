ufiltered=u;
for j=1:100
for i=1:len
    if u(i)>=0
        
    else
        if i>3 && i<len-2
            ufiltered(i-1)=ufiltered(i-1)+0.3*u(i);
            ufiltered(i-2)=ufiltered(i-2)+0.16*u(i);
            ufiltered(i+1)=ufiltered(i+1)+0.3*u(i);
            ufiltered(i+2)=ufiltered(i+2)+0.16*u(i);
            ufiltered(i+3)=ufiltered(i+3)+0.14*u(i);
            ufiltered(i-3)=ufiltered(i-3)+0.14f*u(i);
            ufiltered(i)=ufiltered(i)-1.2*u(i);
        end
    end
end
u=ufiltered;
end