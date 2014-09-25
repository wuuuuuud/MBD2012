for j=19:29
    hindex=floor((j+1)/2);
    uindex=ceil((j+1)/2);
    for i=1:channelN
        H((i-1)*len+1:i*len,:)=generateCM(h((i-1)*kerlen+1:i*kerlen),len,len);
    end
    
    gamma/2*norm(H*uhistory{uindex}-g)^2+sum(abs(Dx*uhistory{uindex}))+delta/2*hhistory{hindex}'*G_Delta*hhistory{hindex}
    gamma/2*norm(H*uhistory{uindex}-g)^2
    sum(abs(Dx*uhistory{uindex}))
end