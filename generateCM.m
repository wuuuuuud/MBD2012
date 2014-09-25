function [ m ] = generateCM( ker,srclen,dstlen )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
ker=fliplr(ker);
ker=flipud(ker);
kerlen=length(ker);
initP=round(0.5*(-(dstlen+kerlen)+srclen+3));
m=zeros(dstlen,srclen);
if kerlen > srclen
    disp('creating a convolve matrix, expecting a dense one');
    for i=1:dstlen
        for j=1:srclen
            index=-(i+initP-1)+j+1;
            if 0<index && index<=kerlen
                m(i,j)=ker(index);
            end
        end
    end
else
    disp('creating a convolve matrix, expecting a sparse one');
    for i=1:kerlen
        for j=1:dstlen
            position=i+initP+j-1-1;
            if position>0 && position<=srclen
                m(j,position)=ker(i);
            end
        end
    end
end

end

