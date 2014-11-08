function [ m ] = getDownSampleMatrix( srclen,dstlen,downSampleFactor,phase )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
m=sparse(1:(1+floor((srclen-phase)/downSampleFactor)),phase:downSampleFactor:srclen,1,dstlen+2,srclen);
m=m(1:dstlen,1:srclen);

end

