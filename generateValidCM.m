function [ m ] = generateValidCM( ker,srclen )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
m= generateCM(ker,srclen,abs(length(ker)-srclen)+1);
return;

end

