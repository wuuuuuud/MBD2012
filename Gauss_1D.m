function ker = Gauss_1D( HWHH )
%根据给定的半高宽的一半，生成一个高斯模糊核，并归一化
%input:  HWHH - half width of the half high
%output: ker - the result of gaussian bulry kernel

xx = [-round(3*HWHH):1:round(3*HWHH)];%离散化，并保证是对称的
ker = 1 / sqrt(2 * pi * HWHH^2) * exp(( -xx.^2)/(2 * HWHH^2 ));
nor = sum(sum( ker ));
ker = ker / nor;
end
