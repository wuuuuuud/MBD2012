function [ result ] = isSingleSensorKernel( sensorKernels )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
sensorKernelsSize = size(sensorKernels);
if sensorKernelsSize(1) == 1
    result = true;
else
    result = false;
end

end

