parameters = struct();
parameters.type = 'ConvolutionAndDownsample';
parameters.channelNumber = 5;
parameters.downsamplePhases = [1,2,3,2,4];
parameters.sensorKernels=normalize([1,2,5,2,1]);
parameters.kernels=cell(5,1);
parameters.kernels{1}=Gauss_1D(20);
parameters.kernels{2}=Gauss_1D(23);
parameters.kernels{3}=Gauss_1D(17);
parameters.kernels{4}=Gauss_1D(19);
parameters.kernels{5}=Gauss_1D(21);
parameters.downsampleCoefficient = 10;
inputdata = load('Entry001.csv');
Y = inputdata(:,2)';
result= getDegradedData(Y,parameters);
plot(result);