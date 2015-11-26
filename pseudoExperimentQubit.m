clear;clc;
n=2;
d = 2^n;
lambda = zeros(1,d);
lambda(1) = 1;
lambda = lambda/sum(lambda);
rho = makeRandomDensityMatrix(lambda);
[A] = makeQubitMeasurements(n,1);






