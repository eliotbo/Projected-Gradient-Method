clear;clc;

rho = makeRandomDensityMatrix([1 0 0 0 0 0 0 0 ]);
trace(rho)
rho = rho*1000;
x = reshape(rho,[8^2 1]);

[A P] = make3QubitHVADRLandPauli(x);
size(A)