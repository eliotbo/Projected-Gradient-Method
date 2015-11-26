function [rho] = makeRandomDensityMatrix(lambda)

  %generate random density matrix
    d = length(lambda);
    rho = zeros(d);
    randM = rand(d).*exp(1i*2*pi*rand(d));
    randU = expm((randM-randM'));
   
    for ii=1:d
        psi = randU(:,ii);
        rho = rho + psi*psi'*lambda(ii);
    end
    
    