function [A] = makeQubitMeasurementsRand(n,N)

A = [];
for ii=1:N
    for i=1:n
        ai{i} = rand(1,2).*exp(1i*2*pi*rand(1,2));
        ai{i} = ai1/sqrt(ai1*ai1');
    end
    
    aq=1;
    for i=1:n
        aq = kron(aq,ai{i});
    end
    ak = kron(aq,conj(aq));
    A = [A;ak];
end

