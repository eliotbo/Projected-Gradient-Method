function [A] = makeRandom2dProjections(d,N)

for ii=1:N
   
            ai1 = rand(1,2).*exp(1i*2*pi*rand(1,2));
            ai1 = ai1/sqrt(ai1*ai1');
            ai2 = rand(1,2).*exp(1i*2*pi*rand(1,2));
            ai2 = ai2/sqrt(ai2*ai2');
            aq1 = zeros(1,sqrt(d));
            aq2 = zeros(1,sqrt(d));
            aq1pos = ceil(rand(1,2)*sqrt(d));
            aq2pos = ceil(rand(1,2)*sqrt(d));
            while aq1pos(1) == aq1pos(2)
                aq1pos = ceil(rand(1,2)*sqrt(d));
            end
            while aq2pos(1) == aq2pos(2)
                aq2pos = ceil(rand(1,2)*sqrt(d));
            end
            aq1(1,aq1pos(1)) = ai1(1);
            aq1(1,aq1pos(2)) = ai1(2);
            aq2(1,aq2pos(1)) = ai2(1);
            aq2(1,aq2pos(2)) = ai2(2);
            aq1 = aq1+1E-10;
            aq1 = aq1/sqrt(aq1*aq1');
            aq2 = aq2+1E-10;
            aq2 = aq2/sqrt(aq2*aq2');
            ai = kron(aq1,aq2);
            A(ii,1:d^2) = kron(ai,conj(ai));
        end