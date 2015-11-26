function [A] = makePauliMeasurements(n,q)

d=2^n;
A = [];

[M] = permn([1 2 3 4], n);

ss = size(M);

s{1} = [1 0;0 1];
s{2} = [1 0; 0 -1];
s{3} = [0 1;1 0];
s{4} = [0 -1i;1i 0];

A = [];

for ii=1:q
    
    if ii == 1
        U = eye(2);
    else
        X = (randn(2) + 1i*randn(2))/sqrt(2);
        [Q,R] = qr(X);
        R = diag(diag(R)./abs(diag(R)));
        U = Q*R;
    end
    
    p{1} = U'*[1 0;0 1]*U;
    p{2} = U'*[1 0; 0 -1]*U;
    p{3} = U'*[0 1;1 0]*U;
    p{4} = U'*[0 -1i;1i 0]*U;
    
    for i=1:ss(1)
        pos = M(i,:);
        
        m=1;
        for ii=1:n
            m=kron(m,p{pos(ii)});
        end
        
        v=  reshape(m,[1 d^2]);
        A = [A; v];
        
    end
end





% figure(1)
% imagesc(abs(AA))

