function [A] = makeQubitMeasurementsHVADRL(n)

d=2^n;
A = [];

[M] = permn([1 2 3 4 5 6], n);

ss = size(M);

% s{1} = [1 0;0 1];
% s{2} = [1 0; 0 -1];
% s{3} = [0 1;1 0];
% s{4} = [0 -1i;1i 0];


W = [1 0; 
    0 1; 
    1/sqrt(2) 1/sqrt(2); 
    1/sqrt(2) -1/sqrt(2); 
    1/sqrt(2) 1i/sqrt(2); 
    1/sqrt(2) -1i/sqrt(2)];


for i=1:6
    v= W(i,:);
   s{i} =  v'*v;
end


A = [];
    
% for k=1:length(s)
%     p{k} = s{k};
% end

for i=1:ss(1)
    pos = M(i,:);
    
    m=1;
    for ii=1:n
        m=kron(m,s{pos(ii)});
    end
    
    v=  reshape(m,[1 d^2]);
    A = [A; v];
end






% figure(1)
% imagesc(abs(AA))

