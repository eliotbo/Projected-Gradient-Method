function A = sicQubits(n,q)



d=2^n;
A = [];

[M] = permn([1 2 3 4], n);

ss = size(M);

% s{1} = [1 0;0 1];
% s{2} = [1 0; 0 -1];
% s{3} = [0 1;1 0];
% s{4} = [0 -1i;1i 0];

W =[ 1.0000 + 0.0000i   0.0000 + 0.0000i;
  -0.5774 + 0.0000i  -0.8165 + 0.0000i;
  -0.2887 - 0.5000i   0.8165 + 0.0000i;
  -0.2887 + 0.5000i   0.8165 + 0.0000i]

% W = [1 0; 
%     0 1; 
%     1/sqrt(2) 1/sqrt(2); 
%     1/sqrt(2) -1/sqrt(2); 
%     1/sqrt(2) 1i/sqrt(2); 
%     1/sqrt(2) -1i/sqrt(2)];


for i=1:4
    v= W(i,:);
   s{i} =  v'*v;
end


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
    
    for k=1:length(s)
        p{k} = U'*s{k}*U;
    end

    
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


% W = [1 0; 
%     1/sqrt(2) 1/sqrt(2); 
%     1/sqrt(2) exp(1i*4*pi/3)/sqrt(2); 
%     1/sqrt(2) exp(-1i*4*pi/3)/sqrt(2)];
% 
% r= 1/sqrt(2);
% t=sqrt(1-r^2);
% W = [1 0; 
%     r t; 
%     r t*exp(1i*4*pi/3); 
%     r t*exp(-1i*4*pi/3)];
% 
% G = zeros(2);
% for i=1:4
%    v = W(i,:);
%    if i==1
%        yo=1;
%    else
%        yo = 2/3;
%    end
%    m = yo*v'*v;
%    M{i} = m;
%    G = m + G;
% end
% 
% R=[];
% for i=1:4
%    Q = G^(-1/2)*M{i}*G^(-1/2);
%    P{i}=Q;
%    v = primaryEigenvector(Q);
%    R = [R;transpose(v)];
% end
% 
% T=zeros(2);
% for i=1:4
%    T=  P{i}+T;
% end
% 
% abs(R)
% R
% 
% figure(134)
% imagesc(abs(R*R').^2);