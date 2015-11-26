function SICpure(d,N)


% d=20;
alpha = 0.8;
gam = 0.01;

 numIt = 1000;
% beta = 10;
% N = beta*d;
% N = d*4-4;

AdN = ['A' num2str(d) '_N' num2str(N)];

 beta = d/N;
 
% ref = ones(d^2)/(d+1);
%
% ref = ref - eye(d^2)/(d+1)+ eye(d^2);

% mu = (beta-1)/(beta*d-1);


mu = (1-beta)/beta/(N-1);
mu = (N-d)/d/(N-1)

% mu = 0.0404;

ref = ones(N)*mu;
ref = ref - eye(N)*mu+ eye(N);


minf = 100;
samples = 10;
for k=1:samples
    A=[];
    for i=1:N;
        vec = rand(1,d).*exp(1i*2*pi*rand(1,d));
        vec = vec/sqrt(vec*vec');
        A = [A;vec];
    end
    
   
    
    
    % load('bestA3','A');
    
    %   A = A + normrnd(0,0.3,[d^2 d])+1i*normrnd(0,0.3,[d^2 d]);
    for n=1:N
        A(n,:) = A(n,:)/sqrt(A(n,:)*A(n,:)');
    end
    
    Ainit = A;
    
    v{1} = A*0;
    iterations = [1:numIt];
    a=0;
    for i=iterations
%         if i<1000
%             gam = 0.001;
%         else
%             gam=0.00005;
%         end
        a=a+1;
        
        AA = A*A';
        f(a) = sum(sum((abs(AA).^2 - ref).^2));
        
        
        %
        fm = (abs(AA).^2 - ref);
        % %     gradient = 4*(A')*((AA).*fm);
        gradient = 4*A'*((AA).*fm);
        
        
        
        v{i+1} = alpha*v{i} - gam*(gradient')/norm(gradient);
        %         v{i+1} = alpha*v{i} - gam*(gradient');
        A = A + v{i+1};
        
        %          for n=1:d^2
        %         A(n,:) = A(n,:)/sqrt(A(n,:)*A(n,:)');
        %         end
    end
    f(end)
    if f(end)<minf
        minf = f(end);
        bestA = A;
    end
end

A =bestA;
minf

save('bestA3','A');

figure(193)
imagesc(abs(gradient'))

for n=1:N
    A(n,:) = A(n,:)/sqrt(A(n,:)*A(n,:)');
end

figure(188)
imagesc(abs(A*A').^2)

M = abs(A*A').^2;

a=0;
for i=2:length(M-1)
    for ii=1:i-1
        a=a+1;
        q(a) = M(i,ii);
    end
end

mq.me = mean(q);
mq.st = std(q)



figure(190)
plot(iterations,(abs(f)))

% figure(191)
% imagesc(abs(grad))

figure(192)
imagesc(abs(Ainit*Ainit').^2)




save(AdN,'A')
% v{1} = A*0;
% iterations = [1:numIt];
% for i=iterations
%     AA = A*A';
%     fm = (abs(AA).^2 - ref);
% %     gradient = 4*(A')*((AA).*fm);
%     gradient = 4*(conj(AA).*fm)*conj(A);
%     f(i) = sum(sum(abs(fm).^2));
%     v{i+1} = alpha*v{i} - gam*(gradient);
%
%
%     A = A + v{i+1};
%
%    if f(i)<0.1
%        break;
%    end
%
% end
%
%  figure(188)
%     imagesc(abs(A*A').^2)
%
%     figure(189)
%     plot(iterations,log(log(abs(f))))
%
%     figure(190)
% plot(iterations,(abs(f)))
% % imagesc(ref)
%
%
%
% sum(abs(A').^2)
%
