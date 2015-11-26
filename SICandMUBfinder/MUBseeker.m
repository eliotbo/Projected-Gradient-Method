clear;clc;
d=3;
alpha = 0.95;
gam = 0.01;



ref = ones((d+1)*d,(d+1)*d)/d;

for i=1:(d+1)
      ref((1:d)+d*(i-1),(1:d)+d*(i-1)) = 0;
end

  ref = ref + eye((d+1)*d,(d+1)*d);
 
figure(111)
imagesc(ref)

A=[];
for i=1:(d+1)*d;
    vec = rand(1,d).*exp(1i*2*pi*rand(1,d));
    vec = vec/sqrt(vec*vec');
    A = [A;vec];
end

numIt = 400;


%    load('bestA3','A');
  
   for i=1:(d+1);
    UN = A((d)*(i-1) +(1:d),1:end) ;
    [Q R] = qr(UN');
    A((d)*(i-1) +(1:d),1:end) =Q';
end
   
%   A = A + normrnd(0,0.1,[d*(d+1) d])+1i*normrnd(0,0.1,[d*(d+1) d]);
   for n=1:d^2
        A(n,:) = A(n,:)/sqrt(A(n,:)*A(n,:)');
   end
 
Ainit = A;

v{1} = A*0;
iterations = [1:numIt];
a=0;
for i=iterations
    
a=a+1;
    
     AA = A*A';
     f(a) = sum(sum((abs(AA).^2 - ref).^2));
     

             fm = (abs(AA).^2 - ref);
% %     gradient = 4*(A')*((AA).*fm);
            gradient = 4*A'*((AA).*fm);
     

     
        v{i+1} = alpha*v{i} - gam*(gradient')/norm(gradient);
%         v{i+1} = alpha*v{i} - gam*(gradient');
        A = A + v{i+1};
        
         for n=1:(d+1)*d
        A(n,:) = A(n,:)/sqrt(A(n,:)*A(n,:)');
     end
end

   f(end)
save('bestA3','A');

figure(193)
imagesc(abs(gradient'))

 for n=1:(d+1)*d
        A(n,:) = A(n,:)/sqrt(A(n,:)*A(n,:)');
 end

figure(188)
imagesc(abs(A*A').^2)

figure(190)
plot(iterations,(abs(f)))

% figure(191)
% imagesc(abs(grad))

figure(192)
imagesc(abs(Ainit*Ainit').^2)

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
