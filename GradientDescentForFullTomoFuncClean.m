%Apply quantum state recovery with the gradient descent method

function [rhor Xhi2 puri fid] = GradientDescentForFullTomoFuncClean(...
    A,counts,init,actualState,actualPrimary,pData,alpha,iMax,...
    gam,method,alpha2,poiss,smolinWay);

s = size(A);
N=s(1);
d = sqrt(s(2));
AA = A'*A;

tau =0.;

%Calulate probabilities and add Poissonian noise to them
P = pData(1:N)*alpha+1E-10;
% Pnoise = P';
Pnoise = reshape(P,[length(P) 1]);

%Prepare initial guess
% init=2;
if init==1    
    %best solution in the general space of matrices
     x = pinv(AA+eye(length(AA))/10*0)*A'*Pnoise;
    % x = rand(d^2,1)*counts;
elseif init==2 
    %random state produces with projector matrix
    x = rand(length(AA),1)*max(Pnoise);
end

%%% start with best pure state
% rhox = reshape(x,[sqrt(length(x)) sqrt(length(x))]);
% [primary maxEig] = primaryEigenvector(rhox);
% x = primary*primary'*maxEig;
% eig(x)
% x = reshape(x,[length(AA) 1]);

%start with closest state in the 2-norm
rhor = reshape(x,[d d]);
  rhorSmo = smolin(rhor);
  x = reshape(rhorSmo,[d^2 1]);
  
  %see what the initial Smolin Fid and pur are on console
%      fidSmolin = fidelityRho(actualState,rhorSmo /trace(rhorSmo ))
%     purityRe = trace((rhorSmo /trace(rhorSmo ))^2);
%     puriSmolin = purityRe



%%%%% try to make the purity higher because the initial guess is always too
%%%%% mixed
%     rhor = reshape(x,[d d]);
%     [eigvec eigval]= eig(rhor);
%     eigval = eigval - sum(sum(eigval))*0.05;
%     eigval(eigval<0) = 0;
%     eigval = eigval/sum(sum(eigval));
%     rhor = eigvec*eigval*eigvec'*counts;
%     x = reshape(rhor,[d^2 1]);

x_nplus1 = x;
previousX = x;
gradFromX = x*0;
y = x;

% Xhi2(1)=1E20;

%Gradient descent method
for i=1:iMax

    rhor = reshape(x_nplus1,[d d]);
    
    %Force legitimacy of density matrix 
    if smolinWay==1
        rhor = rhor/trace(rhor);
        ddd = eig(rhor);
     rhor = smolin(rhor);
%         
%         ddd- eig(rhor)
    else
        [eigvec eigval]= eig(rhor);
        eigval(eigval<tau) = 0;
        sumEig(i) = sum(sum(abs(eigval)));
        eigval = eigval/sum(sum(eigval));
        rhor = eigvec*eigval*eigvec';  
        rhor = rhor/trace(rhor);
    end
   
    fid(i) = fidelityRho(actualState,rhor);
    purityRe = trace(rhor^2);
    puri(i) = purityRe;
    
    %Change the test state in the opposite side of the gradient
    x = reshape(rhor,[d^2 1])*counts;
    Ax = A*x;
    Axb = (Ax-Pnoise);
    
    if poiss==1
        %gradient for poisson noise
        Xhi2(i) = real(sum(Ax-Pnoise.*log(Ax)));
        gradient =  A'*(1-Pnoise./Ax);
    elseif poiss==2
        %gradient for gaussian noise
        Xhi2(i+1) = sum(abs(Axb).^2./abs(Ax))/length(Pnoise);
        gradient = A'*((Axb./Ax).*(2-(Axb./Ax)));
    else
        Xhi2(i+1) = sum(abs(Axb).^2);
        gradient = 2*A'*Axb;
    end

    switch method
        case 'gradientDescent'
            x_nplus1 = x-gam.*gradient; 
        case 'momentum'
            if i==1 v{1} = gradient*eps; 
                    v2{1} = v{1};
            end
            v{i+1} = alpha2*v{i} - gam*gradient ;
            x_nplus1 = x + v{i+1}  ;
    end

    previousGrad = gradient;
    previousX = x;
   
    if Xhi2(i)==min(Xhi2) 
        bestrho=rhor;  
    end
end

rhor = bestrho;

%  rhor = smolin(rhor);
%     fid(i) = fidelityRho(actualState,rhor);
%     purityRe = trace(rhor^2);
%     puri(i) = purityRe;
%     
   
%     [eigvec eigval]= eig(rhor);
%      eigval(eigval<tau) = 0;
%     sumEig(i) = sum(sum(abs(eigval)));
%      eigval = eigval/sum(sum(eigval));
%     rhor = eigvec*eigval*eigvec';  
    
     
    rhor = rhor/trace(rhor);


% fid(i) = fidelityRho(actualState,rhor)
% x = reshape(rhor,[d^2 1])*counts;
% Ax = A*x;
% Axb = (Ax-Pnoise);
% Xhi2(i+1) = sum(abs(Axb).^2);







  
    
