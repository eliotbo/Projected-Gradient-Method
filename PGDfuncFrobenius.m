%Apply quantum state recovery with the gradient descent method

function [rhor Xhi2 puri fid primaryFid] = PGDfuncFrobenius(...
    A,counts,init,actualState,actualPrimary,pData,alpha,iMax,gam,method,alpha2);

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
if init==1    
    %best solution in the general space of matrices
     x = pinv(AA+eye(length(AA))/10)*A'*Pnoise;
    % x = rand(d^2,1)*counts;
elseif init==2 
    %random state produces with projector matrix
    x = rand(length(initialGuess),1)*max(Pnoise);
end

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

Xhi2(1)=1E20;

%Gradient descent method
for i=1:iMax
    rhor = reshape(x_nplus1,[d d]);
    
%     %Force the density matrix to be physical
%     [eigvec eigval]= eig(rhor);
%     eigval(eigval<tau) = 0;
%     sumEig(i) = sum(sum(abs(eigval)));
%      eigval = eigval/sum(sum(eigval));
%     rhor = eigvec*eigval*eigvec';  
%     rhor = rhor/trace(rhor);
%     primary = eigvec(:,find(diag(eigval)==max(diag(eigval))));
    
    rhor = initNearestSPD3(rhor);
    rhor = rhor/trace(rhor);
    [eigvec eigval]= eig(rhor);
    primary = eigvec(:,find(diag(eigval)==max(diag(eigval))));
    
    primaryFid(i) = abs(primary'*actualPrimary);
    fid(i) = fidelityRho(actualState,rhor);
    purityRe = trace(rhor^2);
    puri(i) = purityRe;
    
    %Calculate overlap between actual state and test state
    testPsi = eigvec(:,find(sum(eigval)==max(sum(eigval))));
    
    %Change the test state in the opposite side of the gradient
    x = reshape(rhor,[d^2 1])*counts;
    Ax = A*x;
    Axb = (Ax-Pnoise);
    Xhi2(i+1) = sum(abs(Axb).^2./abs(Ax))/length(Pnoise);
    gradient = A'*((Axb./Ax).*(2-(Axb./Ax)));
         
    switch method
        case 'gradientDescent'
            x_nplus1 = x-gam.*gradient; 
        case 'momentum'
            if i==1 v{1} = gradient*eps; 
                    v2{1} = v{1};
            end

%             vNormed = v{i}/norm(v{i});
%             gradientNormed = gradient/norm(gradient);
%             over(i) = abs(vNormed'*gradientNormed);
%             over2(i) = abs(gradFromX'*gradientNormed);
            v{i+1} = alpha2*v{i} - gam*gradient ;
            x_nplus1 = x + v{i+1}  ;
           
    end

    previousGrad = gradient;
    previousX = x;

%     if mod(i,100)==0; 
%         disp(['Xhi2: ' num2str(Xhi2(end))]); 
%     end
   
    if Xhi2(i+1)==min(Xhi2) 
        bestrho=rhor;  
    end
    
end

rhor = bestrho;
rhor = rhor / trace(rhor);





  
    
