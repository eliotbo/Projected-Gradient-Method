%Full tomography via gradient descent with renormalisation and positivity 
%of eigenvalues.
%
%Parameters:
%experiment:1 if using experimental data. 0 if simulation.
%lambda:    set of eigenvalues of simulated density matrix
%N:         number of projectors used in reconstruction procedure
%offset:    specifies the subset of projectors used in experiment
%initialiseAll:      
%           set to 1 if initial guess is pinv(A'*A)*A'*probs
%           set to 0 if load last guess
%           set to 2 if random initial guess.
%alpha:     multiplies the probabilities. In pratice the probabilities
%           might not be normalised. Alpha is a compensatory degree of freedom.
%counts:    maximum number of counts in an outcome.


function [thisXhi2 thisPurity fidelity A totalPseudoCounts]...
 = GD_FullTomo_init2(N,offset,experiment,alpha,gam,iMax, ...
                     counts,lambda,poiss,measurementType,friction,smolinWay)


method = 'momentum'; 
% method = 'gradientDescent';
%set init to 1 for initial guess to be close to actual state ((A'A)^(-1)*A'*p)
initialiseAll = 1;

alpha2 = friction;
% % alpha2 = 0.5;

d= length(lambda);
loadname = 'bestRhof4.mat';
savename = 'bestRhof4.mat';
%with counts = 1200, offset = 1000, N = 1000, alpha = 0.93,
%the best xhi2 is 1.8641

rho = makeRandomDensityMatrix(lambda);
    rhovv = reshape(rho,[d^2 1]);

if experiment==1
    %retrieve directly measured state vector
    alpha = 0.93 ;
    pureRe = dlmread('pureStateRe');
    pureIm = dlmread('pureStateIm');
    pureM = pureRe + 1i*pureIm;
    pureM = pureM / sqrt(sum(sum(abs(pureM).^2)));
    pureM = pureM*exp(-1i*angle(pureM(1,1)));
    pureV = reshape(pureM,[d 1]);
    rho = pureV*pureV';
    rho = rho / trace(rho);
    load('fullTomoData2.mat');
    A = reimA(offset+1:N+offset,:);
    pData = co([offset+1:N+offset])*counts;
    totalPseudoCounts = sum(pData);
else
    switch measurementType
        case 'mub25'
            [A ]= MUB5GD(1,6);
        case 'rand2d'
            A =   makeRandom2dProjections(d,N);
        case 'mub3d'
             A =   MUB3GD();
        case 'sic3'
            A = sic3();
        case 'qubits'
             [A] = makeQubitMeasurementsHVADRL(log2(d));
             save('A3qubits','A');
           % load('A4qubits'); %this is the matrix for 5 qubits (not 4!)
%            load('A3qubits');
        case 'randQubits'
            [A] = makeQubitMeasurementsRand(log2(d),N);
        case 'sicQubits'
            A = sicQubits(log2(d),1);
        case 'pauli'
            [A] = makePauliMeasurements(log2(d),1);
        case 'pauli3qubits'
            [A pauliNoise] = make3QubitHVADRLandPauli(rho*counts);
    end
%     rho = makeRandomDensityMatrix(lambda);
%     rhovv = reshape(rho,[d^2 1]);
    
    size(A)
    size(rhovv)
    %generate fake data
    pDataExact = real(A*rhovv)*counts;
    pData = poissrnd(real(pDataExact));
    totalPseudoCounts = sum(pDataExact);
%     pData = pDataExact+randn(length(pData),1)*0.;

    length(pData)
    switch measurementType
        case 'pauli'
            r = (counts - pDataExact)/2;
            t = (counts + pDataExact)/2;
            rNoisy = poissrnd(r);
            tNoisy = poissrnd(t);
            pData(2:end) = tNoisy(2:end) - rNoisy(2:end) ;
            totalPseudoCounts = counts*length(r);
%           pData = pDataExact + randn(length(pDataExact),1)*counts/10;
%           (pData-pDataExact)./pDataExact/counts
%           mean(abs((pData-pDataExact*counts)./pDataExact)/counts)
        case 'pauli3qubits'
            pData = pauliNoise;

    end
    
end
% if sum(sum(isnan(pData)))
%     pureM = primaryEigenvector(rho);


pureInit = trace(rho^2);
% disp(['Initial purity:' num2str(pureInit)]);
if initialiseAll==0
    load(loadname); else x = 0; lastXhi2 = 1E20;
end
pureM = primaryEigenvector(rho);

% x = reshape(rho,[d^2 1])*counts;
% % x = rand(d^2,1);

actualState = rho;

[rhof xhi2 puri fid] = GradientDescentForFullTomoFuncClean(A,counts,initialiseAll,...
    actualState,primaryEigenvector(rho),pData,alpha,iMax,gam,method,alpha2,poiss,smolinWay);

% [rhof xhi2 puri fid primaryFid]  = PGDfuncFrobenius(A,counts,initialiseAll,...
%      actualState,primaryEigenvector(rho),pData,alpha,iMax,gam,method,alpha2);

% [rhof xhi2 puri gradXB] = nonconvexFullTomoFunc(A,counts,initialiseAll,...
%     x,pData,alpha,iMax,gam,method,alpha2);



thisXhi2 = xhi2(end);
thisPurity = puri(end);

x = reshape(rhof,[d^2 1])*counts;

% if thisXhi2 < 10*lastXhi2
%     disp('Saving best guess.')
%     lastXhi2 = thisXhi2;
%     save(savename,'x','A','pData','lastXhi2','rho')
% else
%     disp('NOT SAVED because Xhi2 was too high')
% end

disp(' ')
disp(['Final Xhi2: ' num2str(xhi2(end))]);
disp(['Final purity: ' num2str(abs(thisPurity))]);


% primaryVec = primaryEigenvector(rhof);
% fidelityDir = directV'*rhof*directV;

fidelity = fidelityRho(rhof,rho);
% fidelity2 = fidelityRho(rhof2,rho);

% primaryFidelity = abs(pureM'*primaryVec)^2;

disp(['Final Fidelity: ' num2str(abs(fidelity))]);
% disp(['Final primary fidelity: ' num2str(primaryFidelity)]);
disp(['Initial Fidelity: ' num2str(abs(fid(1)),3)])
disp(['Initial Purity: ' num2str(abs(puri(1)),3)])

% disp(['Final Fidelity benchmark: ' num2str(abs(fidelity2))]);

figure(105)
hold off
plot(abs(puri),'r')
ylabel('purity and fidelity')
hold on
plot(abs(fid),'b')
% plot(abs(primaryFid),'g')
% legend('Purity','Fidelity','Primary Fidelity')

% figure(155)
% hold off
% plot(abs(puri2),'r')
% ylabel('purity and fidelity')
% hold on
% plot(abs(fid2),'b')
% plot(abs(primaryFid2),'g')

% figure(109);
% hold off
% plot(abs(primaryVec).^2);
% hold on
% plot(abs(pureM).^2),'r';
% 
% figure(104);
% hold off
% plot(angle(primaryVec));
% hold on
% plot(angle(pureM),'r');


% figure(104);
% subplot(2,1,1);
% imagesc(angle(reshape(primaryVec,[sqrt(d) sqrt(d)])));
% colormap hsv
% caxis([-pi, pi])
% subplot(2,1,2);
% imagesc((angle(pureM)));
% colormap hsv
% caxis([-pi, pi])



figure(2)
plot((xhi2(2:end)))
ylabel('\chi^2')
