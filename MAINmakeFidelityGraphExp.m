%plot the experimental fidelity graph, ''offset'' selects a subset of all
%measurements
clear;clc;

%%%%%%%%%%%% parameters used for demonstration shown in paper %%%%%%%%%%
% N=1000; %%%number of measurements
% experiment = 0; %%%Do not take experimental data
%%%Correction factor if take experimental data (should be about 0.93)
% alpha = 1;
% gam=0.3; %%% descent speed
% iMax =  2000; %%% number of iterations
% counts = 50000; % number of counts -> Tr(rho) when rho is not normlalized
% d=25; %%% dimentionality
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%% parameters used for demonstration shown in paper %%%%%%%%%%
% N=400;
% experiment = 0;
% alpha = 1;
% gam=0.3;
% 
% iMax =  2000;
% counts = 32000;
% d=16;
%%%%%%%%%%%% parameters used for demonstration shown in paper %%%%%%%%%%

% N=1500;
% experiment = 0;
% alpha = 1;
% gam=0.3;
% 
% iMax =  2000;
% counts = 72000;
% d=36;

N=1000;
experiment = 1;
%experiemental correction on counts
alpha = 0.93;
friction = 0.85;
gam=0.05;
numSamples = 50;

iMax =  1000;
counts = 800;
d=25;
poiss = 1;
smolinWay = 1;
% measurementType = 'pauli3qubits';    %1256
% measurementType = 'qubits'; %256
measurementType ='mub3d';
% measurementType ='sic3';

%Generate completely random density matrix with predefined eigenvalues
lambda = zeros(1,9);
x = [1:d];
lambda=0;
purity = 0;

PUR = linspace(1/d,0.999,5);
% PUR = linspace(1/d+0.1,0.999,5);
%     PUR = 0.8;
a=0;
for pur = PUR
    a=a+1;
    while purity<pur
        lambda = lambda + 0.001; %increase std until reach correct purity
        lam = exp(-lambda*x); %exponential distribution of eigenvalues
        lamb = lam/sum(lam);
        purity = sum(lamb.^2);
    end
    lambi(a,1:d) = lamb;
end


a=0;
offset = 0;
% for offset = 0:1000:3000
tic
L= length(PUR);
for i=1:L
a=a+1;
lamb = lambi(a,:);


for k=1:numSamples
%%%%%% Gradient descent %%%%%
tic
    [thisXhi21 thisPurity1 fidelity1  A totalPseudoCounts]...
    = GD_FullTomo_init2(N,offset,experiment,alpha,gam,...
                    iMax,counts,lamb,poiss,measurementType,friction,smolinWay);
                toc

    thisXhi2(k) = thisXhi21 ;
    thisPurity(k) = thisPurity1 ;
    fidelity(k) = fidelity1 ;

end
%%%%%% eigenvector finder %%%%%%%
%     [thisXhi2 thisPurity fidelity primaryFidelity rhoVecs, rhoEig, rho]...
%     = eigenInitialize(N,offset,experiment,alpha,counts,lamb);
    xhi2(a) = mean(thisXhi2);
    purity(a) = mean(thisPurity);
    puritystd(a) = std(thisPurity);
    fid(a) = mean(fidelity);
    fidstd(a) = std(fidelity);
%     prim(a) = mean(primaryFidelity);
%     primstd(a) = std(primaryFidelity);
    
    
    timeSoFar = toc;
    timePerIteration = timeSoFar/a;
    
    timeRemaining= timePerIteration*(L-a);
    timeRemainingMin = timeRemaining/60;
    timeRemainingHours = timeRemaining/3600;
    disp(['Time remaining: ' num2str(timeRemainingMin,3) 'min' ' or ' num2str(timeRemainingHours,3) ' hours']);
%     if fidelity<0.99
%         break
%     end

end

pos = [1:a];

totalPseudoCounts

figure(545)
hold off
plot(PUR,abs(fid(pos)),'-')
hold on
plot(PUR,abs(purity(pos)),'r-')
% plot(PUR,abs(prim(pos)),'g-')
plot(PUR,PUR,'black')
xlabel('Purity');ylabel('pFid->green, fid->blue, pur-red')
% axis([0 8 0.9 1])

figure(198)
hold on
errorbar(PUR,abs(fid(pos)),abs(fidstd),'b--')
hold on
errorbar(PUR,abs(purity(pos)),abs(puritystd),'r--')
% errorbar(PUR,abs(prim(pos)),'g-')
plot(PUR,PUR,'black')
xlabel('Purity');ylabel('pFid->green, fid->blue, pur-red')
