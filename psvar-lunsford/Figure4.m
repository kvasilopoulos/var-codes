%==========================================================================
% This script produces Figure 4 of "Proxy SVARs: Asymptotic Theory, 
% Bootstrap Inference, and the Effects of Income Tax Changes in the United 
% States" by Carsten Jentsch and Kurt G. Lunsford.  
% 
% This script calls the following workspaces:
%   VARData.mat 
%   ProxyData.mat
% 
% This script call the following functions:
%   IdentifyFunc.m
%==========================================================================

clc
clear;
format short g
tic;

%==========================================================================
% Parameters of the Program
%==========================================================================
%model parameters
lags = 4;               %number of lags in the structural VAR
nImp = 20;              %number of steps in the IRF

%bootstrap parameters
nBoot = 10000;          %number of bootstrap replications
BlockSize = 19;         %size of blocks in the bootstrap
seed = 1;               %seed for random number generator
rng(seed);              %iniate the random number generator

%program toggles
toggleBoot = 1;         %1 to use moving block bootstrap
                        %2 to use wild bootstrap


%==========================================================================
% Load Data
%==========================================================================
%all VAR data from Mertens and Ravn (2013)
load VARDataAll.mat

%collect the relevant VAR data
VARData = [VARDataAll(:,1:2),VARDataAll(:,5:7),...
    VARDataAll(:,3),VARDataAll(:,14:15)];

%dimensions of the VAR data
[T,N] = size(VARData);

%tax proxy variables used in Mertens and Ravn (2013)
load ProxyData.mat
M = ProxyData(1+lags:T,:);  %adjust the proxy data for the number of lags
k = size(M,2);              %number of proxy variables


%==========================================================================
% Estimate the VAR
%==========================================================================
%left-hand side variables
LHS = VARData(1+lags:T,:);

%right-hand side variables
RHS = ones(T-lags,1);
for j = 1:lags
    RHS = [RHS,VARData(1+lags-j:T-j,:)];
end

%VAR coefficients and innovations
A_hat = (RHS'*RHS)\(RHS'*LHS);
U_hat = LHS - RHS*A_hat;


%==========================================================================
% Identify the Structural Shocks
%==========================================================================
H1_hat = IdentifyFunc(U_hat,M);


%==========================================================================
% Produce the IRFs
%==========================================================================
%impulse response functions for both the APITR and ACITR shocks
IRF = zeros(N,nImp,2);
for order = 1:2
    %normalize the shock to -1
    shock = -1/H1_hat(order,order);

    %generate the IRFs
    IRF(1:N,1,order) = H1_hat(:,order)*shock;
    history = [IRF(1:N,1,order);zeros((lags-1)*N,1)];
    for j = 2:nImp
        IRF(1:N,j,order) = A_hat(2:N*lags+1,:)'*history;
        history = [IRF(1:N,j,order);history(1:(lags-1)*N,1)];
    end
end


%==========================================================================
% Bootstrap and Confidence Intervals
%==========================================================================
%set up the blocks if using moving block bootstrap
if toggleBoot == 1
    %set number of blocks
    nBlock = ceil((T-lags)/BlockSize);

    %create the blocks
    Blocks = zeros(BlockSize,N,T-lags-BlockSize+1);
    MBlocks = zeros(BlockSize,k,T-lags-BlockSize+1);
    for j = 1:T-lags-BlockSize+1
        Blocks(:,:,j) = U_hat(j:BlockSize+j-1,:);
        MBlocks(:,:,j) = M(j:BlockSize+j-1,:);
    end

    %center the bootstrapped VAR errors
    centering = zeros(BlockSize,N);
    for j = 1:BlockSize
        centering(j,:) = mean(U_hat(j:T-lags-BlockSize+j,:),1);
    end
    centering = repmat(centering,[nBlock,1]);
    centering = centering(1:T-lags,:);

    %center the bootstrapped proxy variables
    Mcentering = zeros(BlockSize,k);
    for j = 1:BlockSize
        subM = M(j:T-lags-BlockSize+j,:);
        Mcentering(j,:) = [mean(subM((subM(:,1) ~= 0),1),1),...
            mean(subM((subM(:,2) ~= 0),2),1)];
    end
    Mcentering = repmat(Mcentering,[nBlock,1]);
    Mcentering = Mcentering(1:T-lags,:);
end

%arrays to store relevant variables
IRF_boot = zeros(N,nImp,nBoot,2);

%loop for the bootstrap replications
for boot = 1:nBoot
    if toggleBoot == 1
        %draw bootstrapped errors
        index = ceil((T - lags - BlockSize + 1)*rand(nBlock,1));
        U_boot = zeros(nBlock*BlockSize,N);
        M_boot = zeros(nBlock*BlockSize,k);
        for j = 1:nBlock
            U_boot(1+BlockSize*(j-1):BlockSize*j,:) = Blocks(:,:,index(j,1));
            M_boot(1+BlockSize*(j-1):BlockSize*j,:) = MBlocks(:,:,index(j,1));
        end
        U_boot = U_boot(1:T-lags,:);
        M_boot = M_boot(1:T-lags,:);

        %center the bootstrapped errors
        U_boot = U_boot - centering;
        for j = 1:k
            M_boot((M_boot(:,j)~=0),j) =...
                M_boot((M_boot(:,j)~=0),j) - Mcentering((M_boot(:,j)~=0),j);
        end

        %recursively generate data
        LHS_boot = zeros(T-lags,N);
        RHS_boot = [ones(T-lags,1),zeros(T-lags,N*lags)];
        RHS_boot(1,:) = RHS(1,:);
        for j = 1:T-lags-1
            LHS_boot(j,:) = RHS_boot(j,:)*A_hat + U_boot(j,:);
            RHS_boot(j+1,:) = [1,LHS_boot(j,:),RHS_boot(j,2:N*(lags-1)+1)];
        end
        LHS_boot(T-lags,:) = RHS_boot(T-lags,:)*A_hat + U_boot(T-lags,:);
    elseif toggleBoot == 2
        %generate the shocks for the wild bootstrap
        E = 2*(rand(T-lags,1) > 0.5) - 1;    

        %bootstrapped proxy
        M_boot = M.*(E*ones(1,k));

        %compute the bootstrapped shocks
        U_boot = U_hat.*(E*ones(1,N));

        %recursively generate data
        LHS_boot = zeros(T-lags,N);
        RHS_boot = [ones(T-lags,1),zeros(T-lags,N*lags)];
        RHS_boot(1,:) = RHS(1,:);
        for j = 1:T-lags-1
            LHS_boot(j,:) = RHS_boot(j,:)*A_hat + U_boot(j,:);
            RHS_boot(j+1,:) = [1,LHS_boot(j,:),RHS_boot(j,2:N*(lags-1)+1)];
        end
        LHS_boot(T-lags,:) = RHS_boot(T-lags,:)*A_hat + U_boot(T-lags,:);
    end

    %bootstrap VAR coefficients and innovations
    A_hatboot = (RHS_boot'*RHS_boot)\(RHS_boot'*LHS_boot);
    U_hatboot = LHS_boot - RHS_boot*A_hatboot;

    %identification
    H1_hatboot = IdentifyFunc(U_hatboot,M_boot);

    for order = 1:2
        %normalize the initial IRF shock to -1
        shock_boot = -1/H1_hatboot(order,order);

        %generate the IRFs
        IRF_boot(1:N,1,boot,order) = H1_hatboot(:,order)*shock_boot;
        history_boot = [IRF_boot(1:N,1,boot,order);zeros((lags-1)*N,1)];
        for j = 2:nImp
            IRF_boot(1:N,j,boot,order) =...
                A_hatboot(2:N*lags+1,:)'*history_boot;
            history_boot =...
                [IRF_boot(1:N,j,boot,order);history_boot(1:(lags-1)*N,1)];
        end
    end
end

%sort the bootstrapped IRFs
IRF_bootSort = zeros(N,nImp,nBoot,2);
for order = 1:2
    for n = 1:N
        for j = 1:nImp
            IRF_bootSort(n,j,:,order) = sort(IRF_boot(n,j,:,order));
        end
    end
end

%percentile confidence intervals
IRF_confidence68 = zeros(N,nImp,2,2);
IRF_confidence68(:,:,1,1) =...
    IRF_bootSort(:,:,round(nBoot*0.16,1));
IRF_confidence68(:,:,2,1) =...
    IRF_bootSort(:,:,round(nBoot*0.84),1);
IRF_confidence68(:,:,1,2) =...
    IRF_bootSort(:,:,round(nBoot*0.16),2);
IRF_confidence68(:,:,2,2) =...
    IRF_bootSort(:,:,round(nBoot*0.84),2);

IRF_confidence90 = zeros(N,nImp,2,2);
IRF_confidence90(:,:,1,1) =...
    IRF_bootSort(:,:,round(nBoot*0.05,1));
IRF_confidence90(:,:,2,1) =...
    IRF_bootSort(:,:,round(nBoot*0.95),1);
IRF_confidence90(:,:,1,2) =...
    IRF_bootSort(:,:,round(nBoot*0.05),2);
IRF_confidence90(:,:,2,2) =...
    IRF_bootSort(:,:,round(nBoot*0.95),2);


%==========================================================================
% Plot IRFs and Confidence Intervals
%==========================================================================
figure
subplot(2,2,1)
plot(zeros(1,20),'k-')
hold on
plot(IRF(7,:,1),'b-')
plot(IRF_confidence68(7,:,1,1),'b--')
plot(IRF_confidence68(7,:,2,1),'b--')
plot(IRF_confidence90(7,:,1,1),'b-.')
plot(IRF_confidence90(7,:,2,1),'b-.')
hold off
xlim([1,nImp])
xlabel('Quarters')
title('Consumption (nondurables and services)')

subplot(2,2,2)
plot(zeros(1,20),'k-')
hold on
plot(IRF(7,:,2),'b-')
plot(IRF_confidence68(7,:,1,2),'b--')
plot(IRF_confidence68(7,:,2,2),'b--')
plot(IRF_confidence90(7,:,1,2),'b-.')
plot(IRF_confidence90(7,:,2,2),'b-.')
hold off
xlim([1,nImp])
xlabel('Quarters')
title('Consumption (nondurables and services)')

subplot(2,2,3)
plot(zeros(1,20),'k-')
hold on
plot(IRF(8,:,1),'b-')
plot(IRF_confidence68(8,:,1,1),'b--')
plot(IRF_confidence68(8,:,2,1),'b--')
plot(IRF_confidence90(8,:,1,1),'b-.')
plot(IRF_confidence90(8,:,2,1),'b-.')
hold off
xlim([1,nImp])
xlabel('Quarters')
title('Durable good purchases')

subplot(2,2,4)
plot(zeros(1,20),'k-')
hold on
plot(IRF(8,:,2),'b-')
plot(IRF_confidence68(8,:,1,2),'b--')
plot(IRF_confidence68(8,:,2,2),'b--')
plot(IRF_confidence90(8,:,1,2),'b-.')
plot(IRF_confidence90(8,:,2,2),'b-.')
hold off
xlim([1,nImp])
xlabel('Quarters')
title('Durable good purchases')

toc;

