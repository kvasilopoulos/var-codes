%==========================================================================
% This script produces Figures 1 and 2 of "Proxy SVARs: Asymptotic Theory, 
% Bootstrap Inference, and the Effects of Income Tax Changes in the United 
% States" by Carsten Jentsch and Kurt G. Lunsford.  
% 
% This script calls the following workspaces:
%   VARData.mat 
%   ProxyData.mat
% 
% This script call the following functions:
%   IdentifyFunc.m
%   acf.m
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
alpha = 0.68;           %level of significance for confidence intervals

%bootstrap parameters
nBoot = 500;          %number of bootstrap replications
BlockSize = 19;         %size of blocks in the bootstrap
seed = 1;               %seed for random number generator
rng(seed);              %iniate the random number generator

%program toggles
toggleShock = 2;        %1 to produce the shock to the APITR
                        %2 to produce the shock to the ACITR
toggleBoot = 1;         %1 to use moving block bootstrap
                        %2 to use wild bootstrap


%==========================================================================
% Load Data
%==========================================================================
%data for benchmark specification in Mertens and Ravn (2013)
load VARData.mat

%dimensions of the VAR data
[T,K] = size(VARData);

%re-order variables depending on shock of interest
if toggleShock == 2
    VARData(:,:) = [VARData(:,2),VARData(:,1),VARData(:,3:K)];
end

%store VAR data for both orderings
VARData(:,:,2) = [VARData(:,2),VARData(:,1),VARData(:,3:K)];

%tax proxy variables used in Mertens and Ravn (2013)
load ProxyData.mat
M = ProxyData(1+lags:T,:);  %adjust the proxy data for the number of lags
r = size(M,2);              %number of proxy variables


%==========================================================================
% Estimate the VAR
%==========================================================================
%storage arrays
LHS = zeros(T-lags,K,2);
RHS = zeros(T-lags,1+K*lags,2);
A_hat = zeros(1+K*lags,K,2);
U_hat = zeros(T-lags,K,2);

%estimate VARs for both APITR and ACITR ordered first
for order = 1:2
    %left-hand side variables
    LHS(:,:,order) = VARData(1+lags:T,:,order);

    %right-hand side variables
    preRHS = ones(T-lags,1);
    for j = 1:lags
        preRHS = [preRHS,VARData(1+lags-j:T-j,:,order)];
    end
    RHS(:,:,order) = preRHS;

    %VAR coefficients and innovations
    A_hat(:,:,order) =...
        (RHS(:,:,order)'*RHS(:,:,order))\(RHS(:,:,order)'*LHS(:,:,order));
    U_hat(:,:,order) = LHS(:,:,order) - RHS(:,:,order)*A_hat(:,:,order);
end


%==========================================================================
% Autocorrelation Functions of VAR Residuals
%==========================================================================
%acf of u(t)
rhoU = acf(U_hat,7);

%acf of |u(t)|
rhoUabs = acf(abs(U_hat),7);

%acf of u(t)^2
rhoUsq = acf((U_hat.^2),7);

%table of autocorrelations for appendix
Table = [rhoU;rhoUabs;rhoUsq];


%==========================================================================
% Identify the Structural Shocks
%==========================================================================
%construct B1 for both the APITR and ACITR ordered first
H1_hat = zeros(K,r,2);
for order = 1:2
    H1_hat(:,:,order) = IdentifyFunc(U_hat(:,:,order),M);
end

%==========================================================================
% Produce the IRFs
%==========================================================================
%impulse response functions for both the APITR and ACITR ordered first
IRF = zeros(K+2,nImp,2);
for order = 1:2
    %normalize the shock to -1
    shock = -1/H1_hat(order,order,order);

    %generate the IRFs
    IRF(1:K,1,order) = H1_hat(:,order,order)*shock;
    history = [IRF(1:K,1,order);zeros((lags-1)*K,1)];
    for j = 2:nImp
        IRF(1:K,j,order) = A_hat(2:K*lags+1,:,order)'*history;
        history = [IRF(1:K,j,order);history(1:(lags-1)*K,1)];
    end

    %personal income tax revenue
    IRF(K+1,:,order) =...
        IRF(order,:,order)/0.1667 + IRF(3,:,order);

    %corporate income tax revenue
    IRF(K+2,:,order) =...
        IRF(order,:,order)/0.2996 + IRF(4,:,order);
end


%==========================================================================
% Bootstrap and Confidence Intervals
%==========================================================================
%set the number of blocks
if toggleBoot == 1
    nBlock = ceil((T-lags)/BlockSize);
end

%arrays to store relevant variables
H1_hatboot = zeros(K,r,nBoot,2);
IRF_boot = zeros(K+2,nImp,nBoot,2);
M_count = zeros(nBoot,r,2);

%loop for ordering of variables
for order = 1:2
    %if using the moving block bootstrap, create the blocks and centerings
    if toggleBoot == 1
        Blocks = zeros(BlockSize,K,T-lags-BlockSize+1);
        MBlocks = zeros(BlockSize,r,T-lags-BlockSize+1);
        for j = 1:T-lags-BlockSize+1
            Blocks(:,:,j) = U_hat(j:BlockSize+j-1,:,order);
            MBlocks(:,:,j) = M(j:BlockSize+j-1,:);
        end

        %center the bootstrapped VAR errors
        centering = zeros(BlockSize,K);
        for j = 1:BlockSize
            centering(j,:) = mean(U_hat(j:T-lags-BlockSize+j,:,order),1);
        end
        centering = repmat(centering,[nBlock,1]);
        centering = centering(1:T-lags,:);

        %center the bootstrapped proxy variables
        Mcentering = zeros(BlockSize,r);
        for j = 1:BlockSize
            subM = M(j:T-lags-BlockSize+j,:);
            Mcentering(j,:) = [mean(subM((subM(:,1) ~= 0),1),1),...
                mean(subM((subM(:,2) ~= 0),2),1)];
        end
        Mcentering = repmat(Mcentering,[nBlock,1]);
        Mcentering = Mcentering(1:T-lags,:);
    end

    %loop for the bootstrap replications
    for boot = 1:nBoot
        if toggleBoot == 1
            %draw bootstrapped residuals and proxies
            index = ceil((T - lags - BlockSize + 1)*rand(nBlock,1));
            U_boot = zeros(nBlock*BlockSize,K);
            M_boot = zeros(nBlock*BlockSize,r);
            for j = 1:nBlock
                U_boot(1+BlockSize*(j-1):BlockSize*j,:) = Blocks(:,:,index(j,1));
                M_boot(1+BlockSize*(j-1):BlockSize*j,:) = MBlocks(:,:,index(j,1));
            end
            U_boot = U_boot(1:T-lags,:);
            M_boot = M_boot(1:T-lags,:);

            %center the bootstrapped residuals and proxies
            U_boot = U_boot - centering;
            for j = 1:r
                M_boot((M_boot(:,j)~=0),j) =...
                    M_boot((M_boot(:,j)~=0),j) - Mcentering((M_boot(:,j)~=0),j);
            end

            %produce the bootstrapped left- and right-hand side data
            LHS_boot = zeros(T-lags,K);
            RHS_boot = [ones(T-lags,1),zeros(T-lags,K*lags)];
            RHS_boot(1,:) = RHS(1,:,order);
            for j = 1:T-lags-1
                LHS_boot(j,:) = RHS_boot(j,:)*A_hat(:,:,order) + U_boot(j,:);
                RHS_boot(j+1,:) = [1,LHS_boot(j,:),RHS_boot(j,2:K*(lags-1)+1)];
            end
            LHS_boot(T-lags,:) =...
                RHS_boot(T-lags,:)*A_hat(:,:,order) + U_boot(T-lags,:);
        elseif toggleBoot == 2
            %generate the shocks for the wild bootstrap
            E = 2*(rand(T-lags,1) > 0.5) - 1;    

            %bootstrapped proxies
            M_boot = M.*(E*ones(1,r));

            %compute the bootstrapped residuals
            U_boot = U_hat(:,:,order).*(E*ones(1,K));

            %produce the bootstrapped left- and right-hand side data
            LHS_boot = zeros(T-lags,K);
            RHS_boot = [ones(T-lags,1),zeros(T-lags,K*lags)];
            RHS_boot(1,:) = RHS(1,:,order);
            for j = 1:T-lags
                LHS_boot(j,:) = RHS_boot(j,:)*A_hat(:,:,order) + U_boot(j,:);
                RHS_boot(j+1,:) = [1,LHS_boot(j,:),RHS_boot(j,2:K*(lags-1)+1)];
            end
            RHS_boot = RHS_boot(1:end-1,:);
        end

        %bootstrap VAR coefficients and innovations
        A_hatboot = (RHS_boot'*RHS_boot)\(RHS_boot'*LHS_boot);
        U_hatboot = LHS_boot - RHS_boot*A_hatboot;

        %identification
        H1_hatboot(:,:,boot,order) = IdentifyFunc(U_hatboot,M_boot);

        %normalize the initial IRF shock to -1
        shock_boot = -1/H1_hatboot(order,order,boot,order);

        %generate the IRFs
        IRF_boot(1:K,1,boot,order) =...
            H1_hatboot(:,order,boot,order)*shock_boot;
        history_boot = [IRF_boot(1:K,1,boot,order);zeros((lags-1)*K,1)];
        for j = 2:nImp
            IRF_boot(1:K,j,boot,order) =...
                A_hatboot(2:K*lags+1,:)'*history_boot;
            history_boot =...
                [IRF_boot(1:K,j,boot,order);history_boot(1:(lags-1)*K,1)];
        end
        
        %personal income tax revenue
        IRF_boot(K+1,:,boot,order) =...
            IRF_boot(order,:,boot,order)/0.1667 + IRF_boot(3,:,boot,order);

        %corporate income tax revenue
        IRF_boot(K+2,:,boot,order) =...
            IRF_boot(order,:,boot,order)/0.2996 + IRF_boot(4,:,boot,order);
        
        %count the number proxy variables not censored to zero
        M_count(boot,:,order) = sum(abs(M_boot) > 0,1);
    end
end

%sort the bootstrapped IRFs
IRF_bootSort = zeros(K+2,nImp,nBoot,2);
for order = 1:2
    for n = 1:K+2
        for j = 1:nImp
            IRF_bootSort(n,j,:,order) = sort(IRF_boot(n,j,:,order));
        end
    end
end

%percentil confidence interval
IRF_confidence = zeros(K+2,nImp,2,2);
IRF_confidence(:,:,1,1) =...
    IRF_bootSort(:,:,round(nBoot*(1 - alpha)/2),1);
IRF_confidence(:,:,2,1) =...
    IRF_bootSort(:,:,round(nBoot*(1 - (1 - alpha)/2)),1);
IRF_confidence(:,:,1,2) =...
    IRF_bootSort(:,:,round(nBoot*(1 - alpha)/2),2);
IRF_confidence(:,:,2,2) =...
    IRF_bootSort(:,:,round(nBoot*(1 - (1 - alpha)/2)),2);

%display the fewest non-censored proxy observations in the bootstrap
[min(M_count(:,:,1));min(M_count(:,:,2))]


%==========================================================================
% Plot IRFs and Confidence Intervals
%==========================================================================
figure
if toggleShock == 1      %plots for shock to APITR
    subplot(3,2,1)
    plot(zeros(1,20),'k-')
    hold on
    plot(IRF(1,:,1),'b-')
    plot(IRF(2,:,2),'r-d')
    plot(IRF_confidence(1,:,1,1),'b--')
    plot(IRF_confidence(1,:,2,1),'b--')
    plot(IRF_confidence(2,:,1,2),'r-.')
    plot(IRF_confidence(2,:,2,2),'r-.')
    hold off
    xlim([1,nImp])
    xlabel('Quarters')
    title('Average Personal Income Tax Rate')

    subplot(3,2,2)
    plot(zeros(1,20),'k-')
    hold on
    plot(IRF(6,:,1),'b-')
    plot(IRF(6,:,2),'r-d')
    plot(IRF_confidence(6,:,1,1),'b--')
    plot(IRF_confidence(6,:,2,1),'b--')
    plot(IRF_confidence(6,:,1,2),'r-.')
    plot(IRF_confidence(6,:,2,2),'r-.')
    hold off
    xlim([1,nImp])
    xlabel('Quarters')
    title('Output')

    subplot(3,2,3)
    plot(zeros(1,20),'k-')
    hold on
    plot(IRF(3,:,1),'b-')
    plot(IRF(3,:,2),'r-d')
    plot(IRF_confidence(3,:,1,1),'b--')
    plot(IRF_confidence(3,:,2,1),'b--')
    plot(IRF_confidence(3,:,1,2),'r-.')
    plot(IRF_confidence(3,:,2,2),'r-.')
    hold off
    xlim([1,nImp])
    xlabel('Quarters')
    title('Personal Income Tax Base')

    subplot(3,2,4)
    plot(zeros(1,20),'k-')
    hold on
    plot(IRF(8,:,1),'b-')
    plot(IRF(8,:,2),'r-d')
    plot(IRF_confidence(8,:,1,1),'b--')
    plot(IRF_confidence(8,:,2,1),'b--')
    plot(IRF_confidence(8,:,1,2),'r-.')
    plot(IRF_confidence(8,:,2,2),'r-.')
    hold off
    xlim([1,nImp])
    xlabel('Quarters')
    title('Personal Income Tax Revenues')

    subplot(3,2,5)
    plot(zeros(1,20),'k-')
    hold on
    plot(IRF(2,:,1),'b-')
    plot(IRF(1,:,2),'r-d')
    plot(IRF_confidence(2,:,1,1),'b--')
    plot(IRF_confidence(2,:,2,1),'b--')
    plot(IRF_confidence(1,:,1,2),'r-.')
    plot(IRF_confidence(1,:,2,2),'r-.')
    hold off
    xlim([1,nImp])
    xlabel('Quarters')
    title('Average Corporate Income Tax Rate')

    subplot(3,2,6)
    plot(zeros(1,20),'k-')
    hold on
    plot(IRF(5,:,1),'b-')
    plot(IRF(5,:,2),'r-d')
    plot(IRF_confidence(5,:,1,1),'b--')
    plot(IRF_confidence(5,:,2,1),'b--')
    plot(IRF_confidence(5,:,1,2),'r-.')
    plot(IRF_confidence(5,:,2,2),'r-.')
    hold off
    xlim([1,nImp])
    xlabel('Quarters')
    title('Government Purchases')
elseif toggleShock == 2      %plots for shock to ACITR
    subplot(3,2,1)
    plot(zeros(1,20),'k-')
    hold on
    plot(IRF(2,:,2),'b-')
    plot(IRF(1,:,1),'r-d')
    plot(IRF_confidence(2,:,1,2),'b--')
    plot(IRF_confidence(2,:,2,2),'b--')
    plot(IRF_confidence(1,:,1,1),'r-.')
    plot(IRF_confidence(1,:,2,1),'r-.')
    hold off
    xlim([1,nImp])
    xlabel('Quarters')
    title('Average Corporate Income Tax Rate')

    subplot(3,2,2)
    plot(zeros(1,20),'k-')
    hold on
    plot(IRF(6,:,2),'b-')
    plot(IRF(6,:,1),'r-d')
    plot(IRF_confidence(6,:,1,2),'b--')
    plot(IRF_confidence(6,:,2,2),'b--')
    plot(IRF_confidence(6,:,1,1),'r-.')
    plot(IRF_confidence(6,:,2,1),'r-.')
    hold off
    xlim([1,nImp])
    xlabel('Quarters')
    title('Output')

    subplot(3,2,3)
    plot(zeros(1,20),'k-')
    hold on
    plot(IRF(4,:,2),'b-')
    plot(IRF(4,:,1),'r-d')
    plot(IRF_confidence(4,:,1,2),'b--')
    plot(IRF_confidence(4,:,2,2),'b--')
    plot(IRF_confidence(4,:,1,1),'r-.')
    plot(IRF_confidence(4,:,2,1),'r-.')
    hold off
    xlim([1,nImp])
    xlabel('Quarters')
    title('Corporate Income Tax Base')

    subplot(3,2,4)
    plot(zeros(1,20),'k-')
    hold on
    plot(IRF(9,:,2),'b-')
    plot(IRF(9,:,1),'r-d')
    plot(IRF_confidence(9,:,1,2),'b--')
    plot(IRF_confidence(9,:,2,2),'b--')
    plot(IRF_confidence(9,:,1,1),'r-.')
    plot(IRF_confidence(9,:,2,1),'r-.')
    hold off
    xlim([1,nImp])
    xlabel('Quarters')
    title('Corporate Income Tax Revenues')

    subplot(3,2,5)
    plot(zeros(1,20),'k-')
    hold on
    plot(IRF(1,:,2),'b-')
    plot(IRF(2,:,1),'r-d')
    plot(IRF_confidence(1,:,1,2),'b--')
    plot(IRF_confidence(1,:,2,2),'b--')
    plot(IRF_confidence(2,:,1,1),'r-.')
    plot(IRF_confidence(2,:,2,1),'r-.')
    hold off
    xlim([1,nImp])
    xlabel('Quarters')
    title('Average Personal Income Tax Rate')

    subplot(3,2,6)
    plot(zeros(1,20),'k-')
    hold on
    plot(IRF(5,:,2),'b-')
    plot(IRF(5,:,1),'r-d')
    plot(IRF_confidence(5,:,1,2),'b--')
    plot(IRF_confidence(5,:,2,2),'b--')
    plot(IRF_confidence(5,:,1,1),'r-.')
    plot(IRF_confidence(5,:,2,1),'r-.')
    hold off
    xlim([1,nImp])
    xlabel('Quarters')
    title('Government Purchases')
end

toc;

