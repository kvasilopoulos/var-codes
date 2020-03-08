%==========================================================================
% This script runs the Monte Carlo simulations in "Proxy SVARs: Asymptotic 
% Theory, Bootstrap Inference, and the Effects of Income Tax Changes in the 
% United States" by Carsten Jentsch and Kurt G. Lunsford.  
% 
% This script does not call any workspaces.
% 
% This script calls the following functions:
%   IdentifyFunc.m
%==========================================================================

clc
clear;
format short g
tic;

%==========================================================================
% Parameters of the Script
%==========================================================================
%parameters of the data generating process
lags = 1;                       %number of lags in the VAR
A = [0.2,0;0.5,0.5];            %the VAR(1) matrix
HH = [1,0.3;0.3,1];             %variance-covariance matrix of u(t)
Psi = 0.5;                      %covariance of m(t) and epsilon_1(t)
gamma1 = 0.05;                  %GARCH(1,1) parameter
gamma2 = 0.90;                  %GARCH(1,1) parameter
gamma0 = 1 - gamma1 - gamma2;   %GARCH(1,1) parameter

%IRF parameters
nImp = 6;                       %number of steps in the IRF

%simulation parameters
nSample = [100,250,500,1000];       %sample size for each simulation
BlockSize =[16,20,24,28];           %size of the blocks used for each sample size
nBoot = 2000;                       %bootstrap replication
nSim = 1000;                        %number of simulations

%toggle to use iid or GARCH(1,1)
toggleGARCH = 0;            %0 is iid
                            %1 is GARCH(1,1)

%toggle to use one standard deviation or normalized IRFs
toggleShock = 0;            %0 is one standard deviation
                            %1 normalizes the initial impulse to -1

%randon number generator
seed = 1;
rng(seed);


%==========================================================================
% Generate the True Impulse Responses
%==========================================================================
%produce the H matrix
[eVec,eVal] = eig(HH);
H = -eVec*sqrt(eVal);

%number of variables in the VAR
N = size(A,2);

%choose the shock for the IRFs
if toggleShock == 0
    shock = 1;
elseif toggleShock == 1
    shock = -1/H(1,1);
end

%compute true IRF
IRF = zeros(N,nImp);
IRF(:,1) = H(:,1)*shock;
for t = 2:nImp
    IRF(:,t) = A*IRF(:,t-1);
end


%==========================================================================
% Simulation
%==========================================================================
%storage matrices for moving block bootstrap
IRF_include68MBB = zeros(N,nImp,size(nSample,2),nSim);
IRF_include95MBB = zeros(N,nImp,size(nSample,2),nSim);

%storage matrices for Rademacher wild bootstrap
IRF_include68WBR = zeros(N,nImp,size(nSample,2),nSim);
IRF_include95WBR = zeros(N,nImp,size(nSample,2),nSim);

%storage matrices for Normal wild bootstrap
IRF_include68WBN = zeros(N,nImp,size(nSample,2),nSim);
IRF_include95WBN = zeros(N,nImp,size(nSample,2),nSim);

%run the simulation
for sim = 1:nSim
    %----------------------------------------------------------------------
    % Produce the simulated data
    %----------------------------------------------------------------------

    %simulate the structural shocks and VAR innovations
    if toggleGARCH == 0         %iid simulation
        %structural shocks
        Eps = randn(1000+nSample(1,end),N);
        Eps1 = Eps(:,1);
        
        %VAR innovations
        U = Eps*H';
    elseif toggleGARCH == 1;    %GARCH(1,1) simulation
        %structural shocks
        w = randn(1000+nSample(1,end),N);
        g = zeros(1000+nSample(1,end),N);
        g(1,:) = ones(1,N);
        Eps = zeros(1000+nSample(1,end),N);
        Eps(1,:) = g(1,:).*w(1,:);
        for t = 2:1000+nSample(1,end)
            g(t,:) = sqrt(gamma0*ones(1,N) +...
                gamma1*Eps(t-1,:).*Eps(t-1,:) + gamma2*g(t-1,:).*g(t-1,:));
            Eps(t,:) = g(t,:).*w(t,:);
        end
        Eps1 = Eps(:,1);
        
        %VAR innovations
        U = Eps*H';
    end

    %simulate the time series Y(t)
    Y = zeros(1000+nSample(1,end),N);
    Y(1,:) = U(1,:);
    for t = 2:1000+nSample(1,end)
        Y(t,:) = Y(t-1,:)*A' + U(t,:);
    end

    %drop the first 1000 observations
    Eps1 = Eps1(1001:1000+nSample(1,end),:);
    X = Y(1000:999+nSample(1,end),:);
    Y = Y(1001:1000+nSample(1,end),:);

    %construct the proxy variable
    V = randn(nSample(1,end),1);
    M = Psi*Eps1 + V;

    %----------------------------------------------------------------------
    % Loop for each sample size
    %----------------------------------------------------------------------
    for loop = 1:size(nSample,2)
        %set the sample size and the sample
        T = nSample(1,loop);
        Xsample = X(1:T,:);
        Ysample = Y(1:T,:);
        Msample = M(1:T,:);

        %estimate the model
        RHS = [Xsample,ones(T,1)];
        A_hat = (Ysample'*RHS)/(RHS'*RHS);
        U_hat = Ysample - RHS*A_hat';

        %------------------------------------------------------------------
        % The Moving Block Bootstrap
        %------------------------------------------------------------------
        %create the blocks for the moving block bootstrap
        L = BlockSize(1,loop);
        J = ceil(T/L);
        UBlocks = zeros(L,N,T-L+1);
        MBlocks = zeros(L,1,T-L+1);
        for ii = 1:T-L+1
            UBlocks(:,:,ii) = U_hat(ii:L+ii-1,:);
            MBlocks(:,:,ii) = M(ii:L+ii-1,:);
        end

        %center the bootstrapped VAR errors and proxy variables
        centering = zeros(L,N);
        Mcentering = zeros(L,1);
        for ii = 1:L
            centering(ii,:) = mean(U_hat(ii:T-L+ii,:),1);
            Mcentering(ii,:) = mean(M(ii:T-L+ii,:),1);
        end
        centering = repmat(centering,[J,1]);
        centering = centering(1:T,:);
        Mcentering = repmat(Mcentering,[J,1]);
        Mcentering = Mcentering(1:T,:);

        %the bootstrap iterations
        IRF_boot = zeros(N,nImp,nBoot);
        for boot = 1:nBoot
            %draw bootstrapped VAR residuals and proxies
            index = ceil((T - L + 1)*rand(J,1));
            U_boot = zeros(J*L,N);
            M_boot = zeros(J*L,1);
            for j = 1:J
                U_boot(1+L*(j-1):L*j,:) = UBlocks(:,:,index(j,1));
                M_boot(1+L*(j-1):L*j,:) = MBlocks(:,:,index(j,1));
            end
            U_boot = U_boot(1:T,:);
            M_boot = M_boot(1:T,:);

            %center the VAR residuals and proxies
            U_boot = U_boot - centering;
            M_boot = M_boot - Mcentering;

            %recursively generate bootstrapped data
            RHS_boot = [[X(1,:),1];[zeros(T-1,N),ones(T-1,1)]];
            Y_boot = zeros(T,N);
            for t = 1:T-1
                Y_boot(t,:) = RHS_boot(t,:)*A_hat' + U_boot(t,:);
                RHS_boot(t+1,1:N) = Y_boot(t,:);
            end
            Y_boot(T,:) = RHS_boot(T,:)*A_hat' + U_boot(T,:);

            %regression
            A_hatboot = (Y_boot'*RHS_boot)/(RHS_boot'*RHS_boot);
            U_hatboot = Y_boot - RHS_boot*A_hatboot';

            %identification
            H1_boot = IdentifyFunc(U_hatboot,M_boot);

            %choose the shock for the IRFs
            if toggleShock == 0
                shock_boot = 1;
            elseif toggleShock == 1
                shock_boot = -1/H1_boot(1,1);
            end

            %estimate the IRFs
            IRF_boot(:,1,boot) = H1_boot*shock_boot;
            for t = 2:nImp
                IRF_boot(:,t,boot) = A_hatboot(:,1:N)*IRF_boot(:,t-1,boot);
            end
        end

        %confidence intervals for IRFs
        IRF_bootSort = zeros(N,nImp,nBoot);
        for n = 1:N
            for i = 1:nImp
                IRF_bootSort(n,i,:) = sort(IRF_boot(n,i,:));
            end
        end
        IRF_confidence68(:,:,1) = IRF_bootSort(:,:,nBoot*0.16);
        IRF_confidence68(:,:,2) = IRF_bootSort(:,:,nBoot*0.84);
        IRF_include68MBB(:,:,loop,sim) = (IRF_confidence68(:,:,2) >= IRF).*...
            (IRF >= IRF_confidence68(:,:,1));

        IRF_confidence95(:,:,1) = IRF_bootSort(:,:,nBoot*0.025);
        IRF_confidence95(:,:,2) = IRF_bootSort(:,:,nBoot*0.975);
        IRF_include95MBB(:,:,loop,sim) = (IRF_confidence95(:,:,2) >= IRF).*...
            (IRF >= IRF_confidence95(:,:,1));

        %------------------------------------------------------------------
        % The Wild Bootstrap with Rademacher Multipliers
        %------------------------------------------------------------------
        %the bootstrap iterations
        IRF_boot = zeros(N,nImp,nBoot);
        for boot = 1:nBoot
            %generate the shocks for the wild bootstrap
            E = 2*(rand(T,1) > 0.5) - 1;

            %compute the bootstrapped shocks and proxy variable
            U_boot = U_hat.*(E*ones(1,N));
            M_boot = Msample.*E;

            %recursively generate bootstrapped data
            RHS_boot = [[X(1,:),1];[zeros(T-1,N),ones(T-1,1)]];
            Y_boot = zeros(T,N);
            for t = 1:T-1
                Y_boot(t,:) = RHS_boot(t,:)*A_hat' + U_boot(t,:);
                RHS_boot(t+1,1:N) = Y_boot(t,:);
            end
            Y_boot(T,:) = RHS_boot(T,:)*A_hat' + U_boot(T,:);

            %regression
            A_hatboot = (Y_boot'*RHS_boot)/(RHS_boot'*RHS_boot);
            U_hatboot = Y_boot - RHS_boot*A_hatboot';

            %identification
            H1_boot = IdentifyFunc(U_hatboot,M_boot);

            %choose the shock for the IRFs
            if toggleShock == 0
                shock_boot = 1;
            elseif toggleShock == 1
                shock_boot = -1/H1_boot(1,1);
            end

            %estimate the IRFs
            IRF_boot(:,1,boot) = H1_boot*shock_boot;
            for t = 2:nImp
                IRF_boot(:,t,boot) = A_hatboot(:,1:N)*IRF_boot(:,t-1,boot);
            end
        end

        %confidence intervals for IRFs
        IRF_bootSort = zeros(N,nImp,nBoot);
        for n = 1:N
            for i = 1:nImp
                IRF_bootSort(n,i,:) = sort(IRF_boot(n,i,:));
            end
        end
        IRF_confidence68(:,:,1) = IRF_bootSort(:,:,nBoot*0.16);
        IRF_confidence68(:,:,2) = IRF_bootSort(:,:,nBoot*0.84);
        IRF_include68WBR(:,:,loop,sim) = (IRF_confidence68(:,:,2) >= IRF).*...
            (IRF >= IRF_confidence68(:,:,1));

        IRF_confidence95(:,:,1) = IRF_bootSort(:,:,nBoot*0.025);
        IRF_confidence95(:,:,2) = IRF_bootSort(:,:,nBoot*0.975);
        IRF_include95WBR(:,:,loop,sim) = (IRF_confidence95(:,:,2) >= IRF).*...
            (IRF >= IRF_confidence95(:,:,1));

        %------------------------------------------------------------------
        % The Wild Bootstrap with Standard Normal Multipliers
        %------------------------------------------------------------------
        %the bootstrap iterations
        IRF_boot = zeros(N,nImp,nBoot);
        for boot = 1:nBoot
            %generate the shocks for the wild bootstrap
            E = randn(T,1);

            %compute the bootstrapped shocks and proxy variable
            U_boot = U_hat.*(E*ones(1,N));
            M_boot = Msample.*E;

            %recursively generate bootstrapped data
            RHS_boot = [[X(1,:),1];[zeros(T-1,N),ones(T-1,1)]];
            Y_boot = zeros(T,N);
            for t = 1:T-1
                Y_boot(t,:) = RHS_boot(t,:)*A_hat' + U_boot(t,:);
                RHS_boot(t+1,1:N) = Y_boot(t,:);
            end
            Y_boot(T,:) = RHS_boot(T,:)*A_hat' + U_boot(T,:);

            %regression
            A_hatboot = (Y_boot'*RHS_boot)/(RHS_boot'*RHS_boot);
            U_hatboot = Y_boot - RHS_boot*A_hatboot';

            %identification
            H1_boot = IdentifyFunc(U_hatboot,M_boot);

            %choose the shock for the IRFs
            if toggleShock == 0
                shock_boot = 1;
            elseif toggleShock == 1
                shock_boot = -1/H1_boot(1,1);
            end

            %estimate the IRFs
            IRF_boot(:,1,boot) = H1_boot*shock_boot;
            for t = 2:nImp
                IRF_boot(:,t,boot) = A_hatboot(:,1:N)*IRF_boot(:,t-1,boot);
            end
        end

        %confidence intervals for IRFs
        IRF_bootSort = zeros(N,nImp,nBoot);
        for n = 1:N
            for i = 1:nImp
                IRF_bootSort(n,i,:) = sort(IRF_boot(n,i,:));
            end
        end
        IRF_confidence68(:,:,1) = IRF_bootSort(:,:,nBoot*0.16);
        IRF_confidence68(:,:,2) = IRF_bootSort(:,:,nBoot*0.84);
        IRF_include68WBN(:,:,loop,sim) = (IRF_confidence68(:,:,2) >= IRF).*...
            (IRF >= IRF_confidence68(:,:,1));

        IRF_confidence95(:,:,1) = IRF_bootSort(:,:,nBoot*0.025);
        IRF_confidence95(:,:,2) = IRF_bootSort(:,:,nBoot*0.975);
        IRF_include95WBN(:,:,loop,sim) = (IRF_confidence95(:,:,2) >= IRF).*...
            (IRF >= IRF_confidence95(:,:,1));

        display(loop)
    end

display(sim)
end

%storage arrays for coverage rates
IRF_Coverage68MBB = zeros(nImp,size(nSample,2)*N);
IRF_Coverage95MBB = zeros(nImp,size(nSample,2)*N);
IRF_Coverage68WBR = zeros(nImp,size(nSample,2)*N);
IRF_Coverage95WBR = zeros(nImp,size(nSample,2)*N);
IRF_Coverage68WBN = zeros(nImp,size(nSample,2)*N);
IRF_Coverage95WBN = zeros(nImp,size(nSample,2)*N);

for loop = 1:size(nSample,2)
    %coverage rates for moving block bootstrap
    IRF_Coverage68MBB(:,1+N*(loop-1):N*loop) =...
        mean(IRF_include68MBB(:,:,loop,:),4)';
    IRF_Coverage95MBB(:,1+N*(loop-1):N*loop) =...
        mean(IRF_include95MBB(:,:,loop,:),4)';

    %coverage rates for wild bootstrap with (-1,1) shocks
    IRF_Coverage68WBR(:,1+N*(loop-1):N*loop) =...
        mean(IRF_include68WBR(:,:,loop,:),4)';
    IRF_Coverage95WBR(:,1+N*(loop-1):N*loop) =...
        mean(IRF_include95WBR(:,:,loop,:),4)';

    %coverage rates for wild bootstrap with standard normal shocks
    IRF_Coverage68WBN(:,1+N*(loop-1):N*loop) =...
        mean(IRF_include68WBN(:,:,loop,:),4)';
    IRF_Coverage95WBN(:,1+N*(loop-1):N*loop) =...
        mean(IRF_include95WBN(:,:,loop,:),4)';
end

toc;

