function VARci = doMSWwivrobust(VAR,clevel,NWlags)

[~,What,~] = CovAhat_Sigmahat_Gamma(VAR.p,VAR.X(:,[end end-size(VAR.DET,2)+1:end-1 1:VAR.n*VAR.p]),VAR.m,VAR.res',NWlags);                

AL          = VAR.bet(1:end-size(VAR.DET,2),:)';
%Sigma       = VAR.Sigma;
Gamma       = VAR.res'*VAR.m./VAR.T; 
Caux        = [eye(VAR.n),MARep(AL,VAR.p,VAR.irhor-1)]; 

C        = reshape(Caux,[VAR.n,VAR.n,VAR.irhor]); 
[G,~]    = Gmatrices(AL,Caux,VAR.p,VAR.irhor,VAR.n);

%% Label the submatrices of the asy var of (vecA,Gamma)  
W1          = What(1:(VAR.n^2)*VAR.p,1:(VAR.n^2)*VAR.p);
W12         = What(1:(VAR.n^2)*VAR.p,1+(VAR.n^2)*VAR.p:end);
W2          = What(1+(VAR.n^2)*VAR.p:end,1+(VAR.n^2)*VAR.p:end);

%% 4) Definitions to apply the formulae in MSW for noncumulative IRFs

%a) Definitions to compute the MSW confidence interval for $\lambda_{k,i}$
for jj = 1:length(clevel)
    critval=norminv(1-((1-clevel(jj)/100)/2),0,1)^2;

    e         = eye(VAR.n);
    ahat      = zeros(VAR.n,VAR.irhor); 
    bhat      = zeros(VAR.n,VAR.irhor);
    chat      = zeros(VAR.n,VAR.irhor);
    Deltahat  = zeros(VAR.n,VAR.irhor);
    MSWlbound = zeros(VAR.n,VAR.irhor);
    MSWubound = zeros(VAR.n,VAR.irhor);
    casedummy = zeros(VAR.n,VAR.irhor);

    scale = -1;

    for j =1:VAR.n

        for ih=1:VAR.irhor
        ahat(j,ih)     = (VAR.T*(Gamma(1,1)^2))-(critval*W2(1,1));
        bhat(j,ih)     = -2*VAR.T*scale*(e(:,j)'*C(:,:,ih)*Gamma)*Gamma(1,1)...
            + 2*critval*scale*(kron(Gamma',e(:,j)'))*G(:,:,ih)*W12(:,1)...
            + 2*critval*scale*e(:,j)'*C(:,:,ih)*W2(:,1);
        chat(j,ih)     = ((VAR.T^.5)*scale*e(:,j)'*C(:,:,ih)*Gamma).^2 ...
            -critval*(scale^2)*(kron(Gamma',e(:,j)'))*G(:,:,ih)*W1*...
            ((kron(Gamma',e(:,j)'))*G(:,:,ih))' ...
            -2*critval*(scale^2)*(kron(Gamma',e(:,j)'))*G(:,:,ih)*W12*C(:,:,ih)'*e(:,j)...
            -critval*(scale^2)*e(:,j)'*C(:,:,ih)*W2*C(:,:,ih)'*e(:,j); 
        Deltahat(j,ih) = bhat(j,ih).^2-(4*ahat(j,ih)*chat(j,ih));

        if ahat(j,ih)>0 && Deltahat(j,ih)>0
           casedummy(j,ih) = 1;
           MSWlbound(j,ih) = (-bhat(j,ih) - (Deltahat(j,ih)^.5))/(2*ahat(j,ih));
           MSWubound(j,ih) = (-bhat(j,ih) + (Deltahat(j,ih)^.5))/(2*ahat(j,ih));
        elseif ahat(j,ih)<0 && Deltahat(j,ih)>0
           casedummy(j,ih) = 2;
           MSWlbound(j,ih) = (-bhat(j,ih) + (Deltahat(j,ih)^.5))/(2*ahat(j,ih));
           MSWubound(j,ih) = (-bhat(j,ih) - (Deltahat(j,ih)^.5))/(2*ahat(j,ih));
        elseif ahat(j,ih)>0 && Deltahat(j,ih)<0
           casedummy(j,ih) = 3;
           MSWlbound(j,ih) = NaN;
           MSWubound(j,ih) = NaN;
        else 
           casedummy(j,ih) = 4;
           MSWlbound(j,ih) = -inf;
           MSWubound(j,ih) = inf;
        end

        end
    end
        MSWlbound(1,1)=scale;
        MSWubound(1,1)=scale;
        
VARci.irsL(:,:,jj) = MSWubound';
VARci.irsH(:,:,jj) = MSWlbound';        
end
