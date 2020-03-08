function VARci = doDeltaMethod(VAR,clevel,NWlags)

[~,What,~] = CovAhat_Sigmahat_Gamma(VAR.p,VAR.X(:,[end end-size(VAR.DET,2)+1:end-1 1:VAR.n*VAR.p]),VAR.m,VAR.res',NWlags);                

AL          = VAR.bet(1:end-size(VAR.DET,2),:)';
%Sigma       = VAR.Sigma;
Gamma       = VAR.res'*VAR.m./VAR.T; 

Caux        = [eye(VAR.n),MARep(AL,VAR.p,VAR.irhor-1)]; 
%The function MARep uses the reduced-form estimator vec(A) to 
%compute the implied MA coefficients. You can replace this function
%by your own routine, but make sure that the dimensions match. 

C        = reshape(Caux,[VAR.n,VAR.n,VAR.irhor]); 
[G,~]    = Gmatrices(AL,Caux,VAR.p,VAR.irhor,VAR.n); 

scale = -1;
e = eye(VAR.n);

for jj = 1:length(clevel)
    critval=norminv(1-((1-clevel(jj)/100)/2),0,1)^2;
            
        lambdahat=zeros(VAR.n,VAR.irhor);
        DmethodVar=zeros(VAR.n,VAR.irhor);
        Dmethodlbound=zeros(VAR.n,VAR.irhor);
        Dmethodubound=zeros(VAR.n,VAR.irhor);

    for ih=1:VAR.irhor
        for ivar=1:VAR.n
            lambdahat(ivar,ih)=scale*e(:,ivar)'*C(:,:,ih)*Gamma./Gamma(1,1);
            d1=(kron(Gamma',e(:,ivar)')*scale*G(:,:,ih));
            d2=(scale*e(:,ivar)'*C(:,:,ih))-(lambdahat(ivar,ih)*e(:,1)');                                    
            d=[d1,d2]';         
            DmethodVar(ivar,ih)=d'*What*d;           
            Dmethodlbound(ivar,ih)= lambdahat(ivar,ih)-...
                ((critval./VAR.T)^.5)*(DmethodVar(ivar,ih)^.5)/abs(Gamma(1,1));
            Dmethodubound(ivar,ih)= lambdahat(ivar,ih)+...
                ((critval./VAR.T)^.5)*(DmethodVar(ivar,ih)^.5)/abs(Gamma(1,1));  
        end    
    end
VARci.irsL(:,:,jj) = Dmethodubound';
VARci.irsH(:,:,jj) = Dmethodlbound';
end

end





