function VARci = doMSWbootstrap(VAR,nboot,clevel,NWlags)
[Whatall,~,V] = CovAhat_Sigmahat_Gamma(VAR.p,VAR.X(:,[end end-size(VAR.DET,2)+1:end-1 1:VAR.n*VAR.p]),VAR.m,VAR.res',NWlags);                

AL          = VAR.bet(1:end-size(VAR.DET,2),:)';
Sigma       = VAR.Sigma;
Gamma       = VAR.res'*VAR.m./VAR.T; 

% Make sure that Whatall is symmetric and positive semidefinite
Whatall     = (Whatall + Whatall')/2;     
[aux1,aux2] = eig(Whatall);     
Whatall     = aux1*max(aux2,0)*aux1'; 

% Generate draws from the Gaussian asy. dist. 
% (centered at point estimators)
dall        = size(Whatall,1);
gvar    = [mvnrnd(zeros(nboot,dall),(Whatall)/VAR.T)'];
           %Added an extra column of zeros to access point estimators       
Draws   = bsxfun(@plus,gvar,[AL(:);V*Sigma(:);Gamma(:)]);
 
 
%%  Evaluate the parameter of interest 
% (which is allowed to depend on the full vector vecA, vechSigma, Gamma)
pdSigma    = zeros(1,nboot);
IRFs       = zeros(VAR.n, VAR.irhor, nboot);
    
 for idraws = 1:nboot   
       AL   = reshape(Draws(1:(VAR.n^2)*VAR.p,idraws),[VAR.n,VAR.n*VAR.p]);   % Generate the draws for AL    
  vechSigma = Draws((VAR.n^2)*VAR.p+1:(VAR.n^2)*VAR.p+(VAR.n*(VAR.n+1)/2),idraws); % Generate the draws from Sigma
      Sigma = tril(ones(VAR.n),0); Sigma(Sigma==1) = vechSigma';       
      Sigma = Sigma + tril(Sigma,-1)';      
      %Check if the draws are positive definite
      if min(eig(Sigma))>0
          pdSigma(1,idraws) = 1;      
          Gamma = reshape(Draws(((VAR.n^2)*VAR.p)+(VAR.n*(VAR.n+1)/2)+1:end,idraws),[VAR.n,VAR.k]); %Draws from Gamma                        
           IRS = doIRFs(AL,Sigma,Gamma,VAR,-1);
          for i=1:size(VAR.irs,3) 
           IRFs(:,:,idraws,i) = IRS(:,:,i)';   
          end
      end 
 end
    
%% Implement "Standard" Bootstrap Inference
 
 aux        = reshape(pdSigma,[1,1,nboot]);
  
for jj = 1:length(clevel)
    for i=1:size(VAR.irs,3)
        bootsIRFs  = quantile(IRFs(:,:,aux==1,i),...
                           [((1-clevel(jj)/100)/2),1-((1-clevel(jj)/100)/2)],3); 
        VARci.irsH(:,:,i,jj) = bootsIRFs(:,:,1)';
        VARci.irsL(:,:,i,jj) = bootsIRFs(:,:,2)';
        
        bootsIRFsH  = quantile(IRFs(:,:,aux==1,i)-repmat(VAR.irs(:,:,i)',[1 1 size(IRFs(:,:,aux==1,i),3)]),...
                           [((1-clevel(jj)/100)/2),1-((1-clevel(jj)/100)/2)],3);
                       
        VARci.irsHhall(:,:,i,jj) = VAR.irs(:,:,i)-bootsIRFsH(:,:,1)';
        VARci.irsLhall(:,:,i,jj) = VAR.irs(:,:,i)-bootsIRFsH(:,:,2)';
    end
end

