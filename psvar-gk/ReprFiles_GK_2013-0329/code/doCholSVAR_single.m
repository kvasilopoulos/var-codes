function VAR = doCholSVAR_single(VAR)
%Cholesky factorization

 X      = lagmatrix(VAR.vars,1:VAR.p);
 X      = X(VAR.p+1:end,:);
 Y      = VAR.vars(VAR.p+1:end,:);
  
% Run VAR
%%%%%%%%%%%%
VAR.bet=[X ones(length(X),1)]\Y; 
VAR.res = Y-[X ones(length(X),1)]*VAR.bet;
VAR.Sigma = (VAR.res'*VAR.res)/(VAR.T-VAR.n*VAR.p-1);

% Identification
%%%%%%%%%%%%%%%%%

VAR.B   =   chol(VAR.Sigma,'lower');

% Impulse Responses
%%%%%%%%%%%%%%%%%%%%
% initial shock: eps(1,1)=1
irs(VAR.p+1,:) = VAR.B(:,3);
 for jj=2:VAR.irhor
    lvars = (irs(VAR.p+jj-1:-1:jj,:))';
    irs(VAR.p+jj,:) = lvars(:)'*VAR.bet(1:VAR.p*VAR.n,:);     
 end
VAR.irs = irs(VAR.p+1:end,:); 