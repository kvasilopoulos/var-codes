function VAR = doPVAR(VAR)

 VARSLAGS = lagmatrix(VAR.vars,1:VAR.p);
 VARSLAGS = VARSLAGS(VAR.p+1:end,:);
 VARS     = VAR.vars(VAR.p+1:end,:);
 MSHOCKS = VAR.mshocks(VAR.p+1:end,:);

 [VAR.T,VAR.n]    = size(VARS);
  VAR.k =size(MSHOCKS,2);
 
% Run Reduced-form VAR
%%%%%%%%%
 VAR.bet=[VARSLAGS VAR.DET(VAR.p+1:end,:)]\VARS; 
 res = VARS-[VARSLAGS VAR.DET(VAR.p+1:end,:)]*VAR.bet;
 VAR.Sigma = (res'*res)/(VAR.T-VAR.n*VAR.p-1);
      
% Narrative Identification
%%%%%%%%%%%%%%%%%%%%%%%%%%
LAMb = MSHOCKS\res;
LAMb11  = LAMb(1:VAR.k,1:VAR.k);
LAMb21  = LAMb(1:VAR.k,VAR.k+1:VAR.n);
b21ib11 = (LAMb11\LAMb21)';

Sig11   = VAR.Sigma(1:VAR.k,1:VAR.k);
Sig21   = VAR.Sigma(VAR.k+1:VAR.n,1:VAR.k);
Sig22   = VAR.Sigma(VAR.k+1:VAR.n,VAR.k+1:VAR.n);
ZZp     = b21ib11*Sig11*b21ib11'-(Sig21*b21ib11'+b21ib11*Sig21')+Sig22;
b12b12p = (Sig21- b21ib11*Sig11)'*(ZZp\(Sig21- b21ib11*Sig11));
b11b11p = Sig11-b12b12p;
b22b22p = ZZp+(Sig21- b21ib11*Sig11)*b21ib11'+b21ib11*(Sig21- b21ib11*Sig11)'+...
    b21ib11*b12b12p*b21ib11';

b22     = chol(b22b22p,'lower');
b12     = ((Sig21- b21ib11*Sig11)'+b12b12p*b21ib11')/(b22');
b12ib22 = b12/b22;

VAR.sigmaG  = b22(1,1);
VAR.thetaG  = b12ib22(:,1);
VAR.thetaY  = b12ib22(:,2);
VAR.thetaW  = b12ib22(:,3);
b11iSig = eye(VAR.k)/(eye(VAR.k)-b12ib22*b21ib11);
b21iSig = b21ib11*b11iSig;

VAR.gammaT  = b21iSig(1,:)';

F2     = (b22(2,1)/VAR.sigmaG-b21iSig(2,:)*VAR.thetaG);
F1     = (-F2*VAR.gammaT'+(1-VAR.gammaT'*VAR.thetaG)*b21iSig(2,:))';
VAR.zetaT = (eye(VAR.k)*(1-VAR.gammaT'*VAR.thetaG)+F1*VAR.thetaY')\F1;
VAR.zetaG = (1-VAR.zetaT'*VAR.thetaY)/(1-VAR.gammaT'*VAR.thetaG)*F2;
VAR.sigmaY= (1-VAR.zetaT'*VAR.thetaY)*b22(2,2);
SigmaTSigmaTp =b11iSig\b11b11p/b11iSig';
     
% Impulse Response to a monetary Shock
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:VAR.k
   ind = 1:VAR.k; ind(VAR.k) =j; ind(j) =VAR.k;
SSp = SigmaTSigmaTp(ind,ind);
SigmaT=chol(SSp,'lower');
VAR.SigmaT(:,:,j) = SigmaT(ind,ind);
VAR.D(:,:,j) = [b11iSig*VAR.SigmaT(:,:,j) b12;
         b21iSig*VAR.SigmaT(:,:,j) b22];
     
irs(VAR.p+1,:) = VAR.D(:,j,j)/VAR.D(j,j,j)*VAR.mshocksize(j);
 for jj=2:VAR.irhor
 lvars = (irs(VAR.p+jj-1:-1:jj,:))';
 irs(VAR.p+jj,:) = lvars(:)'*VAR.bet(1:VAR.p*VAR.n,:);     
 end
 
VAR.irs(:,:,j) = irs(VAR.p+1:end,:); 
if VAR.k==1
VAR.irsTRY = VAR.irs(:,1)-VAR.irs(:,3);
end
end


% Impulse Response to a Spending Shock
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
irsg(VAR.p+1,:) = 100*VAR.D(:,VAR.k+1,1)/VAR.D(VAR.k+1,VAR.k+1,1)*VAR.gshocksize;

 for jj=2:VAR.irhor
 lvarsg = (irsg(VAR.p+jj-1:-1:jj,:))';
 irsg(VAR.p+jj,:) = lvarsg(:)'*VAR.bet(1:VAR.p*VAR.n,:);     
 end
 
VAR.irsg(:,:) = irsg(VAR.p+1:end,:); 


e             = (VAR.D(:,:,1)\res')';
VAR.et        =  e(:,1:VAR.k);
VAR.e         =  e(:,VAR.k+1:VAR.n);   
VAR.b11eT     = (VAR.D(1:VAR.k,1:VAR.k,1)*e(:,1:VAR.k,1)')';
VAR.b21ib11   = b21ib11;
VAR.MSHOCKS = MSHOCKS;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reliability
%%%%%%%%%%%%%%

m = MSHOCKS;
ED    = eye(VAR.k)*sum(sum(m,2)~=0)/VAR.T;
mu1   = m'*res(:,1:VAR.k)/VAR.T;
b11 = sqrt(b11b11p);

Bi = [1/b11+(Sig21-Sig11*b21ib11)'*inv(ZZp)/b11*b21ib11 -(Sig21-Sig11*b21ib11)'*inv(ZZp)/b11];
VAR.et = (Bi*res')'; % Estimated monetary Shocks

PHI = mu1/b11;
GAM = inv(ED)*PHI;
E  = GAM*VAR.et(sum(m,2)~=0);
V  = m(sum(m,2)~=0)-E;
VAR.RM = inv(E'*E+V'*V)*E'*E;



