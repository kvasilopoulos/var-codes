function VARbs = doPVAR(VAR,DAT,nboot,clevel)

% Wild Bootstrap
  jj=1;
     while jj<nboot+1
       rr = 1-2*(rand(VAR.T,1)>0.5);
       eb     =  VAR.e.*(rr*ones(1,VAR.n-VAR.k));
       b11eTb = VAR.b11eT.*(rr*ones(1,VAR.k));
       resb =  [eye(VAR.k); VAR.b21ib11]*b11eTb'+VAR.D(:,VAR.k+1:VAR.n)*eb';
       
       varsb = zeros(VAR.p+VAR.T,VAR.n);
       varsb(1:VAR.p,:)=VAR.vars(1:VAR.p,:);
        for j=VAR.p+1:VAR.p+VAR.T
        lvars = (varsb(j-1:-1:j-VAR.p,:))';
        varsb(j,:) = lvars(:)'*VAR.bet(1:VAR.p*VAR.n,:)+VAR.DET(j,:)*VAR.bet(VAR.p*VAR.n+1:end,:)+resb(:,j-VAR.p)';     
        end
              
        VARBS = VAR;
        VARBS.vars = varsb;
        VARBS.mshocks = [VAR.MSHOCKS(1:VAR.p,:); VAR.MSHOCKS.*(rr*ones(1,VAR.k))];
        VARBS = doPVAR(VARBS);
        
        for j=1:VAR.k;
        irs = VARBS.irs(:,:,j);
        VARbs.irs(:,jj,j) = irs(:);
        end
        
        VARbs.irsg(:,jj)   = VARBS.irsg(:);
        VARbs.sigmaG(:,jj) = VARBS.sigmaG;
        VARbs.sigmaY(:,jj) = VARBS.sigmaY;
        VARbs.SigmaT(:,jj) = VARBS.SigmaT(:);
        VARbs.thetaY(:,jj) = VARBS.thetaY;
        VARbs.thetaG(:,jj) = VARBS.thetaG;
        VARbs.gammaT(:,jj) = VARBS.gammaT;
        VARbs.zetaT(:,jj)  = VARBS.zetaT;
        VARbs.zetaG(:,jj)  = VARBS.zetaG;
        VARbs.RM(:,jj)     = VARBS.RM(:);
         
        if VAR.k==1
        VARbs.irsTRY(:,jj) = irs(:,1)-irs(:,3);
        VARbs.BS_mshocks(:,jj) = VARBS.mshocks;
        truemshocks = [VAR.MSHOCKS(1:VAR.p,:); VAR.SigmaT*(VAR.D(1,1)*100*VAR.et.*rr/DAT.TRY) ];

        truemshocks(VARBS.mshocks==0)=0;
        VARbs.BS_mshocks_true(:,jj) = truemshocks;          
        VARbs.BS_VARS(:,jj) = varsb(:);
        end
      jj=jj+1;   
     end  
    
 for j=1:VAR.k    
 VARbs.irsH(:,:,j)=reshape(quantile(VARbs.irs(:,:,j)',(1-clevel/100)/2),VAR.irhor,VAR.n);
 VARbs.irsL(:,:,j)=reshape(quantile(VARbs.irs(:,:,j)',1-(1-clevel/100)/2),VAR.irhor,VAR.n);
 end
 if VAR.k==1
  VARbs.irsTRYH=reshape(quantile(VARbs.irsTRY',(1-clevel/100)/2),VAR.irhor,1);
  VARbs.irsTRYL=reshape(quantile(VARbs.irsTRY',1-(1-clevel/100)/2),VAR.irhor,1);
    
 end

  VARbs.irsgH=reshape(quantile(VARbs.irsg',(1-clevel/100)/2),VAR.irhor,VAR.n);
  VARbs.irsgL=reshape(quantile(VARbs.irsg',1-(1-clevel/100)/2),VAR.irhor,VAR.n);
   
   VARbs.sigmaGci = [quantile(VARbs.sigmaG',(1-clevel/100)/2) quantile(VARbs.sigmaG',1-(1-clevel/100)/2)];
   VARbs.sigmaYci = [quantile(VARbs.sigmaY',(1-clevel/100)/2) quantile(VARbs.sigmaY',1-(1-clevel/100)/2)];
   VARbs.SigmaTci = [quantile(VARbs.SigmaT',(1-clevel/100)/2) quantile(VARbs.SigmaT',1-(1-clevel/100)/2)];
   VARbs.thetaYci = [quantile(VARbs.thetaY',(1-clevel/100)/2)' quantile(VARbs.thetaY',1-(1-clevel/100)/2)'];
   VARbs.thetaGci  = [quantile(VARbs.thetaG',(1-clevel/100)/2)' quantile(VARbs.thetaG',1-(1-clevel/100)/2)'];
   VARbs.gammaTci  = [quantile(VARbs.gammaT',(1-clevel/100)/2)' quantile(VARbs.gammaT',1-(1-clevel/100)/2)'];
   VARbs.zetaGci   = [quantile(VARbs.zetaG',(1-clevel/100)/2) quantile(VARbs.zetaG',1-(1-clevel/100)/2)];
   VARbs.zetaTci   = [quantile(VARbs.zetaT',(1-clevel/100)/2)' quantile(VARbs.zetaT',1-(1-clevel/100)/2)'];
   VARbs.RMci      = [quantile(VARbs.RM',(1-clevel/100)/2)' quantile(VARbs.RM',1-(1-clevel/100)/2)'];





