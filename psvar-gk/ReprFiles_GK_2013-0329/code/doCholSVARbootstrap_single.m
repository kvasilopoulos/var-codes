function VARbs = doCholSVARbootstrap_single(VAR,nboot,clevel)

% Wild Bootstrap
  jj=1;
    % removes a constant
    res = detrend(VAR.res,'constant');
     while jj<nboot+1
       %T x 1 -1 or 1 random numbers with prob 0.5-0.5
       rr = 1-2*(rand(VAR.T,1)>0.5);
       %T x n randomly choose the sign of the time T shocks (all)
       resb = (res.*(rr*ones(1,VAR.n)))';
       
       varsb = zeros(VAR.p+VAR.T,VAR.n);
       % Initial values
       varsb(1:VAR.p,:)=VAR.vars(1:VAR.p,:);
        for j=VAR.p+1:VAR.p+VAR.T
            lvars = (varsb(j-1:-1:j-VAR.p,:))';
            varsb(j,:) = lvars(:)'*VAR.bet(1:VAR.p*VAR.n,:)+VAR.bet(VAR.p*VAR.n+1,:)+resb(:,j-VAR.p)';
        end
        VARBS = VAR;
        VARBS.vars = varsb;        
        VARBS = doCholSVAR_single(VARBS);
        
        for j=1:VAR.k;
        irs = VARBS.irs(:,:,j);
        VARbs.irs(:,jj,j) = irs(:);
        end
                 
      jj=jj+1;   

     end  
    
 for j=1:VAR.k  
 VARbs.irsL(:,:,j)=reshape(quantile(VARbs.irs(:,:,j)',(1-clevel/100)/2),VAR.irhor,size(VARBS.irs,2));
 VARbs.irsH(:,:,j)=reshape(quantile(VARbs.irs(:,:,j)',1-(1-clevel/100)/2),VAR.irhor,size(VARBS.irs,2));
 end


