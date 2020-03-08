function VARbs = doProxySVARbootstrap_single(VAR,nboot,clevel)

% Wild Bootstrap
  jj=1;
    % removes a constant
    res = detrend(VAR.res,'constant');
    if VAR.switch_extern==1
        for gg=1:VAR.n_e
            res_e{1,gg} = detrend(VAR.res_e{1,gg},'constant');
        end;
    end;
     while jj<nboot+1
       %T x 1 -1 or 1 random numbers with prob 0.5-0.5
       rr = 1-2*(rand(VAR.T,1)>0.5);
       %T x n randomly choose the sign of the time T shocks (all)
       resb = (res.*(rr*ones(1,VAR.n)))';
       varsb = zeros(VAR.p+VAR.T,VAR.n);

       if VAR.switch_extern==1
           for gg=1:VAR.n_e
                resb_e{1,gg} = (res_e{1,gg}.*(rr(VAR.T-VAR.T_e(1,gg)+1:end,:)*ones(1,VAR.n+1)))';
                varsb_e{1,gg} = zeros(VAR.p+VAR.T_e(1,gg),1);               
           end;
       end;
       
       % Initial values
       varsb(1:VAR.p,:)=VAR.vars(1:VAR.p,:);
        for j=VAR.p+1:VAR.p+VAR.T
            lvars = (varsb(j-1:-1:j-VAR.p,:))';
            varsb(j,:) = lvars(:)'*VAR.bet(1:VAR.p*VAR.n,:)+VAR.bet(VAR.p*VAR.n+1:end,:)+resb(:,j-VAR.p)';
        end
        
        if VAR.switch_extern==1
        for gg=1:VAR.n_e
            varsb_e{1,gg}(1:VAR.p,:)=VAR.extern{1,gg}(1:VAR.p,:);
            clear varsb_temp;
            varsb_temp(1:VAR.p,:)=[VAR.vars(1+VAR.T-VAR.T_e(1,gg):VAR.p+VAR.T-VAR.T_e(1,gg),:) VAR.extern{1,gg}(1:VAR.p,:)];   %Create VARs with beta_e parameters; note: not the same as beta VAR
            for j=VAR.p+1:VAR.p+VAR.T_e(1,gg)
                ii=j+VAR.T-VAR.T_e(1,gg);
%Separate VAR for all auxiliary regressions                            
                lvars = (varsb_temp(j-1:-1:j-VAR.p,:))';
                varsb_temp(j,:) = lvars(:)'*VAR.bet_e{1,gg}(1:VAR.p*(VAR.n+1),:)+VAR.bet_e{1,gg}(VAR.p*(VAR.n+1)+1:end,:)+resb_e{1,gg}(:,j-VAR.p)';
                varsb_e{1,gg}(j,:) = varsb_temp(j,VAR.n+1); 
            end;
        end;
        end;
              
        VARBS = VAR;
        VARBS.vars = varsb;
        if VAR.switch_extern==1
            VARBS.extern = varsb_e;
        end;
        VARBS.proxies = [VAR.m.*(rr(VAR.T-VAR.T_m-VAR.T_m_end+1:VAR.T-VAR.T_m_end,1)*ones(1,size(VAR.proxies,2)))];
        VARBS = doProxySVAR_single(VARBS);
        
        for j=1:VAR.k;
        irs = VARBS.irs(:,:,j);
        VARbs.irs(:,jj,j) = irs(:);

        string_cell = {''};
        no_string_cell=length(string_cell);
        
        if VAR.switch_extern==1
            for ii=1:no_string_cell
                eval(['irs' string_cell{1,ii} '_e = VARBS.irs' string_cell{1,ii} '_e(:,:,j);']);
                eval(['irs' string_cell{1,ii} '_e_matur = VARBS.irs' string_cell{1,ii} '_e_matur(:,:,j);']);        
                eval(['irs' string_cell{1,ii} '_ts = VARBS.irs' string_cell{1,ii} '_ts(:,:,j);']);        
                eval(['irs' string_cell{1,ii} '_er = VARBS.irs' string_cell{1,ii} '_er(:,:,j);']);        
                eval(['irs' string_cell{1,ii} '_rr = VARBS.irs' string_cell{1,ii} '_rr(:,:,j);']);        
                eval(['irs' string_cell{1,ii} '_bkeven = VARBS.irs' string_cell{1,ii} '_bkeven(:,:,j);']);        
                if VAR.switch_exp==1
                    eval(['irs' string_cell{1,ii} '_exp = VARBS.irs' string_cell{1,ii} '_exp(:,:,j);']);        
                end;
                for gg=1:VAR.n_e
                    eval(['irs' string_cell{1,ii} '_e_all{1,' num2str(gg) '} = VARBS.irs' string_cell{1,ii} '_e_all{1,' num2str(gg) '}(:,:,j);']);
                end;
                eval(['VARbs.irs' string_cell{1,ii} '_e(:,jj,j) = irs' string_cell{1,ii} '_e(:);']);
                eval(['VARbs.irs' string_cell{1,ii} '_e_matur(:,jj,j) = irs' string_cell{1,ii} '_e_matur(:);']);
                eval(['VARbs.irs' string_cell{1,ii} '_ts(:,jj,j) = irs' string_cell{1,ii} '_ts(:);']);
                eval(['VARbs.irs' string_cell{1,ii} '_er(:,jj,j) = irs' string_cell{1,ii} '_er(:);']);
                eval(['VARbs.irs' string_cell{1,ii} '_rr(:,jj,j) = irs' string_cell{1,ii} '_rr(:);']);
                eval(['VARbs.irs' string_cell{1,ii} '_bkeven(:,jj,j) = irs' string_cell{1,ii} '_bkeven(:);']);
                if VAR.switch_exp==1
                    eval(['VARbs.irs' string_cell{1,ii} '_exp(:,jj,j) = irs' string_cell{1,ii} '_exp(:);']);
                end;
                for gg=1:VAR.n_e
                    eval(['VARbs.irs' string_cell{1,ii} '_e_all{1,' num2str(gg) '}(:,jj,j) = irs' string_cell{1,ii} '_e_all{1,' num2str(gg) '}(:);']);
                end;                
            end;
        end;
        end
        %bsRM(:,jj) = VARBS.RM ;
                 
      jj=jj+1;   

     end  
    
 for j=1:VAR.k  
 VARbs.irsL(:,:,j)=reshape(quantile(VARbs.irs(:,:,j)',(1-clevel/100)/2),VAR.irhor,size(VARBS.irs,2));
 VARbs.irsH(:,:,j)=reshape(quantile(VARbs.irs(:,:,j)',1-(1-clevel/100)/2),VAR.irhor,size(VARBS.irs,2));
    if VAR.switch_extern==1
        for ii=1:no_string_cell
            eval(['VARbs.irs' string_cell{1,ii} 'L_e(:,:,j)=reshape(quantile(VARbs.irs' string_cell{1,ii} '_e(:,:,j)'',(1-clevel/100)/2),VAR.irhor,size(VARBS.irs' string_cell{1,ii} '_e,2));']);
            eval(['VARbs.irs' string_cell{1,ii} 'H_e(:,:,j)=reshape(quantile(VARbs.irs' string_cell{1,ii} '_e(:,:,j)'',1-(1-clevel/100)/2),VAR.irhor,size(VARBS.irs' string_cell{1,ii} '_e,2));']);
            eval(['VARbs.irs' string_cell{1,ii} 'L_e_matur(:,:,j)=reshape(quantile(VARbs.irs' string_cell{1,ii} '_e_matur(:,:,j)'',(1-clevel/100)/2),VAR.irhor,size(VARBS.irs' string_cell{1,ii} '_e_matur,2));']);
            eval(['VARbs.irs' string_cell{1,ii} 'H_e_matur(:,:,j)=reshape(quantile(VARbs.irs' string_cell{1,ii} '_e_matur(:,:,j)'',1-(1-clevel/100)/2),VAR.irhor,size(VARBS.irs' string_cell{1,ii} '_e_matur,2));']);
            eval(['VARbs.irs' string_cell{1,ii} 'L_ts(:,:,j)=reshape(quantile(VARbs.irs' string_cell{1,ii} '_ts(:,:,j)'',(1-clevel/100)/2),VAR.irhor,size(VARBS.irs' string_cell{1,ii} '_ts,2));']);
            eval(['VARbs.irs' string_cell{1,ii} 'H_ts(:,:,j)=reshape(quantile(VARbs.irs' string_cell{1,ii} '_ts(:,:,j)'',1-(1-clevel/100)/2),VAR.irhor,size(VARBS.irs' string_cell{1,ii} '_ts,2));']);
            eval(['VARbs.irs' string_cell{1,ii} 'L_er(:,:,j)=reshape(quantile(VARbs.irs' string_cell{1,ii} '_er(:,:,j)'',(1-clevel/100)/2),VAR.irhor,size(VARBS.irs' string_cell{1,ii} '_er,2));']);
            eval(['VARbs.irs' string_cell{1,ii} 'H_er(:,:,j)=reshape(quantile(VARbs.irs' string_cell{1,ii} '_er(:,:,j)'',1-(1-clevel/100)/2),VAR.irhor,size(VARBS.irs' string_cell{1,ii} '_er,2));']);
            eval(['VARbs.irs' string_cell{1,ii} 'L_rr(:,:,j)=reshape(quantile(VARbs.irs' string_cell{1,ii} '_rr(:,:,j)'',(1-clevel/100)/2),VAR.irhor,size(VARBS.irs' string_cell{1,ii} '_rr,2));']);
            eval(['VARbs.irs' string_cell{1,ii} 'H_rr(:,:,j)=reshape(quantile(VARbs.irs' string_cell{1,ii} '_rr(:,:,j)'',1-(1-clevel/100)/2),VAR.irhor,size(VARBS.irs' string_cell{1,ii} '_rr,2));']);
            eval(['VARbs.irs' string_cell{1,ii} 'L_bkeven(:,:,j)=reshape(quantile(VARbs.irs' string_cell{1,ii} '_bkeven(:,:,j)'',(1-clevel/100)/2),VAR.irhor,size(VARBS.irs' string_cell{1,ii} '_bkeven,2));']);
            eval(['VARbs.irs' string_cell{1,ii} 'H_bkeven(:,:,j)=reshape(quantile(VARbs.irs' string_cell{1,ii} '_bkeven(:,:,j)'',1-(1-clevel/100)/2),VAR.irhor,size(VARBS.irs' string_cell{1,ii} '_bkeven,2));']);
            if VAR.switch_exp==1
                eval(['VARbs.irs' string_cell{1,ii} 'L_exp(:,:,j)=reshape(quantile(VARbs.irs' string_cell{1,ii} '_exp(:,:,j)'',(1-clevel/100)/2),VAR.irhor,size(VARBS.irs' string_cell{1,ii} '_exp,2));']);
                eval(['VARbs.irs' string_cell{1,ii} 'H_exp(:,:,j)=reshape(quantile(VARbs.irs' string_cell{1,ii} '_exp(:,:,j)'',1-(1-clevel/100)/2),VAR.irhor,size(VARBS.irs' string_cell{1,ii} '_exp,2));']);
            end;
            for gg=1:VAR.n_e
                    eval(['VARbs.irs' string_cell{1,ii} 'L_e_all{1,' num2str(gg) '}(:,:,j)=reshape(quantile(VARbs.irs' string_cell{1,ii} '_e_all{1,' num2str(gg) '}(:,:,j)'',(1-clevel/100)/2),VAR.irhor,size(VARBS.irs' string_cell{1,ii} '_e_all{1,' num2str(gg) '},2));']);
                    eval(['VARbs.irs' string_cell{1,ii} 'H_e_all{1,' num2str(gg) '}(:,:,j)=reshape(quantile(VARbs.irs' string_cell{1,ii} '_e_all{1,' num2str(gg) '}(:,:,j)'',1-(1-clevel/100)/2),VAR.irhor,size(VARBS.irs' string_cell{1,ii} '_e_all{1,' num2str(gg) '},2));']);
            end;                
        end;
    end;
 end