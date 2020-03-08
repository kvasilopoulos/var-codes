function VAR = doProxySVAR_single(VAR)
%Mertens-Ravn 2013 external instrument, modified to 
%1) handle different length of VAR and factors
%2) Multiple instruments to explain the same variable

 X      = lagmatrix(VAR.vars,1:VAR.p);
 X      = X(VAR.p+1:end,:);
 Y      = VAR.vars(VAR.p+1:end,:);
 [VAR.T,VAR.n] = size(Y);
 [VAR.T_m,VAR.n_m] = size(VAR.proxies);

 if VAR.switch_extern==1
      VAR.n_e= length(VAR.extern_vars);
      for gg=1:VAR.n_e
        Y_e{1,gg}    = [VAR.vars(VAR.smpl_min_VAR_e(1,gg)-VAR.smpl_min_VAR+VAR.p+1:end,:) VAR.extern{1,gg}(VAR.p+1:end,:)];
        X_e{1,gg}    = lagmatrix([VAR.vars(VAR.smpl_min_VAR_e(1,gg)-VAR.smpl_min_VAR+1:end,:) VAR.extern{1,gg}],1:VAR.p);
        X_e{1,gg}    = X_e{1,gg}(VAR.p+1:end,:);
        VAR.T_e(1,gg)   = size(X_e{1,gg},1);
        VAR.T_m_e(1,gg) = min(VAR.T_m,VAR.T_e(1,gg));
      end;
      VAR.n_e= length(VAR.extern_vars);
      VAR.n_ts  =   length(VAR.term_spreads);
      VAR.n_er  =   length(VAR.excess_return);
      VAR.n_rr  =   length(VAR.real_rates);
      if VAR.switch_exp==1
        VAR.n_exp =   length(VAR.exp_rates);
      end;
 end;
 
 % number of proxies
 VAR.k  = 1;
 %Assuming proxies start at least p periods later
 VAR.m  = VAR.proxies(1:end,:);
 
% Run VAR
%%%%%%%%%%%%
VAR.bet=[X ones(length(X),1)]\Y; 
VAR.res = Y-[X ones(length(X),1)]*VAR.bet;
        
if VAR.switch_extern==1
    for gg=1:VAR.n_e
        VAR.bet_e{1,gg}=[X_e{1,gg} ones(size(X_e{1,gg},1),1)]\Y_e{1,gg}; 
        VAR.res_e{1,gg}= Y_e{1,gg}-[X_e{1,gg} ones(size(X_e{1,gg},1),1)]*VAR.bet_e{1,gg};
    end;
end;
      
% Identification
%%%%%%%%%%%%%%%%%

% Only the restricted sample is used for identification
VAR.Sigma_m = (VAR.res(VAR.T-VAR.T_m-VAR.T_m_end+1:VAR.T-VAR.T_m_end,:)'*VAR.res(VAR.T-VAR.T_m-VAR.T_m_end+1:VAR.T-VAR.T_m_end,:))/(VAR.T_m-VAR.n*VAR.p-1);
if VAR.switch_extern==1
    for gg=1:VAR.n_e
        VAR.Sigma_m_e{1,gg} = (VAR.res_e{1,gg}(VAR.T_e(1,gg)-VAR.T_m_e(1,gg)-VAR.T_m_end+1:VAR.T_e(1,gg)-VAR.T_m_end,:)'*...
            VAR.res_e{1,gg}(VAR.T_e(1,gg)-VAR.T_m_e(1,gg)-VAR.T_m_end+1:VAR.T_e(1,gg)-VAR.T_m_end,:))/(VAR.T_m_e(1,gg)-(VAR.n+1)*VAR.p-1);           
    end;
end;

Phib = [ones(VAR.T_m,1) VAR.m]\VAR.res(VAR.T-VAR.T_m-VAR.T_m_end+1:VAR.T-VAR.T_m_end,:);
if VAR.switch_extern==1
    for gg=1:VAR.n_e
        Phib_e{1,gg} = [ones(VAR.T_m_e(1,gg),1) VAR.m(VAR.T_m-VAR.T_m_e(1,gg)+1:end,:)]\VAR.res_e{1,gg}(VAR.T_e(1,gg)-VAR.T_m_e(1,gg)-VAR.T_m_end+1:VAR.T_e(1,gg)-VAR.T_m_end,:);
    end;
end;
Res_m = VAR.res(VAR.T-VAR.T_m-VAR.T_m_end+1:VAR.T-VAR.T_m_end,:)-[ones(VAR.T_m,1) VAR.m]*Phib;
Res_const   =   VAR.res(VAR.T-VAR.T_m-VAR.T_m_end+1:VAR.T-VAR.T_m_end,1)-ones(VAR.T_m,1)*(ones(VAR.T_m,1)\VAR.res(VAR.T-VAR.T_m-VAR.T_m_end+1:VAR.T-VAR.T_m_end,1));
XX_m    =   [ones(VAR.T_m,1) VAR.m];
%Calculate robust standard errors
SS_m  =   zeros(VAR.n_m+1,VAR.n_m+1);
for ii=1:VAR.T_m
    SS_m  =   SS_m+1/VAR.T_m*XX_m(ii,:)'*XX_m(ii,:)*Res_m(ii,1)^2;
end;
Avarb_m     =   inv(1/VAR.T_m*XX_m'*XX_m)*SS_m*inv(1/VAR.T_m*XX_m'*XX_m);
RR_m     =   [zeros(VAR.n_m,1) eye(VAR.n_m)];
WW_m    =   VAR.T_m*(RR_m*Phib(:,1))'*inv(RR_m*Avarb_m*RR_m')*(RR_m*Phib(:,1));
VAR.F_m_rob     =   WW_m/VAR.n_m;

SST_m = Res_const'*Res_const;
SSE_m = Res_m(:,1)'*Res_m(:,1);
VAR.F_m = ((SST_m-SSE_m)/VAR.n_m)/(SSE_m/(length(VAR.res(VAR.T-VAR.T_m-VAR.T_m_end+1:VAR.T-VAR.T_m_end,1))-(VAR.n_m+1)));
VAR.R2_m = (1-SSE_m/(SST_m));
VAR.R2adj_m = VAR.R2_m-(VAR.n_m/((length(VAR.res(VAR.T-VAR.T_m-VAR.T_m_end+1:VAR.T-VAR.T_m_end,1))-(VAR.n_m+1))))*(1-VAR.R2_m); 

%TSLS, get a forecast of u1
uhat1           =   [ones(VAR.T_m,1) VAR.m]*Phib(:,1);
b21ib11_TSLS    =   [ones(VAR.T_m,1) uhat1]\VAR.res(VAR.T-VAR.T_m-VAR.T_m_end+1:VAR.T-VAR.T_m_end,2:end);  
b21ib11_TSLS    =   b21ib11_TSLS(2:end,:)';

b21ib11     =   b21ib11_TSLS;
if VAR.switch_extern==1
    for gg=1:VAR.n_e
           uhat1_temp = [ones(VAR.T_m_e(1,gg),1) VAR.m(VAR.T_m-VAR.T_m_e(1,gg)+1:end,:)]*Phib_e{1,gg}(:,1);
           b21ib11_TSLS_e{1,gg}  =   [ones(VAR.T_m_e(1,gg),1) uhat1_temp]\VAR.res_e{1,gg}(VAR.T_e(1,gg)-VAR.T_m_e(1,gg)-VAR.T_m_end+1:VAR.T_e(1,gg)-VAR.T_m_end,2:end);
           b21ib11_TSLS_e{1,gg}  =   b21ib11_TSLS_e{1,gg}(2:end,:)';
           b21ib11_e{1,gg} = b21ib11_TSLS_e{1,gg};
    end;
end;

% Identification of b11 and b12 from the covariance matrix of the VAR
Sig11   = VAR.Sigma_m(1:VAR.k,1:VAR.k);
Sig21   = VAR.Sigma_m(VAR.k+1:VAR.n,1:VAR.k);
Sig22   = VAR.Sigma_m(VAR.k+1:VAR.n,VAR.k+1:VAR.n);
ZZp     = b21ib11*Sig11*b21ib11'-(Sig21*b21ib11'+b21ib11*Sig21')+Sig22;
b12b12p = (Sig21- b21ib11*Sig11)'*(ZZp\(Sig21- b21ib11*Sig11));
b11b11p = Sig11-b12b12p;
b11 = sqrt(b11b11p);
VAR.b1 = [b11; b21ib11*b11];
VAR.Phib = Phib;
if VAR.switch_extern==1
    for gg=1:VAR.n_e
            Sig11_temp   = VAR.Sigma_m_e{1,gg}(1:VAR.k,1:VAR.k);
            Sig21_temp   = VAR.Sigma_m_e{1,gg}(VAR.k+1:end,1:VAR.k);
            Sig22_temp   = VAR.Sigma_m_e{1,gg}(VAR.k+1:end,VAR.k+1:end);
            ZZp_temp     = b21ib11_e{1,gg}*Sig11_temp*b21ib11_e{1,gg}'-(Sig21_temp*b21ib11_e{1,gg}'+b21ib11_e{1,gg}*Sig21_temp')+Sig22_temp;
            b12b12p_temp = (Sig21_temp- b21ib11_e{1,gg}*Sig11_temp)'*(ZZp_temp\(Sig21_temp- b21ib11_e{1,gg}*Sig11_temp));
            b11b11p_temp = Sig11_temp-b12b12p_temp;
            b11_e{1,gg} = sqrt(b11b11p_temp);
            VAR.b1_e{1,gg} = [b11_e{1,gg}; b21ib11_e{1,gg}*b11_e{1,gg}];
    end;
end;

% Impulse Responses
%%%%%%%%%%%%%%%%%%%%
% initial shock: eps(1,1)=1
irs(VAR.p+1,:) = VAR.b1(:,1);

if VAR.switch_extern==1
    for gg=1:VAR.n_e
        irs_e_all{1,gg}(VAR.p+1,:) = VAR.b1_e{1,gg}(:,1);
        irs_e(VAR.p+1,gg) = VAR.b1_e{1,gg}(VAR.n+1,1);
    end;
end;
for jj=2:VAR.irhor+max(max(VAR.term_spreads_matur),max(VAR.real_rates_init+VAR.real_rates_matur-1))
    lvars = (irs(VAR.p+jj-1:-1:jj,:))';
    irs(VAR.p+jj,:) = lvars(:)'*VAR.bet(1:VAR.p*VAR.n,:);     
end
if VAR.switch_extern==1
    for gg=1:VAR.n_e
        for jj=2:VAR.irhor+max(max(VAR.term_spreads_matur),max(VAR.real_rates_init+VAR.real_rates_matur-1))
            lvars_temp = irs_e_all{1,gg}(VAR.p+jj-1:-1:jj,:)';
            irs_e_all{1,gg}(VAR.p+jj,:) = lvars_temp(:)'*VAR.bet_e{1,gg}(1:VAR.p*(VAR.n+1),:);
            irs_e(VAR.p+jj,gg)=irs_e_all{1,gg}(VAR.p+jj,VAR.n+1);
        end;
    end;
end;
VAR.irs   = irs(VAR.p+1:VAR.p+VAR.irhor,:);
if VAR.switch_extern==1
    VAR.irs_e = irs_e(VAR.p+1:VAR.p+VAR.irhor,:);
    for gg=1:VAR.n_e
        VAR.irs_e_all{1,gg}=irs_e_all{1,gg}(VAR.p+1:VAR.p+VAR.irhor,:);
    end;
    
    VAR.irs_e_matur = zeros(VAR.irhor,VAR.n_ts);
    for jj=1:VAR.n_ts
        nn=VAR.term_spreads_matur(1,jj)/12;
        weight_ts_vec=weight_vec(0.0567,nn);  %calculates weighting for risk neutral yield estimation, average FF between 1979-2012: 5.67%
        for ii=1:VAR.irhor
            if (VAR.term_spreads(1,jj)-VAR.n)>0
                gg=VAR.term_spreads(1,jj)-VAR.n;
                switch VAR.monpol_FF
                    case 'no'    %if FF is not the the monpolshock, FF is the first external
                        VAR.irs_e_matur(ii,jj) = weight_ts_vec*(irs_e(VAR.p+ii+VAR.term_spreads_init(1,jj)-1:VAR.p+ii+VAR.term_spreads_matur(1,jj)-1,1));
                    case 'yes' %pick FF from the VAR
                        VAR.irs_e_matur(ii,jj) = weight_ts_vec*(irs_e_all{1,gg}(VAR.p+ii+VAR.term_spreads_init(1,jj)-1:VAR.p+ii+VAR.term_spreads_matur(1,jj)-1,1));
                end;
            else %term spread from the original VAR
                switch VAR.monpol_FF
                    case 'no'    %if FF is not the the monpolshock, FF is the first external
                        VAR.irs_e_matur(ii,jj) = weight_ts_vec*(irs_e(VAR.p+ii+VAR.term_spreads_init(1,jj)-1:VAR.p+ii+VAR.term_spreads_matur(1,jj)-1,1));
                    case 'yes' %pick FF from the VAR
                        VAR.irs_e_matur(ii,jj) = weight_ts_vec*(irs(VAR.p+ii+VAR.term_spreads_init(1,jj)-1:VAR.p+ii+VAR.term_spreads_matur(1,jj)-1,1));
                end;
            end;
        end
    end;
%Pick the ones that I want to calculate as a term spread
    irs_all     =   [VAR.irs VAR.irs_e];
    if ~isempty(VAR.term_spreads)
        VAR.irs_ts=irs_all(:,VAR.term_spreads)-VAR.irs_e_matur;
    end;
%Pick the ones to calculate excess premia
    if ~isempty(VAR.term_spreads)
        VAR.irs_er=irs_all(:,VAR.excess_return)+VAR.irs_ts(:,VAR.excess_return_matur);
    end;
    VAR.irs_infl = zeros(VAR.irhor,VAR.n_rr);
%Calculate average expected inflation rates
    for jj=1:VAR.n_rr
        nn=VAR.real_rates_matur(1,jj)/12;
        weight_infl_vec=weight_vec(0.034968,nn);  %calculates weighting for risk neutral yield estimation, average infl between 1979-2012: 3.4968%
        for ii=1:VAR.irhor
            if (VAR.real_rates(1,jj)-VAR.n)>0
                gg=VAR.real_rates(1,jj)-VAR.n;
                VAR.irs_bkeven(ii,jj) = weight_infl_vec*(irs_e_all{1,gg}(VAR.p+ii+VAR.real_rates_init(1,jj):VAR.p+ii+VAR.real_rates_init(1,jj)-1+VAR.real_rates_matur(1,jj),2)-irs_e_all{1,gg}(VAR.p+ii+VAR.real_rates_init(1,jj)-1:VAR.p+ii+VAR.real_rates_init(1,jj)-1+VAR.real_rates_matur(1,jj)-1,2))*12;  %Forward looking Annualized CPI is the second
            else %pick infl from the original VAR
                VAR.irs_bkeven(ii,jj) = weight_infl_vec*(irs(VAR.p+ii+VAR.real_rates_init(1,jj):VAR.p+ii+VAR.real_rates_init(1,jj)-1+VAR.real_rates_matur(1,jj),2)-irs(VAR.p+ii+VAR.real_rates_init(1,jj)-1:VAR.p+ii+VAR.real_rates_init(1,jj)-1+VAR.real_rates_matur(1,jj)-1,2))*12;  %Forward looking Annualized CPI is the second
            end;
        end;
    end
%Pick the ones to calculate real rates
    if ~isempty(VAR.real_rates)
        VAR.irs_rr=irs_all(:,VAR.real_rates)-VAR.irs_bkeven;
    end;
    if VAR.switch_exp==1
%Create expectations
    irs_all_ext     =   [irs irs_e];
    VAR.irs_exp = zeros(VAR.irhor,VAR.n_exp);
    for ii=1:VAR.irhor
        for jj=1:VAR.n_exp
             VAR.irs_exp(ii,jj) = mean(irs_all_ext(VAR.p+ii:12:VAR.p+ii+VAR.exp_rates_matur(1,jj)-12,VAR.exp_rates(1,jj)));
        end;
    end
    end;
end;