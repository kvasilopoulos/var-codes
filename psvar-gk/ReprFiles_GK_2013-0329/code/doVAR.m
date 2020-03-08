function [VAR,VARChol,VARbs,VARCholbs]=doVAR(VAR,DATASET_VAR,DATASET_FACTORS,nboot,clevel)

% VAR specification
%%%%%%%%%%%%%%%%%%%%
 VAR.vars             = DATASET_VAR.TSERIES(VAR.smpl_min_VAR:VAR.smpl_max_VAR,cell2mat(values(DATASET_VAR.MAP,VAR.select_vars)));
 VAR.proxies          = DATASET_FACTORS.TSERIES(VAR.smpl_min_FACTORS:VAR.smpl_max_FACTORS,cell2mat(values(DATASET_FACTORS.MAP,VAR.select_factors)));
 VAR.T                = length(VAR.vars)-VAR.p;
 VAR.n                = size(VAR.vars,2);
 [VAR.T_m,VAR.T_n]    = size(VAR.proxies);
 VAR.T_m_end          = VAR.smpl_max_VAR-VAR.smpl_max_VAR_factors;      %The difference between where the VAR and the factors end in the VAR sample
 VAR.year             = DATASET_VAR.TSERIES(VAR.smpl_min_VAR:VAR.smpl_max_VAR,cell2mat(values(DATASET_VAR.MAP,{'YEAR'})));
 VAR.month            = DATASET_VAR.TSERIES(VAR.smpl_min_VAR:VAR.smpl_max_VAR,cell2mat(values(DATASET_VAR.MAP,{'MONTH'})));
 
 if VAR.switch_extern==1
    VAR.n_e              = length(VAR.extern_vars);
    for gg=1:VAR.n_e 
       VAR.extern{:,gg}           = DATASET_VAR.TSERIES(VAR.smpl_min_VAR_e(1,gg):VAR.smpl_max_VAR,cell2mat(values(DATASET_VAR.MAP,VAR.extern_vars(1,gg))));
    end;
    VAR.smpl_min_FACTORS_e  =   max(VAR.smpl_min_VAR_e,VAR.smpl_min_FACTORS*ones(1,VAR.n_e));
 end;

 VAR=doProxySVAR_single(VAR);
 VARbs=doProxySVARbootstrap_single(VAR,nboot,clevel);
 VARChol=VAR;
 VARChol.vars   =   VAR.vars(:,VAR.chol_order);
 VARChol=doCholSVAR_single(VARChol);
 VARCholbs=doCholSVARbootstrap_single(VARChol,nboot,clevel);
