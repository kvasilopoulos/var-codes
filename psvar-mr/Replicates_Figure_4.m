clear all; close all; addpath('auxfiles');
 
load DATASET;

nboot  = 1000;  % Number of Bootstrap Samples (equals 10000 in the paper)
clevel = 95;    % Bootstrap Percentile Shown

% VAR specification
%%%%%%%%%%%%%%%%%%%%
 VAR.p      = 4;                                 % Number of Lags
 VAR.irhor  = 20;                                % Impulse Response Horizon
 VAR.select_vars      = {'APITR','ACITR','GOV','RGDP','DEBT','FF','PLEVEL','NBR'};
 VAR.vars             = DATASET.TSERIES(:,cell2mat(values(DATASET.MAP,VAR.select_vars)));
 VAR.MAP              = containers.Map([VAR.select_vars,{'INFL'}],[1:size(VAR.vars,2) 9]);
 VAR.proxies          = DATASET.TSERIES(:,cell2mat(values(DATASET.MAP,{'m_PI','m_CI'})));

 % Estimation
 %%%%%%%%%%%%%
 % APITR Ordered first
 VARO1      = VAR;
 VARO1.ord  = [1 2] ;
 VARO1      = doProxySVAR(VARO1,DATASET);
 VARbsO1    = doProxySVARbootstrap(VARO1,nboot,clevel,DATASET);
 
 % ACITR Ordered first
 VARO2      = VAR;
 VARO2.ord  = [2 1] ;
 VARO2      = doProxySVAR(VARO2,DATASET);
 VARbsO2    = doProxySVARbootstrap(VARO2,nboot,clevel,DATASET);

 % Figure 4: APITR Cut
 %%%%%%%%%%%%%%%%%%%%%
 FIG.axes(:,cell2mat(values(VARO1.MAP,{'RGDP'})))   = [-1;3];
 FIG.axes(:,cell2mat(values(VARO1.MAP,{'DEBT'})))   = [-4;3];
 FIG.axes(:,cell2mat(values(VARO1.MAP,{'FF'})))     = [-1.5;1.5];
 FIG.axes(:,cell2mat(values(VARO1.MAP,{'INFL'})))   = [-2;2]; 
 plotdisplay = {'RGDP','DEBT','INFL','FF'};
 doFigureA(VARO1,VARbsO1,VARO2,VARbsO2,FIG,plotdisplay,DATASET,1);
 
 % Figure 4: ACITR Cut
 %%%%%%%%%%%%%%%%%%%%%
 FIG.axes(:,cell2mat(values(VARO1.MAP,{'RGDP'})))   = [-0.5;1.5];
 FIG.axes(:,cell2mat(values(VARO1.MAP,{'DEBT'})))   = [-1.5;1];
 FIG.axes(:,cell2mat(values(VARO1.MAP,{'FF'})))     = [-1;1];
 FIG.axes(:,cell2mat(values(VARO1.MAP,{'INFL'})))   = [-0.8;0.4];
 plotdisplay = {'RGDP','DEBT','INFL','FF'};
 doFigureA(VARO1,VARbsO1,VARO2,VARbsO2,FIG,plotdisplay,DATASET,2);
