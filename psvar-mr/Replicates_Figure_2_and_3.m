clear all; close all; addpath('auxfiles');
 
load DATASET;

nboot  = 1000;  % Number of Bootstrap Samples (equals 10000 in the paper)
clevel = 95;    % Bootstrap Percentile Shown

% VAR specification
%%%%%%%%%%%%%%%%%%%%
 VAR.p      = 4;                                 % Number of Lags
 VAR.irhor  = 20;                                % Impulse Response Horizon
 VAR.select_vars      = {'APITR','ACITR','PITB','CITB','GOV','RGDP','DEBT'};
 VAR.vars             = DATASET.TSERIES(:,cell2mat(values(DATASET.MAP,VAR.select_vars)));
 VAR.MAP              = containers.Map([VAR.select_vars,{'PITREV','CITREV'}],[1:size(VAR.vars,2) 8 9]);
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

 % Figure 2
 %%%%%%%%%%%
 FIG.axes(:,cell2mat(values(VARO1.MAP,{'RGDP'})))   = [-1;3.5];
 FIG.axes(:,cell2mat(values(VARO1.MAP,{'GOV'})))    = [-5;4];
 FIG.axes(:,cell2mat(values(VARO1.MAP,{'APITR'})))  = [-1.2;0.5];
 FIG.axes(:,cell2mat(values(VARO1.MAP,{'ACITR'})))  = [-1.5;3];
 FIG.axes(:,cell2mat(values(VARO1.MAP,{'PITB'})))   = [-1;3];
 FIG.axes(:,cell2mat(values(VARO1.MAP,{'PITREV'}))) = [-7;4];
 plotdisplay = {'APITR','RGDP','PITB','PITREV','ACITR','GOV'};
 doFigureA(VARO1,VARbsO1,VARO2,VARbsO2,FIG,plotdisplay,DATASET,1);
 
 % Figure 3
 %%%%%%%%%%%
 FIG.axes(:,cell2mat(values(VARO1.MAP,{'RGDP'})))   = [-0.5;2];
 FIG.axes(:,cell2mat(values(VARO1.MAP,{'GOV'})))    = [-4;4];
 FIG.axes(:,cell2mat(values(VARO1.MAP,{'APITR'})))  = [-0.4;0.6];
 FIG.axes(:,cell2mat(values(VARO1.MAP,{'ACITR'})))  = [-1.2;1];
 FIG.axes(:,cell2mat(values(VARO1.MAP,{'CITB'})))   = [-3;8];
 FIG.axes(:,cell2mat(values(VARO1.MAP,{'CITREV'}))) = [-5;8];
 plotdisplay = {'ACITR','RGDP','CITB','CITREV','APITR','GOV'};
 doFigureA(VARO1,VARbsO1,VARO2,VARbsO2,FIG,plotdisplay,DATASET,2);
 
 
