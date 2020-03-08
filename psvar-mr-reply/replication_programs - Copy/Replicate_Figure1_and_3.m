clear all; close all; 
addpath('auxiliary files')

% Create DATASET
%%%%%%%%%%%%%%%%
DATASET.TSERIES = xlsread('datasets/MR2013_AER_DATASET.xlsx','Quarterly');
DATASET.LABEL   = {'DATES','T_PI','T_CI','m_PI','m_CI','APITR','ACITR','PITB','CITB','GOV'};
DATASET.VALUE   = [  1,       2,     3,     4,    5 ,    6 ,     7 ,     8  ,   9 ,    10 ];
DATASET.UNIT    = [  0,       2,     2,     2,    2 ,    2,      2 ,     1  ,   1 ,    1  ];

DATASET.LABEL   = [DATASET.LABEL,{'RGDP','DEBT','FF','PLEVEL','NBR','EMP','HperW','LF','CNDSV','CD','INR','IR'}];
DATASET.VALUE   = [DATASET.VALUE,   11,    12  ,  13 ,  14 ,   15  , 16  ,   17 ,  18 ,  19  ,  20 ,  21,   22 ];
DATASET.UNIT    = [DATASET.UNIT,     1  ,   1  ,   2 ,  1  ,   1  ,  1  ,   1   ,  1  ,  1  ,   1 ,   1 ,  1  ]; 
    
DATASET.LABEL   = [DATASET.LABEL,{'UNR','INFL','PITREV','CITREV'}];
DATASET.VALUE   = [DATASET.VALUE,   23,    24 ,  25   ,  26 ];
DATASET.UNIT    = [DATASET.UNIT,     2  ,   2 ,  1    ,   1 ]; 

DATASET.MAP = containers.Map(DATASET.LABEL,DATASET.VALUE);

% VAR specification
%%%%%%%%%%%%%%%%%%%%
VAR.p                = 4;  % Number of Lags
VAR.irhor            = 20; % Impulse Response Horizon
VAR.select_vars      = {'APITR','ACITR','PITB','CITB','GOV','RGDP','DEBT'};
VAR.vars             = DATASET.TSERIES(:,cell2mat(values(DATASET.MAP,VAR.select_vars)));
VAR.MAP              = containers.Map(VAR.select_vars,1:size(VAR.vars,2));
VAR.proxies          = DATASET.TSERIES(:,cell2mat(values(DATASET.MAP,{'m_PI','m_CI'})));
VAR.DET              = ones(length(VAR.vars),1); % Deterministic Terms

VAR        = doProxySVAR(VAR);

% Inference:
%%%%%%%%%%%%
% method  1: Mertens and Ravn (2013) wild bootstrap 
%         2: Montiel-Olea Stock Watson (2016) parametric bootstrap
%         3: Delta Method 
%         4: Montiel-Olea Stock Watson (2016) asy weak IV
%         5: Jentsch and Lunsford Moving Block Bootstrap
%         6: Jentsch and Lunsford Moving Block Bootstrap (adjusted to allow non zero-mean proxies)

nboot     = 5000;        % Number of Bootstrap Samples (Paper does 5000)
clevel    = 68;          % Bootstrap Percentile Shown
BlockSize = 19;          % size of blocks in the MBB bootstrap
seed      = 2;           % seed for random number generator
rng(seed);               % iniate the random number generator

VARci_wildbs   = doProxySVARci(VAR,clevel,1,nboot); 
VARci_mswbs    = doProxySVARci(VAR,clevel,2,nboot);
VARci_delta    = doProxySVARci(VAR,clevel,3);
VARci_mbb      = doProxySVARci(VAR,clevel,5,nboot,BlockSize);

% Make Plots
legendflag = 1;

plotdisplay = {'RGDP'};
FIGLABELS = {'Output'};

for shock = [1 4]
    if shock ==1
        FIG.AXIS = [-1 4];
        name = 'Figure1_left'
    elseif shock ==4 
        FIG.AXIS = [-0.5 2.5];
        name = 'Figure1_right'
    end       
    VARci.irsL  = VARci_delta.irsL(:,:,shock);
    VARci.irsH  = VARci_delta.irsH(:,:,shock);
    VARci.irsL2 = VARci_wildbs.irsL(:,:,shock);
    VARci.irsH2 = VARci_wildbs.irsH(:,:,shock);
    VARci.irsL3 = VARci_mswbs.irsL(:,:,shock);
    VARci.irsH3 = VARci_mswbs.irsH(:,:,shock);
    VARci.irsL4 = VARci_mbb.irsL(:,:,shock);
    VARci.irsH4 = VARci_mbb.irsH(:,:,shock);
    legendlabels = {'Delta Method','Wild Bootstrap','Par. Bootstrap','Block Bootstrap'};
    do_IRS_Figure1
    str=strcat('figures/',name,'_',plotdisplay{nvar});
             saveas(gcf,str,'epsc');
end

for shock = [1 4]
    if shock ==1
        FIG.AXIS = [-1 4];
        name = 'Figure3_left'
    elseif shock ==4 
        FIG.AXIS = [-0.5 2.5];
        name = 'Figure3_right'
    end       
    VARci.irsL  = VARci_delta.irsL(:,:,shock);
    VARci.irsH  = VARci_delta.irsH(:,:,shock);
    VARci.irsL2 = VARci_wildbs.irsLhall(:,:,shock);
    VARci.irsH2 = VARci_wildbs.irsHhall(:,:,shock);
    VARci.irsL3 = VARci_mswbs.irsLhall(:,:,shock);
    VARci.irsH3 = VARci_mswbs.irsHhall(:,:,shock);
    VARci.irsL4 = VARci_mbb.irsLhall(:,:,shock);
    VARci.irsH4 = VARci_mbb.irsHhall(:,:,shock);
    legendlabels = {'Delta Method','Wild Bootstrap','Par. Bootstrap','Block Bootstrap'};
    do_IRS_Figure1
    str=strcat('figures/',name,'_',plotdisplay{nvar});
             saveas(gcf,str,'epsc');
end



