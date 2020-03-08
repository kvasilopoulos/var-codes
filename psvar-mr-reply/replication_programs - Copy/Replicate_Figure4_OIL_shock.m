clear all; close all; 
addpath('auxiliary files')

% Create DATASET
%%%%%%%%%%%%%%%%

mdata = xlsread('datasets/Kilian_data.xlsx');
mdata = mdata(2:end,:);

DATASET.TSERIES   = [mdata(:,1:6)]; 
DATASET.LABEL     = {'YEAR','MONTH','DPROD','REA','RPO','EXOSHOCK'};
DATASET.UNIT      = [0         0        2      1     1   ]; 
DATASET.FIGLABELS = [{'YEAR','MONTH','Oil Production Growth','Global Economic Activity','Oil Price','Supply Shock'}];  
DATASET.MAP = containers.Map(DATASET.LABEL,1:size(DATASET.TSERIES,2));


% VAR specification
%%%%%%%%%%%%%%%%%%%%
 VAR.p      = 12;                                 % Number of Lags
 VAR.irhor  = 48;                                % Impulse Response Horizon
 VAR.select_vars      = {'DPROD','REA','RPO'};
 VAR.vars             = DATASET.TSERIES(:,cell2mat(values(DATASET.MAP,VAR.select_vars))); 
 VAR.MAP              = containers.Map([VAR.select_vars],[1:size(VAR.vars,2)]);
 VAR.proxies          = DATASET.TSERIES(:,cell2mat(values(DATASET.MAP,{'EXOSHOCK'})));
 VAR.DET              = ones(length(VAR.vars),1);

 VAR            = doProxySVAR(VAR); 

% Inference:
%%%%%%%%%%%%
% method  1: Mertens and Ravn (2013) wild bootstrap 
%         2: Montiel-Olea Stock Watson (2016) parametric bootstrap
%         3: Delta Method 
%         4: Montiel-Olea Stock Watson (2016) asy weak IV
%         5: Jentsch and Lunsford Moving Block Bootstrap
%         6: Jentsch and Lunsford Moving Block Bootstrap (adjusted to allow non zero-mean proxies)

nboot     = 5000;         % Number of Bootstrap Samples (Paper does 5000)
clevel    = 68;          % Bootstrap Percentile Shown
BlockSize = floor(5.03*length(DATASET.TSERIES).^0.25); % size of blocks in the MBB bootstrap
seed      = 2;           % seed for random number generator
rng(seed);               % iniate the random number generator

VARci_wildbs   = doProxySVARci(VAR,clevel,1,nboot); 
VARci_mswbs    = doProxySVARci(VAR,clevel,2,nboot);
VARci_delta    = doProxySVARci(VAR,clevel,3);
% VARci_msw_wiv  = doProxySVARci(VAR,clevel,4);
VARci_mbb      = doProxySVARci(VAR,clevel,6,nboot,BlockSize);


shocksize = 20;
VAR.irs   = shocksize*VAR.irs;

legendflag = 1;
shock = 1;
plotdisplay = {'REA','RPO'};

VARci.irsL  = shocksize*VARci_delta.irsL(:,:,shock);
VARci.irsH  = shocksize*VARci_delta.irsH(:,:,shock);
VARci.irsL2 = shocksize*VARci_wildbs.irsL(:,:,shock);
VARci.irsH2 = shocksize*VARci_wildbs.irsH(:,:,shock);
VARci.irsL3 = shocksize*VARci_mswbs.irsL(:,:,shock);
VARci.irsH3 = shocksize*VARci_mswbs.irsH(:,:,shock);
VARci.irsL4 = shocksize*VARci_mbb.irsL(:,:,shock);
VARci.irsH4 = shocksize*VARci_mbb.irsH(:,:,shock);
legendlabels = {'Delta Method','Wild Bootstrap','Par. Bootstrap','Block Bootstrap'};
name = 'Figure4_OIL_shock';
LEGENDPANEL = 2;
xlab_text='horizon (months)';
do_IRS_Figure4m





