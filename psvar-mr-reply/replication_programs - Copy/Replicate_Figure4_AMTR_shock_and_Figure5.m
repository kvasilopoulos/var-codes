clear all; close all; 
addpath('auxiliary files')

% Create DATASET
%%%%%%%%%%%%%%%%
AMTR                = xlsread('datasets/DATA_MMO.xlsx','AMTR (Figure I)'); 
CONTROLS            = xlsread('datasets/DATA_MMO.xlsx','CONTROLS');
AINC                = xlsread('datasets/DATA_MMO.xlsx','LOG AVG INCOME');
NARR                = xlsread('datasets/DATA_MMO.xlsx','Narrative Shocks (Table IV)'); 
PROX                = NARR(:,6);
PROX(isnan(PROX))   = 0;

DATASET.TSERIES = [AMTR(:,1) -log(1-AMTR(:,2)) CONTROLS(:,2:8) AINC(:,2) PROX]; 
DATASET.LABEL   = {'DATES','LNAMTR1','RGDP','UNRATE','INFL','FFR','GOV_TOT','STPRICE','DEBT','AINC','PROX'};
DATASET.UNIT    = [0         1         1        2      2      2      1         1        1       1     2 ]; 
DATASET.FIGLABELS = [{'year' ,'1/(1-AMTR) All Tax Units','Output','Unemployment Rate','Inflation', ...
                    'Federal Funds Rate','Government Spending','Real Stock Prices', 'Debt','Reported Income','AMTR proxy'}];  
DATASET.MAP = containers.Map(DATASET.LABEL,1:size(DATASET.TSERIES,2));

% VAR specification
%%%%%%%%%%%%%%%%%%%%
VAR.sb              = 1946;  % Sample start year
VAR.se              = 2012;  % Sample end year
TSERIES             = DATASET.TSERIES((DATASET.TSERIES(:,1)>=VAR.sb)&(DATASET.TSERIES(:,1)<=VAR.se),:);
VAR.dates           = TSERIES(:,1);
VAR.DET             = [VAR.dates==1949 VAR.dates==2008 ones(size(TSERIES,1),1)]; % Deterministic Terms (Put constant last)
VAR.p               = 2; % Lag length
VAR.NWlags          = 8; % Newey-West Lags
VAR.irhor           = 6; % Length of IR horizon 
VAR.select_variables= {'LNAMTR1','AINC','RGDP','UNRATE','INFL','FFR','GOV_TOT','STPRICE','DEBT'};
VAR.vars            = TSERIES(:,cell2mat(values(DATASET.MAP,VAR.select_variables)));
VAR.proxies         = TSERIES(:,cell2mat(values(DATASET.MAP,{'PROX'})));
VAR.MAP             = containers.Map([VAR.select_variables],[1:size(VAR.vars,2)]);   

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
VARci_mbb      = doProxySVARci(VAR,clevel,5,nboot,BlockSize);

shocksize = 1;
VAR.irs   = shocksize*VAR.irs;

legendflag = 1;
shock = 1;
plotdisplay = {'AINC','RGDP'};
VARci.irsL  = shocksize*VARci_delta.irsL(:,:,shock);
VARci.irsH  = shocksize*VARci_delta.irsH(:,:,shock);
VARci.irsL2 = shocksize*VARci_wildbs.irsL(:,:,shock);
VARci.irsH2 = shocksize*VARci_wildbs.irsH(:,:,shock);
VARci.irsL3 = shocksize*VARci_mswbs.irsL(:,:,shock);
VARci.irsH3 = shocksize*VARci_mswbs.irsH(:,:,shock);
VARci.irsL4 = shocksize*VARci_mbb.irsL(:,:,shock);
VARci.irsH4 = shocksize*VARci_mbb.irsH(:,:,shock);
legendlabels = {'Delta Method','Wild Bootstrap','Par. Bootstrap','Block Bootstrap'};
name = 'Figure4_AMTR_shock';
LEGENDPANEL = 1;
xlab_text='horizon (years)';
do_IRS_Figure4m


figure 
hold on 
box on
histogram(VARci_mbb.Waldstat_bs)
%ti=title('MBB F-stat Distribution');
xl=xlabel('value');
yl=ylabel('frequency');
%vline(SVARIVci.Waldstat)
str=strcat('figures/',name,'_AER_','Fdistr');
                    saveas(gcf,str,'epsc');
%set([ti], 'FontName', 'AvantGarde','FontSize',16);  
set([xl,yl], 'FontName', 'AvantGarde','FontSize',14);

str=strcat('figures/','Figure5');
            saveas(gcf,str,'epsc');
