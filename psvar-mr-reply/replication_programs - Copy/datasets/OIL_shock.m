clear all; close all;
addpath('auxiliary files')

mdata = xlsread('Kilian_data.xlsx');
mdata(isnan(mdata(:,6)),6)=0;
mdata(:,6) = sum(lagmatrix(mdata(:,6),[0 1 2]),2);
mdata = mdata(2:end,:);

load kilianoilshock
monthlyshock = data.exoshock2(2:end-9);
DATASET.TSERIES   = [mdata(:,1:5)]; 
DATASET.LABEL     = {'YEAR','MONTH','DPROD','REA','RPO'};
DATASET.UNIT      = [0         0        2      1     1   ]; 
DATASET.FIGLABELS = [{'YEAR','MONTH','Oil Production Growth','Global Economic Activity','Oil Price'}];  
DATASET.MAP = containers.Map(DATASET.LABEL,1:size(DATASET.TSERIES,2));

nboot  = 5000;          % Number of Bootstrap Samples (equals 10000 in the paper)
clevel = 68;            % Bootstrap Percentile Shown
BlockSize = floor(5.03*length(DATASET.TSERIES).^0.25);      % size of blocks in the MBB bootstrap

% VAR specification
%%%%%%%%%%%%%%%%%%%%
 VAR.p      = 12;                                 % Number of Lags
 VAR.irhor  = 48;                                % Impulse Response Horizon

 VAR.select_vars      = {'DPROD','REA','RPO'};
 VAR.vars             = DATASET.TSERIES(:,cell2mat(values(DATASET.MAP,VAR.select_vars)));
 
 VAR.MAP              = containers.Map([VAR.select_vars],[1:size(VAR.vars,2)]);
 VAR.proxies          = monthlyshock;
 %VAR.proxies  = mdata(:,6);
 VAR.proxies(isnan(VAR.proxies))=0;
 VAR.DET    = ones(length(VAR.vars),1);
 VAR.NWlags = floor(4*(((length(VAR.vars)-VAR.p)/100)^(2/9)));
 VAR            = doProxySVARwF(VAR); 

% Inference:
% method  1: Mertens and Ravn (2013) wild bootstrap 
%         2: Montiel-Olea Stock Watson (2016) bootstrap
%         3: Delta Method 
%         4: Montiel-Olea Stock Watson (2016) asy weak IV
%         5: Jentsch and Lunsford Moving Block Bootstrap

[VARci_wildbs,VARci_wildbs_cent]  = doProxySVARci(VAR,clevel,1,nboot); 
[VARci_mswbs,VARci_mswbs_cent]    = doProxySVARci(VAR,clevel,2,nboot);
 VARci_delta    = doProxySVARci(VAR,clevel,3);
%VARci_msw_wiv  = doProxySVARci(VAR,clevel,4);
[VARci_mbb,VARci_mbb_cent]         = doProxySVARci(VAR,clevel,6,nboot,BlockSize);

shocksize = 20;
VAR.irs   = shocksize*VAR.irs;

legendflag = 1;
shock = 1;
plotdisplay = {'DPROD','REA','RPO'};
VARci.irsL  = shocksize*VARci_delta.irsL(:,:,shock);
VARci.irsH  = shocksize*VARci_delta.irsH(:,:,shock);
VARci.irsL2 = shocksize*VARci_wildbs.irsL(:,:,shock);
VARci.irsH2 = shocksize*VARci_wildbs.irsH(:,:,shock);
VARci.irsL3 = shocksize*VARci_mswbs.irsL(:,:,shock);
VARci.irsH3 = shocksize*VARci_mswbs.irsH(:,:,shock);
VARci.irsL4 = shocksize*VARci_mbb.irsL(:,:,shock);
VARci.irsH4 = shocksize*VARci_mbb.irsH(:,:,shock);

% VARci.irsL  = shocksize*VARci_delta.irsL(:,:,shock);
% VARci.irsH  = shocksize*VARci_delta.irsH(:,:,shock);
% VARci.irsL2 = shocksize*VARci_wildbs.irsLhall(:,:,shock);
% VARci.irsH2 = shocksize*VARci_wildbs.irsHhall(:,:,shock);
% VARci.irsL3 = shocksize*VARci_mswbs.irsLhall(:,:,shock);
% VARci.irsH3 = shocksize*VARci_mswbs.irsHhall(:,:,shock);
% VARci.irsL4 = shocksize*VARci_mbb.irsLhall(:,:,shock);
% VARci.irsH4 = shocksize*VARci_mbb.irsHhall(:,:,shock);
legendlabels = {'Delta Method','Wild Bootstrap','Par. Bootstrap','Block Bootstrap'};
name = 'OIL_shock';
SAVE = 1;
LEGENDPANEL = 3;
do_IRS_Figure_comp_monthly



figure 
hold on 
box on
histogram(VARci_mbb.Waldstat_bs)
ti=title('MBB F-stat Distribution');
xl=xlabel('value');
yl=ylabel('frequency');
%vline(SVARIVci.Waldstat)
str=strcat('figures/',name,'_Fdistr');
                    saveas(gcf,str,'epsc');
set([ti], 'FontName', 'AvantGarde','FontSize',16);  
set([xl,yl], 'FontName', 'AvantGarde','FontSize',14);
sum(VARci_mbb.Waldstat_bs<[10 3.84],1)/length(VARci_mbb.Waldstat_bs)
