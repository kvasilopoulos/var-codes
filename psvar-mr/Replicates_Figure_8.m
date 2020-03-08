clear all; close all; addpath('auxfiles');
 
load DATASET_A; 

nboot  = 1000;  % Number of Bootstrap Samples (equals 10000 in the paper)
clevel = 95;    % Bootstrap Percentile Shown

% VAR with Average PI  Tax Rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 VAR1.p      = 2;                                 % Number of Lags
 VAR1.irhor  = 5;                                % Impulse Response Horizon
 VAR1.select_vars      = {'APITR','PITB','GOV','RGDP','DEBT'};
 VAR1.vars             = DATASET_A.TSERIES(:,cell2mat(values(DATASET_A.MAP,VAR1.select_vars)));
 VAR1.MAP              = containers.Map([VAR1.select_vars],[1:size(VAR1.vars,2)]);
 VAR1.proxies          = DATASET_A.TSERIES(:,cell2mat(values(DATASET_A.MAP,{'m_PI'})));
 VAR1      = doProxySVAR_single(VAR1,DATASET_A);
 VAR1bs    = doProxySVARbootstrap_single(VAR1,nboot,clevel,DATASET_A);

% VAR with Marginal PI  Tax Rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 VAR2.p      = 2;                                 % Number of Lags
 VAR2.irhor  = 5;                                 % Impulse Response Horizon
 VAR2.select_vars      = {'MPITR','PITB','GOV','RGDP','DEBT'};
 VAR2.vars             = DATASET_A.TSERIES(:,cell2mat(values(DATASET_A.MAP,VAR2.select_vars)));
 VAR2.MAP              = containers.Map([VAR2.select_vars],[1:size(VAR2.vars,2)]);
 VAR2.proxies          = DATASET_A.TSERIES(:,cell2mat(values(DATASET_A.MAP,{'m_PI'})));
 VAR2      = doProxySVAR_single(VAR2,DATASET_A);
 VAR2bs    = doProxySVARbootstrap_single(VAR2,nboot,clevel,DATASET_A);
 
% Plot Figure 
%%%%%%%%%%%%%
FIG.axes(:,cell2mat(values(VAR1.MAP,{'RGDP'})))   = [-1;4.5];
FIG.axes(:,cell2mat(values(VAR1.MAP,{'APITR'})))  = [-1.5;1];
 
plotdisplay = {'APITR','RGDP'};

display1= cell2mat(values(VAR1.MAP,plotdisplay));
display2= cell2mat(values(VAR1.MAP,plotdisplay));
for nvar = 1:length(display1)
                  
        f=figure;    
        box on

            p1=plot(VAR1.irs(:,display1(nvar)),'-','MarkerSize',4,'LineWidth',2,'Color', [0 0 0.5]);
            hold on
            plot(VAR1bs.irsH(:,display1(nvar)),'LineWidth',1,'Color', [0 0 0.5],'LineStyle','--');
            hold on
            plot(VAR1bs.irsL(:,display1(nvar)),'LineWidth',1,'Color', [0 0 0.5],'LineStyle','--');
            hold on
            
            p2=plot(VAR2.irs(:,display2(nvar)),'-d','MarkerSize',4,'LineWidth',2,'Color', [0.9 0 0]);
            hold on
            plot(VAR2bs.irsH(:,display2(nvar)),'LineWidth',1,'Color', [0.9 0 0],'LineStyle','-.');
            hold on
            plot(VAR2bs.irsL(:,display2(nvar)),'LineWidth',1,'Color', [0.9 0 0],'LineStyle','-.');

            axis([0.75 VAR1.irhor FIG.axes(1,cell2mat(values(VAR1.MAP,{plotdisplay{nvar}}))) FIG.axes(2,cell2mat(values(VAR1.MAP,{plotdisplay{nvar}})))]);
            hline(0,'k-')
            ti=title( DATASET_A.FIGLABELS{cell2mat(values(DATASET_A.MAP,{plotdisplay{nvar}}))});
            xl=xlabel('years');
             set(gca,'XTick',[1 2 3 4 5]);
            if cell2mat(values(VAR1.MAP,{'APITR'}))==display1(nvar)
            l=legend([p1,p2],'Average Tax Rate','Marginal Tax Rate');
            set([l], 'FontName', 'AvantGarde','FontSize',14,'Location','NorthWest');
            ti=title('Personal Income Tax Rate');
            end
            
            if DATASET_A.UNIT(cell2mat(values(DATASET_A.MAP,{plotdisplay{nvar}})))==1
            yl=ylabel('percent');
            else 
            yl=ylabel('percentage points'); 
            end
               
            set([xl,yl], 'FontName', 'AvantGarde','FontSize',14);
            set([ti], 'FontName', 'AvantGarde','FontSize',16);
          
end

