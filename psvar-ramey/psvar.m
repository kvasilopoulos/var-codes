% This file estimates the IRF of proxy SVAR for a monetary SVAR
% This program is adapted from replication files of Mertens and Ravn, 2014, 'A
% Reconciliation of SVAR and Narrative Estimates of Tax Mulitpliers', JME

% Note: the number of bootstrap replications is set to 100 for speed (the paper uses 10000)
% In my case, ffr is tax, lip is G, lcpi is Y, lcons is W.
clear all; close all;

% Load the data 

% Load the data
DATA = xlsread('econ214_monetarydat.xlsx',2);
% line 1: Dates
% line 2: log(industry production)
% line 3: unemployment rate
% line 4: log(CPI)
% line 5: federal funds rate
% line 6: log(commodity price)
% line 7: log real durable consumption
% line 8: log real nondurable consumption
% line 9: log real services consumption
% line 10: Romer shock, estimated using full sample
% line 11: Romer shock, full, cumulated
% line 12: Romer shock, estimated 1983-
% line 13: Romer shock, estimated 1983-, cumulated
% line 14: Romer shock, GARCH estimated
% line 15: Romer shock, GARCH estimated, cumulated
% line 16: Gertler-Karadi shock, 3-month ahead fed funds futures
% line 17: Gertler-Karadi shock, 6-month ahead euro dollar futures
% line 18: One-year Treasury bill rate


%  Proxy SVAR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 dates = DATA(:,1);

 VARNA.vars  = DATA(:,[5,2,4,6]);
 VARNA.p = 12;
 [T,n]  = size(VARNA.vars);
 VARNA.DET = ones(T,1);
 VARNA.irhor = 48;

 VARNA.mshocks = DATA(:,10);
 
 % Demean the narrative shocks
  for j=1:size(VARNA.mshocks,2)
     VARNA.mshocks(VARNA.mshocks(:,j)~=0,j)=(VARNA.mshocks(VARNA.mshocks(:,j)~=0,j)-mean(VARNA.mshocks(VARNA.mshocks(:,j)~=0,j)));
  end
  
DAT.TRY = 0.1746; % Average ratio of federal tax revenues to GDP
DAT.GY  = 0.0981; % Average ratio of federal expenditures to GDP

VARNA.mshocksize =  1;
VARNA.gshocksize =  0.01/DAT.GY; % left in because one of the other programs needs it
VARNA = doPVAR(VARNA);
nboot = 500; % Note:  the paper uses 10000
clevel = 95; % Confidence Level
VARNAbs = doPVARbs(VARNA,DAT,nboot,clevel);

% Parameter Estimates
fprintf('thetaG = %f \n',VARNA.thetaG);
fprintf('thetaY = %f \n',VARNA.thetaY);
fprintf('gammaT = %f \n',VARNA.gammaT);
fprintf('zetaT = %f \n',VARNA.zetaT);
fprintf('zetaG = %f \n',VARNA.zetaG);

% Plot Impulse Response to a Monetary Shock
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FIG.axes = [-2 3;-6 4;-2 2;-4 3];

f=figure;
    
    box on
        plot(100*VARNA.irs(:,4),'LineWidth',2,'Color', [0 0 0.5] );
        hold on
        plot(100*VARNAbs.irsH(:,4),'LineWidth',1,'Color',[0 0 0.5],'LineStyle','--' ); 
        hold on
        plot(100*VARNAbs.irsL(:,4),'LineWidth',1,'Color',[0 0 0.5],'LineStyle','--' ); 
        hold on
        axis([0 48 FIG.axes(2,1) FIG.axes(2,2)])
        hline(0,'k-')
        ti=title('Commodity Prices');
        xl=xlabel('months');
        yl=ylabel('percent');
    
        set([xl,yl], 'FontName', 'AvantGarde','FontSize',14);
        set([ti], 'FontName', 'AvantGarde','FontSize',16);
        

f=figure;
    
    box on
        plot(100*VARNA.irs(:,3),'LineWidth',2,'Color', [0 0 0.5] );
        hold on
        plot(100*VARNAbs.irsH(:,3),'LineWidth',1,'Color',[0 0 0.5],'LineStyle','--' ); 
        hold on
        plot(100*VARNAbs.irsL(:,3),'LineWidth',1,'Color',[0 0 0.5],'LineStyle','--' ); 
        hold on
        axis([0 48 FIG.axes(1,1) FIG.axes(1,2)])
        hline(0,'k-')
        ti=title('CPI');
        xl=xlabel('months');
        yl=ylabel('percent');
    
        set([xl,yl], 'FontName', 'AvantGarde','FontSize',14);
        set([ti], 'FontName', 'AvantGarde','FontSize',16); 
 
          
   
    f=figure;
    
    box on
        plot(100*VARNA.irs(:,2),'LineWidth',2,'Color', [0 0 0.5]);
        hold on       
        plot(100*VARNAbs.irsH(:,2),'LineWidth',1,'Color',[0 0 0.5],'LineStyle','--' ); 
        hold on
        plot(100*VARNAbs.irsL(:,2),'LineWidth',1,'Color',[0 0 0.5],'LineStyle','--' ); 
        hold on
        axis([0 48 FIG.axes(4,1) FIG.axes(4,2)])
        hline(0,'k-')
        ti=title('Industrial Production');
        xl=xlabel('months');
        yl=ylabel('percent');
    
        set([xl,yl], 'FontName', 'AvantGarde','FontSize',14);
        set([ti], 'FontName', 'AvantGarde','FontSize',16); 
         
   f=figure;
    
    box on
        plot(VARNA.irs(:,1),'LineWidth',2,'Color', [0 0 0.5]);
        hold on       
        plot(VARNAbs.irsH(:,1),'LineWidth',1,'Color',[0 0 0.5],'LineStyle','--' ); 
        hold on
        plot(VARNAbs.irsL(:,1),'LineWidth',1,'Color',[0 0 0.5],'LineStyle','--' ); 
        hold on
        axis([0 48 FIG.axes(3,1) FIG.axes(3,2)])
        hline(0,'k-')
        ti=title('Fed Fund Rate');
        xl=xlabel('months');
        yl=ylabel('percent');
    
        set([xl,yl], 'FontName', 'AvantGarde','FontSize',14);
        set([ti], 'FontName', 'AvantGarde','FontSize',16); 
        


