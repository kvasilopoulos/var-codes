
%
%
%              SVAR IDENTIFICATION WITH EXTERNAL INSTRUMENTS
%
%
% OLS VAR with bootstrapped confidence bands
%
% written for: "Unsurprising Shocks: Information, Premia and the Monetary
%               Transmission"
%
%
%
% content:      computes responses of US macroeconomic variables to monetary
%               policy shock identified using:
%               (1) cholesky ordering; 
%               (2) Gertler and Karadi (2015) instrument;
%               (3) Orthogonal Proxy (as in Miranda-Agrippino 2016) 
% 
% 
%
%  US economy is a VAR(p)
% =========================
%
% Y(t) = A * Y(t-1) + u(t);     E[u(t)u(t)'] = Sigma;            (1)
%
% u(t) = B0 * e(t)                                               (2)
%
% IDENTIFICATION:
% ---------------
% E[z(t)eM(t)']   =R                                             (3)
% E[z(t)eO(t+j)'] =0 \forall j \neq 0                            (4)
%
% z(t)= external instrument for identification
%
%
%
% miranda (2016) silvia.miranda-agrippino@bankofengland.co.uk
% ----------------------------------------------------------------------- %

clear
clc

addpath([pwd '/datafolder/'])  %mac
addpath([pwd '/subroutines/']) %mac


%plot utils
plotIRFs=true; 
saveCharts=true;

%declare dataset features
%the dataset is a superset of those then used in each VAR below
dataSetSpec.sourceFile          ='FRED-MD_2015SM'; %Mc-Cracken & Ng data

dataSetSpec.useStoredList       =false;            %char with stored data list -- see buildDataSet.m
dataSetSpec.dataList            ={'INDPRO';'UNRATE';'CPIAUCSL';'CRBPI';'FEDFUNDS';'GS1'}; %if stored list leave blank

dataSetSpec.beginSet            =datenum(1969,1,1);
dataSetSpec.endSet              =datenum(2014,12,1);

dataSetSpec.plotData            =false;            %plots individual(!) charts 
dataSetSpec.interpolateMissing  =false;

%load data and build dataset
dataStructure                   =buildDataSet(dataSetSpec);

%common model specification
modelSpec.nLags                 =12;            %VAR lags
modelSpec.nHorizons             =48;            %max horizon for IRFs
modelSpec.cLevel                =[68 90];       %size of error bands
modelSpec.bootSize              =1000;%0;         %size of bootstrap sample for error bands

%-1-CHOLESKY--------------------------------------------------------------%

modelSpec.identification        ='CHOL';

%select relevant data
cholSet                         ={'INDPRO';'UNRATE';'CPIAUCSL';'CRBPI';'FEDFUNDS'};
[~,dataSelection]               =ismember(cholSet,dataStructure.varname);

modelSpec.dataSelection         =dataSelection;

%declare shock variable and shock size
shockVar                        ='FEDFUNDS';
modelSpec.shockVar              =ismember(cholSet,shockVar);
modelSpec.shockSize             =1*double(ismember(cholSet,shockVar));

%compute responses
irfCholesky                     =computeIRFs(dataStructure,modelSpec);

%-2-PSVAR GK--------------------------------------------------------------%

modelSpec.identification        ='PSVAR';
modelSpec.selectedInstrument    ='FF4GK'; 

%select relevant data
psvarSet                        ={'INDPRO';'UNRATE';'CPIAUCSL';'CRBPI';'GS1'};
[~,dataSelection]               =ismember(psvarSet,dataStructure.varname);

modelSpec.dataSelection         =dataSelection;

%declare shock variable and shock size
shockVar                        ='GS1';
modelSpec.shockVar              =ismember(psvarSet,shockVar);
modelSpec.shockSize             =1*double(ismember(psvarSet,shockVar));

%compute responses
irfSVARIV_1                     =computeIRFs(dataStructure,modelSpec);

%-3-PSVAR ROBUST----------------------------------------------------------%

modelSpec.identification        ='PSVAR';
modelSpec.selectedInstrument    ='FF4MA'; 


%select relevant data
psvarSet                        ={'INDPRO';'UNRATE';'CPIAUCSL';'CRBPI';'GS1'};
[~,dataSelection]               =ismember(psvarSet,dataStructure.varname);

modelSpec.dataSelection         =dataSelection;


%declare shock variable and shock size
shockVar                        ='GS1';
modelSpec.shockVar              =ismember(psvarSet,shockVar);
modelSpec.shockSize             =1*double(ismember(psvarSet,shockVar));

%compute responses
irfSVARIV_2                     =computeIRFs(dataStructure,modelSpec);






% -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  %
%SAVE OUTPUT
runid =['IRFs_' datestr(dataSetSpec.beginSet,'YY') datestr(dataSetSpec.endSet,'YY') ]; 

save(runid);
% -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  %




%%
%-plot irfs---------------------------------------------------------------%
Ccolor  =[ 1 .6 .2]; 
GKcolor =[.0 .0 .4]; 
Scolor  =[ 1 .0 .0];    bandFillColor1 =[.85 .85 .85];  bandFillColor2 =[.7 .7 .7]; 

%plot labels
varname                 =dataStructure.varLongName; 
varname{end-1}          ='Policy Rate';
varname{end}            ='Policy Rate';


figure; 
plotColumns =3; pln=1; 
nHorizon    =modelSpec.nHorizons; 
n           =length(varname)-1;

for j=1:n
    
    subplot(ceil(n/plotColumns),plotColumns,pln)
    
    hold on
    
    %bands 90
    fill([0:1:nHorizon, fliplr(0:1:nHorizon)],...
        [squeeze(irfSVARIV_2.irfs_u(j,:,2)) fliplr(squeeze(irfSVARIV_2.irfs_l(j,:,2)))],...
        bandFillColor1,'EdgeColor','none');

    %bands 68
    fill([0:1:nHorizon, fliplr(0:1:nHorizon)],...
        [squeeze(irfSVARIV_2.irfs_u(j,:,1)) fliplr(squeeze(irfSVARIV_2.irfs_l(j,:,1)))],...
        bandFillColor2,'EdgeColor','none');
    
    
    %zero line
    plot(0:nHorizon,zeros(size(0:nHorizon)),'k')    
    hold on
    
    %irfs
    plot(0:nHorizon,irfCholesky.irfs(j,:),'--','LineWidth',1.5,'color',Ccolor)

    %irfs 
    plot(0:nHorizon,irfSVARIV_1.irfs(j,:),'-.','LineWidth',1.5,'color',GKcolor)

    %irfs
    plot(0:nHorizon,irfSVARIV_2.irfs(j,:),'-','LineWidth',1.5,'color',Scolor)

    
    hold off; axis tight

    grid on; xlim([0 nHorizon]);
    set(gca,'XTick',0:6:nHorizon,'XTickLabel',cellstr(num2str((0:6:nHorizon)')),'layer','top')
    title(varname{j},'FontSize',9,'FontWeight','normal')

    if j==1
    
        ylabel('% points','FontSize',9)
    end
    
    if j==n
        
        xlabel('horizon','FontSize',9)

        pln =pln+1;
        pl  =subplot(ceil(n/plotColumns),plotColumns,pln);
        
        %nicer legend
        hold on
        p0 =plot(1:nHorizon+1,nan(size(1:nHorizon+1)),'s','MarkerFaceColor',Ccolor,'MarkerEdgeColor',Ccolor,'MarkerSize',7);
        p1 =plot(1:nHorizon+1,nan(size(1:nHorizon+1)),'s','MarkerFaceColor',GKcolor,'MarkerEdgeColor',GKcolor,'MarkerSize',7);
        p2 =plot(1:nHorizon+1,nan(size(1:nHorizon+1)),'s','MarkerFaceColor',Scolor,'MarkerEdgeColor',Scolor,'MarkerSize',7);
        hold off
        
        legend([p0 p1 p2],{'Cholesky';'Raw Market Surprise';'Orthogonal Proxy'},'Location','NorthWest','FontSize',10);
        legend boxoff
        set(pl,'vis','off'); 

    end

    pln=pln+1;
end


set(gcf,'PaperUnits','centimeters','PaperSize',[22 10]) %[x y]
set(gcf,'PaperPosition',[-2 0 25 10]) %[left bottom width height]
print(gcf,'-dpdf',[ runid '.pdf']); 
        

