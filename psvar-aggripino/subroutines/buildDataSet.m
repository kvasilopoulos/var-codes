

function res=buildDataSet(dataSetSpec)
%
% loads up raw data and applies required transformations
% 
% inputs (structure)
% sourceFile         = mat file name (see content below)
% useStoredList      = can be either *false* or a string identifier;
%                      stored list include:
%                      BGMcee:    EMPL, CPI, PCMDY, FFR, MTWO, TOTR, NBR
%                      BGMsmall:  EMPL, CPI, FFR
%                      BGMmedium: EMPL, CPI, PCMDY, PINC, CONS, IPROD, CAPU,
%                                 URATE, HSTART, PPI, PCED, WAGE, FFR, MONE, 
%                                 MTWO, TOTR, NBR, S&P, 10YBY, EXR
%                      FULL:      uses all
% seriesList         = if useStoredList is false, allows to define a new 
%                      one, comes in the form of a cell string {'';'';''}
% beginSet           = date
% endSet             = date
% takeFirstDiff      = true/false, applies to all
% plotData           = true/false if true plots original and transformed
% interpolateMissing = true/false, fills up missing values with centered MA
%                      does not apply to beginning and end of series
%
% output (structure)
% data        = [Txn] matrix of transformed data
% dates       = [Tx1] vector of reference dates
% varname     = {1xn} cell of data identifiers
% varLongName = {1xn} cell of variables names
% 
% miranda 2015 smirandaagrippino@london.edu


listCEE  = {'PAYEMS';'CPIAUCSL';'PPICMM';'FEDFUNDS';'M2SL';'TOTRESNS';'NONBORRES'};
listBGMs = {'PAYEMS';'CPIAUCSL';'FEDFUNDS'};
listBGMm = {'PAYEMS';'CPIAUCSL';'PPICMM';'W875RX1';'DPCERA3M086SBEA';...
    'INDPRO';'CAPUTLB00004S';'UNRATE';'HOUST';'PPIFGS';'PCEPI';'CES3000000008';...
    'FEDFUNDS';'M1SL';'M2SL';'TOTRESNS';'NONBORRES';'S&P 500';'GS10';'TWEXMMTH'};


%load raw data from file
load(dataSetSpec.sourceFile);

%sourceFile is in mat format; contains:
% data            = [TxN] matrix of raw data
% dates           = [Tx1] vector of dates
% logTransform    = [1xN] logical for log transformations
% dataName        = {1xN} cell of data identifiers
% dataDescription = {1xN} cell of variables names

%detect data list
if ischar(dataSetSpec.useStoredList)
    switch dataSetSpec.useStoredList
        
        case 'CEE'
            dataList = listCEE;
            
        case 'BGMs'
            dataList = listBGMs;
            
        case 'BGMm'
            dataList = listBGMm;
            
        case 'FULL'
            dataList = dataName;
    end
    
elseif ~dataSetSpec.useStoredList
    
    dataList=dataSetSpec.dataList;
end

[~,dataSelect]=ismember(dataList,dataName); dataSelect(dataSelect==0)=[];

%trim relevant items
dataRaw         = data(:,dataSelect);
datesRaw        = dates;
dataName        = dataName(dataSelect);
dataDescription = dataDescription(dataSelect);
logTransform    = logTransform(dataSelect);

%correct unit for Non-Borrowed Reserves
if ismember('NONBORRES',dataList)
    
    dataRaw(:,ismember(dataList,'NONBORRES'))=dataRaw(:,ismember(dataList,'NONBORRES'))./1000;
    
end


data=dataRaw; dates=datesRaw;

%take logs
logTransform(and(any(data<0),logTransform))=false; %remove transform if <0
data(:,logTransform)=log(data(:,logTransform))*100;


%fill up NaNs
if dataSetSpec.interpolateMissing
    
    nanOpt.method  = 2;
    nanOpt.winsize = 1;
    
    data=removeNaNs(data,nanOpt);
    
end

%select relevant time span
timeSelect= dates >= dataSetSpec.beginSet & dates <= dataSetSpec.endSet;

data  = data(timeSelect,:);
dates = dates(timeSelect);

%remove rigged edges
edges = any(isnan(data),2);

data  = data(~edges,:);
dates = dates(~edges);

%plot data
if dataSetSpec.plotData
    
    for j=1:size(data,2)
        
        figure;
        subplot(2,1,1)
        plot(datesRaw,dataRaw(:,j)); axis tight; grid on; dateaxis('x',10);
        set(gca,'FontSize',9)
        title([dataDescription{j}, ' all obs'],'FontSize',11,'FontWeight','normal')
        %
        subplot(2,1,2)
        plot(dates,data(:,j)); axis tight; grid on; dateaxis('x',10);
        set(gca,'FontSize',9)
        title('series over selected sample','FontWeight','normal')
        %
        pause;
        close(gcf)

        
    end
    
end

%load output
res.data        = data;
res.dates       = dates;
res.varname     = dataName;
res.varLongName = dataDescription;