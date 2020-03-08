clear all; close all;

% Quarterly
DATASET.TSERIES=xlsread('MR_AER_DATASET.xlsx','Quarterly');
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

DATASET.FIGLABELS{1,1}  = 'Quarter';                                   
DATASET.FIGLABELS{2,1}  = 'Narrative PI tax shock';                                   
DATASET.FIGLABELS{3,1}  = 'Narrative CI tax shock';  
DATASET.FIGLABELS{4,1}  = 'PI tax shock proxy';  
DATASET.FIGLABELS{5,1}  = 'CI tax shock proxy';  
DATASET.FIGLABELS{6,1}  = 'Average Personal Income Tax Rate';   
DATASET.FIGLABELS{7,1}  = 'Average Corporate Income Tax Rate';  
DATASET.FIGLABELS{8,1}  = 'Personal Income Tax Base';           
DATASET.FIGLABELS{9,1}  = 'Corporate Income Tax Base';  
DATASET.FIGLABELS{10,1} = 'Government Purchases';  
DATASET.FIGLABELS{11,1} = 'Output';                             
DATASET.FIGLABELS{12,1} = 'Government Debt';               
DATASET.FIGLABELS{13,1} = 'Federal Funds Rate';               
DATASET.FIGLABELS{14,1} = 'Price Level';               
DATASET.FIGLABELS{15,1} = 'Nonborrowed Reserves';               
DATASET.FIGLABELS{16,1} = 'Employment/Population';               
DATASET.FIGLABELS{17,1} = 'Hours Per Worker';               
DATASET.FIGLABELS{18,1} = 'Labor Force/Population';               
DATASET.FIGLABELS{19,1} = 'Consumption (Nondurables and Services)';               
DATASET.FIGLABELS{20,1} = 'Durable Good Purchases';               
DATASET.FIGLABELS{21,1} = 'Nonresidential Investment'; 
DATASET.FIGLABELS{22,1} = 'Residential Investment';               
DATASET.FIGLABELS{23,1} = 'Unemployment Rate';               
DATASET.FIGLABELS{24,1} = 'Inflation';               
DATASET.FIGLABELS{25,1} = 'Personal Income Tax Revenues';           
DATASET.FIGLABELS{26,1} = 'Corporate Income Tax Revenues';
save('DATASET','DATASET');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Annual
DATASET_A.TSERIES =xlsread('MR_AER_DATASET.xlsx','Annual');
DATASET_A.LABEL   = {'DATES','T_PI','m_PI','APITR','MPITR','PITB','GOV','RGDP','DEBT'};
DATASET_A.VALUE   = [  1,       2,     3,     4,    5 ,    6 ,     7 ,     8  ,   9 ];
DATASET_A.UNIT    = [  0,       2,     2,     2,    2 ,    1 ,     1 ,     1  ,   1 ];

DATASET_A.MAP = containers.Map(DATASET_A.LABEL,DATASET_A.VALUE);

DATASET_A.FIGLABELS{1,1}  = 'Quarter';                                   
DATASET_A.FIGLABELS{2,1}  = 'Narrative PI tax shock';                                   
DATASET_A.FIGLABELS{3,1}  = 'PI tax shock proxy';  
DATASET_A.FIGLABELS{4,1}  = 'Average Personal Income Tax Rate';   
DATASET_A.FIGLABELS{5,1}  = 'Marginal Personal Income Tax Rate';  
DATASET_A.FIGLABELS{6,1}  = 'Personal Income Tax Base';           
DATASET_A.FIGLABELS{7,1}  = 'Government Purchases';  
DATASET_A.FIGLABELS{8,1}  = 'Output';                             
DATASET_A.FIGLABELS{9,1}  = 'Government Debt';               
save('DATASET_A','DATASET_A');
