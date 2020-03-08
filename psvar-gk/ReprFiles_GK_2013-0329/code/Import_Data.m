%Imports time-series data 

% VAR
DATASET_VAR.TSERIES=csvread('../data/VAR_data.csv',1,0);          %Reads data, skips names
DATASET_VAR.LABEL   = {'YEAR','MONTH','LCPI','LIP','FF','GS1','GS2','CM5YR','CM10YR','CM5F5','EBP'};
DATASET_VAR.VALUE   = [   1,     2,      3,     4,   5 ,  6  ,  7  ,   8   ,    9   ,  10   ,  11];

DATASET_VAR.LABEL   = [DATASET_VAR.LABEL,{'MORTG_SPREAD','CP3M_SPREAD', 'FF_EXP1YR'}];
DATASET_VAR.VALUE   = [DATASET_VAR.VALUE,        12     ,   13        ,      14 ];
 
% % 
% % Creates a mapping between labels and values
DATASET_VAR.MAP = containers.Map(DATASET_VAR.LABEL,DATASET_VAR.VALUE);
% 
%Factors
DATASET_FACTORS.TSERIES=csvread('../data/factor_data.csv',1,0);          %Reads data, skips names
DATASET_FACTORS.LABEL   = {'YEAR','MONTH','MP1_TC','FF4_TC','ED2_TC','ED3_TC','ED4_TC'};
DATASET_FACTORS.VALUE   = [   1,     2   ,  3    ,    4   ,    5   ,   6    ,  7     ];

% % 
% % Creates a mapping between labels and values
DATASET_FACTORS.MAP = containers.Map(DATASET_FACTORS.LABEL,DATASET_FACTORS.VALUE);
% 
save('../data/DATASET.mat','DATASET_VAR','DATASET_FACTORS');
