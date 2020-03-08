% Replication files for Gertler, M. and P. Karadi, 2014, Monetary Policy
% Surprises, Credit Costs and Economic Activity, AEJMacro
%
% The code is based on the code of Mertens-Ravn, (AER, 2013) external
% instrument SVAR code, modified to 
%1) handle different length of VAR and factors
%2) Multiple instruments to explain the same variable
%3) Calculates real rate, term premia and excess premia responses to shocks

clear all; close all; clc;

%tic;

%Import and load time-series data
Import_Data; 
load ../data/DATASET DATASET_VAR DATASET_FACTORS;

nboot  = 1000;  % Number of Bootstrap Samples (equals 10000 in the paper)
clevel = 95;    % Bootstrap Percentile Shown
VAR.fontsize        =   14; %fontsize in figures

%Setting up possible specifications
smpl_min_VAR_vec    =   [1979 7]; %year month    
smpl_max_VAR_vec    =   [2012 6; 2008 6]; %       
monpol_vars_cell     =   {'FF','GS1','GS2'};                                %Monetary policy variables  
monpol_vars_label_cell =  {'Fed Funds Rate','1 year rate','2 year rate'};   %Names of the monetary policy variables    
monpol_vars_label_short_cell=  {'FF','1YR','2YR'};                          %Short names of the monetary policy variables    

factors1_cell_FF     =   {'MP1_TC'};    %First factor for monpol variable FF
factors1_label_cell_FF=   {'MP1_TC'};   %First factor name
factors2_cell_FF        =   {''};       %Further factors
factors2_label_cell_FF  =   {''};       %Names of further factors;

factors1_cell_GS1     =   {'FF4_TC'};   
factors1_label_cell_GS1=   {'FF4_TC'};
factors2_cell_GS1        =   {''};    
factors2_label_cell_GS1  =   {''};      

factors1_cell_GS2     =   {'MP1_TC'};  
factors1_label_cell_GS2=   {'MP1_TC'}; 
factors2_cell_GS2        =   {'FF4_TC','ED2_TC','ED3_TC','ED4_TC'}; 
factors2_label_cell_GS2  =   {'FF4_TC','ED2_TC','ED3_TC','ED4_TC'};

smpl_min_factors_vec =   [ones(1,1)*[1991 1]];          %Starting date for all of the factors
smpl_max_factors_vec =   [ones(1,1)*[2012 6]];          %Finishing date for all of the factors
figure_name         =     'ALL';
spreads_cell        =    {'EBP','MORTG_SPREAD','CP3M_SPREAD'};  
spreads_label_cell  =    {'Excess Bond Premium','Mortgage spread','Commercial Paper spread (3 months)'}; %
spreads_label_short_cell  =    {'EBP','MORTG.','CP3M'}; 
no_spread_single    =   [1];                            % run regressions with a particular spread
switch_extern       =   1;                              % run VARs with particular variables added one-by-one (time consuming)
extern_vars_cell    =    {'FF','GS1','GS2','CM5YR','CM10YR','CM5F5','FF_EXP1YR'};  %(FF needs to be the first)
extern_vars_label_cell   =   {'Federal Funds Rate','1 year rate','2 year rate','5 year rate','10 year rate','5x5 forward','1 year expectations (FF)'};
extern_vars_label_short_cell   =   {'FF','1YR','2YR','5YR','10YR','5F5','EXP1YR(FF)'}; 
extern_vars_smpl_min_vec =  [ones(6,1)*[1979 7];[1983 3]]; %
term_spreads_cell   =   {'GS1','GS2','CM5YR','CM10YR'}; %Variables which term-spreads are plotted
term_spreads_label_cell   =   {'1 year rates','2 year rates','5 year rates','10 year rates'};  %
term_spreads_label_short_cell   =   {'1YR','2YR','5YR','10YR'};  %
VAR.n_ts_1          =   4;                              %The number of term spreads variables plotted before the excess_returns
term_spreads_init   =   [1 1 1 1];    
term_spreads_matur  =   [12 24 60 120];                 %Length of term spread maturities in months
excess_return_cell  =   {'EBP','MORTG_SPREAD'};         %Variables whose excess returns are calculated
excess_return_label_cell  =   {'Corporate bond','Mortgage'};
excess_return_label_short_cell  =   {'EBP','MORTG.'};
excess_return_matur_cell  =   {'CM5YR','CM10YR'};       %Safe returns used for excess return calculations
real_rates_cell   =   {'GS2','CM5YR','CM10YR'};         %Variables whose real rate is calculated
real_rates_label_cell   =   {'2 year rate','5 year rate','10 year rate'};  
real_rates_label_short_cell   =   {'2YR','5YR','10YR'}; 
real_rates_init  =   [1 1 1];   %
real_rates_matur  =   [24 60 120];                      %Maturity of yields whose real rate is calculated
bkeven_label_cell   =   {'2 year breakeven inflation rate','5 year breakeven inflation rate','10 year breakeven inflation rate'};  
bkeven_label_short_cell   =   {'2YR','5YR','10YR'};  
VAR.switch_exp  =   1;                                  %0 or 1 to calculate expectation based VAR
exp_rates_cell  =   {'FF_EXP1YR','FF_EXP1YR','FF_EXP1YR','FF_EXP1YR'};
exp_rates_label_cell  =   {'1 year rate (FF)','2 year rate (FF)','5 year rate (FF)','10 year rate (FF)'};  %
exp_rates_label_short_cell  =   {'1YR','2YR','5YR','10YR'};  
exp_rates_comp_cell = {'GS1','GS2','CM5YR','CM10YR'};
exp_rates_matur =   [12 24 60 120];    

VAR.irhor  = 48;                                        % Impulse Response Horizon
VAR.p      = 12;                                        % VAR lag length

no_smpl_min_VAR     =   size(smpl_min_VAR_vec,1);
no_smpl_max_VAR     =   size(smpl_max_VAR_vec,1);
no_monpol_vars      =   length(monpol_vars_cell);
no_extern_vars      =   length(extern_vars_cell);
no_term_spreads     =   length(term_spreads_cell);
no_excess_return    =   length(excess_return_cell);
no_real_rates       =   length(real_rates_cell);
no_exp_rates        =   length(exp_rates_cell);
no_spreads_all      =   length(spreads_cell);
no_spreads          =   length(no_spread_single);

ii=cell(no_monpol_vars,1);
        %Counting the figures
for ii_monpol=1:no_monpol_vars
        if strcmp(monpol_vars_cell{ii_monpol},'FF')
            VAR.monpol_FF = 'yes';
        else
            VAR.monpol_FF = 'no';
        end;
        ii{ii_monpol} = 1000*(ii_monpol-1)+1;
        eval(['factors1_cell=factors1_cell_' monpol_vars_cell{ii_monpol} ';']);
        eval(['factors2_cell=factors2_cell_' monpol_vars_cell{ii_monpol} ';']);
        eval(['factors1_label_cell=factors1_label_cell_' monpol_vars_cell{ii_monpol} ';']);
        eval(['factors2_label_cell=factors2_label_cell_' monpol_vars_cell{ii_monpol} ';']);
        no_factors          =   length(factors1_cell);
        no_factors1         =   length(factors1_cell);
        no_factors2         =   length(factors2_cell);
    
        for ii_smpl = [1 2]             
            switch ii_smpl
                case 1
                    ii_smpl_min     =   1;
                    ii_smpl_max     =   1;
                case 2
                    ii_smpl_min     =   1;
                    ii_smpl_max     =   2;
            end;                    
            for ii_factors=1:no_factors
                for ii_spreads=1:no_spreads+1;
                    tic;
                    VAR.switch_extern=switch_extern;
                    if ii_spreads<=no_spreads
                        VAR.switch_extern   =    0;    %Switch off extended VAR for simple VARs
                    end;
                    smpl_min_VAR = smpl_min_VAR_vec(ii_smpl_min,:);

                    %factors sample starts minimum p periods after the VAR sample
                    if smpl_min_factors_vec(ii_factors,1)>=smpl_min_VAR(1,1)+VAR.p/12
                        smpl_min_FACTORS = smpl_min_factors_vec(ii_factors,:);
                    else
                        smpl_min_FACTORS(1,1) = smpl_min_VAR(1,1)+floor(VAR.p/12);
                        smpl_min_FACTORS(1,2) = smpl_min_VAR(1,2)+(VAR.p-12*floor(VAR.p/12));
                    end;
                    smpl_max_VAR = smpl_max_VAR_vec(ii_smpl_max,:);
                    if (smpl_max_factors_vec(ii_factors,1)>smpl_max_VAR(1,1))
                        smpl_max_FACTORS    =   smpl_max_VAR;
                    else
                        smpl_max_FACTORS = smpl_max_factors_vec(ii_factors,:); %smpl_max_VAR; %[2012 6]; %this has crisis in it!
                    end

                    %print what is being calculated
                    fprintf('\n\n#%3.0f\n',ii{ii_monpol});
                    if ii_spreads==no_spreads+1
                        fprintf(['MONPOL: ' monpol_vars_cell{1,ii_monpol} '\nSPREADS:' figure_name ...
                            '\nFACTORS: ' factors1_label_cell{1,ii_factors} ', ' factors2_label_cell{1,ii_factors} ...
                            ', ' num2str(smpl_min_FACTORS(1,1)) '-' ...
                            '\nSAMPLE: ' num2str(smpl_min_VAR(1,1)) '-' num2str(smpl_max_VAR(1,1)) '\n']);                    
                    elseif ii_spreads>=1
                        fprintf(['MONPOL: ' monpol_vars_cell{1,ii_monpol} '\nSPREADS: ' spreads_cell{1,no_spread_single(1,ii_spreads)} ...
                            '\nFACTORS: ' factors1_label_cell{1,ii_factors} ', ' factors2_label_cell{1,ii_factors} ...
                            ' ' num2str(smpl_min_FACTORS(1,1)) '-' ...
                            '\nSAMPLE: ' num2str(smpl_min_VAR(1,1)) '-' num2str(smpl_max_VAR(1,1)) '\n']);
                    else
                        fprintf(['MONPOL: ' monpol_vars_cell{1,ii_monpol} ...
                            '\nFACTORS: ' factors1_label_cell{1,ii_factors} ', ' factors2_label_cell{1,ii_factors} ...
                            ', ' num2str(smpl_min_FACTORS(1,1)) '-' ...
                            '\nSAMPLE: ' num2str(smpl_min_VAR(1,1)) '-' num2str(smpl_max_VAR(1,1)) '\n']);
                    end;

                    %Find the dates in the sample
                    VAR.smpl_min_VAR = find(and(DATASET_VAR.TSERIES(:,cell2mat(values(DATASET_VAR.MAP,{'YEAR'})))==smpl_min_VAR(1,1), ...
                        DATASET_VAR.TSERIES(:,cell2mat(values(DATASET_VAR.MAP,{'MONTH'})))==smpl_min_VAR(1,2)));
                    VAR.smpl_max_VAR = find(and(DATASET_VAR.TSERIES(:,cell2mat(values(DATASET_VAR.MAP,{'YEAR'})))==smpl_max_VAR(1,1), ...
                        DATASET_VAR.TSERIES(:,cell2mat(values(DATASET_VAR.MAP,{'MONTH'})))==smpl_max_VAR(1,2)));
                    VAR.smpl_max_VAR_factors = find(and(DATASET_VAR.TSERIES(:,cell2mat(values(DATASET_VAR.MAP,{'YEAR'})))==smpl_max_FACTORS(1,1), ...
                        DATASET_VAR.TSERIES(:,cell2mat(values(DATASET_VAR.MAP,{'MONTH'})))==smpl_max_FACTORS(1,2)));        %Maximum place of factors in the VAR dataset
                    
                    VAR.smpl_min_FACTORS = find(and(DATASET_FACTORS.TSERIES(:,cell2mat(values(DATASET_FACTORS.MAP,{'YEAR'})))==smpl_min_FACTORS(1,1), ...
                            DATASET_FACTORS.TSERIES(:,cell2mat(values(DATASET_FACTORS.MAP,{'MONTH'})))==smpl_min_FACTORS(1,2)));
                    VAR.smpl_max_FACTORS = find(and(DATASET_FACTORS.TSERIES(:,cell2mat(values(DATASET_FACTORS.MAP,{'YEAR'})))==smpl_max_FACTORS(1,1), ...
                            DATASET_FACTORS.TSERIES(:,cell2mat(values(DATASET_FACTORS.MAP,{'MONTH'})))==smpl_max_FACTORS(1,2)));

                    %Select the variables in the VAR
                    VAR.select_vars      = {monpol_vars_cell{1,ii_monpol}};
                    VAR.select_vars_label= {monpol_vars_label_cell{1,ii_monpol}};
                    VAR.select_vars_label_short= {monpol_vars_label_short_cell{1,ii_monpol}};
                    ii_vars             = 1;  

                    VAR.select_vars      =  [VAR.select_vars,{'LCPI','LIP'}];   %Add prices and industrial production
                    VAR.select_vars_label=  [VAR.select_vars_label,{'CPI','IP'}];                 
                    VAR.select_vars_label_short=  [VAR.select_vars_label_short,{'CPI','IP'}];
                    
                    ii_vars     =   ii_vars+2;                    
                    if ii_spreads==no_spreads+1
                        for jj_spreads=1:no_spreads_all         %Adding all spreads
                            VAR.select_vars     =  [VAR.select_vars,{spreads_cell{1,jj_spreads}}];
                            VAR.select_vars_label=  [VAR.select_vars_label,{spreads_label_cell{1,jj_spreads}}];
                            VAR.select_vars_label_short=  [VAR.select_vars_label_short,{spreads_label_short_cell{1,jj_spreads}}];                            
                        end;
                        VAR.chol_order       = [2 3 1 3+cumsum(ones(1,no_spreads_all))];        %Cholesky ordering of the selected variables                        
                    else                                        %Add a single spread 
                        VAR.select_vars     =  [VAR.select_vars,{spreads_cell{1,no_spread_single(1,ii_spreads)}}];
                        VAR.select_vars_label=  [VAR.select_vars_label,{spreads_label_cell{1,no_spread_single(1,ii_spreads)}}];
                        VAR.select_vars_label_short=  [VAR.select_vars_label_short,{spreads_label_short_cell{1,no_spread_single(1,ii_spreads)}}];                        
                        VAR.chol_order       = [2 3 1 4];       %Cholesky ordering of the selected variables
                    end;                    
                    
                    if strcmp(factors2_cell{1,1},'')
                        VAR.select_factors   = {factors1_cell{1,ii_factors}};
                        VAR.select_factors_label   = {factors1_label_cell{1,ii_factors}};
                    else
                        VAR.select_factors   = {factors1_cell{1,ii_factors}};
                        VAR.select_factors_label   = {factors1_label_cell{1,ii_factors}};
                        for jj_factors2=1:no_factors2
                            VAR.select_factors   = [VAR.select_factors,{factors2_cell{1,jj_factors2}}];
                            VAR.select_factors_label   = [VAR.select_factors_label,{factors2_label_cell{1,ii_factors}}];
                        end;
                        VAR.select_factors1     =   VAR.select_factors;
                        VAR.select_factors2     =   VAR.select_factors;
                    end;
                    
                    
                    
                    VAR.extern_vars     =   {};
                    VAR.extern_vars_label=  {};
                    VAR.extern_vars_label_short=  {};
                    zz=1;
                    for uu=1:no_extern_vars
                        if ~strcmp(extern_vars_cell{1,uu},monpol_vars_cell{ii_monpol}) %drop if it is the monetary policy variable
                                VAR.extern_vars     =   [VAR.extern_vars,{extern_vars_cell{1,uu}}];
                                VAR.extern_vars_label=  [VAR.extern_vars_label,{extern_vars_label_cell{1,uu}}];
                                VAR.extern_vars_label_short=  [VAR.extern_vars_label_short,{extern_vars_label_short_cell{1,uu}}];
                                VAR.extern_vars_smpl_min_vec(zz,:)=extern_vars_smpl_min_vec(uu,:);
                                if VAR.extern_vars_smpl_min_vec(zz,1)<smpl_min_VAR(1,1)
                                    VAR.extern_vars_smpl_min_vec(zz,:)=smpl_min_VAR(1,:);
                                end;
                                VAR.smpl_min_VAR_e(1,zz) = find(and(DATASET_VAR.TSERIES(:,cell2mat(values(DATASET_VAR.MAP,{'YEAR'})))==VAR.extern_vars_smpl_min_vec(zz,1), ...
                                    DATASET_VAR.TSERIES(:,cell2mat(values(DATASET_VAR.MAP,{'MONTH'})))==VAR.extern_vars_smpl_min_vec(zz,2)));
                                zz=zz+1;
                        end;
                    end;
                    
                    %Create a vector of variables that are plotted as term
                    %spreads
                    term_spreads_vec = NaN(1,no_term_spreads);
                    for uu=1:no_term_spreads
                        term_spreads_temp=find(strcmp(term_spreads_cell{1,uu},VAR.select_vars),1);
                        if ~isempty(term_spreads_temp)
                            term_spreads_vec(1,uu)=term_spreads_temp;
                        end;
                    end;
                    for uu=1:no_term_spreads
                        term_spreads_temp=find(strcmp(term_spreads_cell{1,uu},VAR.extern_vars),1);
                        if ~isempty(term_spreads_temp)
                            term_spreads_vec(1,uu)=term_spreads_temp+length(VAR.select_vars);   %external variables after main variables
                        end;
                    end;

                    VAR.term_spreads    =   term_spreads_vec;
                    VAR.term_spreads_init = term_spreads_init;
                    VAR.term_spreads_matur = term_spreads_matur;
                    VAR.term_spreads_label_cell = term_spreads_label_cell;
                    VAR.term_spreads_label_short_cell = term_spreads_label_short_cell;
                    
                    %Create a vector of variables that are plotted as
                    %excess premia
                    excess_return_vec = NaN(1,no_excess_return);
                    for uu=1:no_excess_return
                        excess_return_temp=find(strcmp(excess_return_cell{1,uu},VAR.select_vars),1);
                        if ~isempty(excess_return_temp)
                            excess_return_vec(1,uu)=excess_return_temp;
                        end;
                    end;
                    for uu=1:no_excess_return
                        excess_return_temp=find(strcmp(excess_return_cell{1,uu},VAR.extern_vars),1);
                        if ~isempty(excess_return_temp)
                            excess_return_vec(1,uu)=excess_return_temp+length(VAR.select_vars);   %external variables after main variables
                        end;
                    end;

                    excess_return_matur_vec = NaN(1,no_excess_return);
                    for uu=1:no_excess_return
                        excess_return_matur_temp=find(strcmp(excess_return_matur_cell{1,uu},term_spreads_cell),1);
                        if ~isempty(excess_return_matur_temp)
                            excess_return_matur_vec(1,uu)=excess_return_matur_temp;   %external variables after main variables
                        end;
                    end;
                    
                    VAR.excess_return    =   excess_return_vec;
                    VAR.excess_return_matur = excess_return_matur_vec;
                    VAR.excess_return_label_cell = excess_return_label_cell;
                    VAR.excess_return_label_short_cell = excess_return_label_short_cell;
                    
                    %Create a vector of variables that are plotted as real
                    %rates
                    real_rates_vec = NaN(1,no_real_rates);
                    for uu=1:no_real_rates
                        real_rates_temp=find(strcmp(real_rates_cell{1,uu},VAR.select_vars),1);
                        if ~isempty(real_rates_temp)
                            real_rates_vec(1,uu)=real_rates_temp;
                        end;
                    end;
                    for uu=1:no_real_rates
                        real_rates_temp=find(strcmp(real_rates_cell{1,uu},VAR.extern_vars),1);
                        if ~isempty(real_rates_temp)
                            real_rates_vec(1,uu)=real_rates_temp+length(VAR.select_vars);   %external variables after main variables
                        end;
                    end;

                    VAR.real_rates    =   real_rates_vec;
                    VAR.real_rates_init = real_rates_init;                    
                    VAR.real_rates_matur = real_rates_matur;
                    VAR.real_rates_label_cell = real_rates_label_cell;
                    VAR.bkeven_label_cell = bkeven_label_cell;
                    VAR.real_rates_label_short_cell = real_rates_label_short_cell;
                    VAR.bkeven_label_short_cell = bkeven_label_short_cell;

                    if VAR.switch_exp==1
                    %Create a vector of expectations 
                    exp_rates_vec = NaN(1,no_exp_rates);
                    for uu=1:no_exp_rates
                        exp_rates_temp=find(strcmp(exp_rates_cell{1,uu},VAR.select_vars),1);
                        if ~isempty(exp_rates_temp)
                            exp_rates_vec(1,uu)=exp_rates_temp;
                        end;
                    end;
                    for uu=1:no_exp_rates
                        exp_rates_temp=find(strcmp(exp_rates_cell{1,uu},VAR.extern_vars),1);
                        if ~isempty(exp_rates_temp)
                            exp_rates_vec(1,uu)=exp_rates_temp+length(VAR.select_vars);   %external variables after main variables
                        end;
                    end;

                    exp_rates_comp_vec = NaN(1,no_exp_rates);
                    for uu=1:no_exp_rates
                        exp_rates_comp_temp=find(strcmp(exp_rates_comp_cell{1,uu},VAR.select_vars),1);
                        if ~isempty(exp_rates_comp_temp)
                            exp_rates_comp_vec(1,uu)=exp_rates_comp_temp;
                        end;
                    end;
                    for uu=1:no_exp_rates
                        exp_rates_comp_temp=find(strcmp(exp_rates_comp_cell{1,uu},VAR.extern_vars),1);
                        if ~isempty(exp_rates_comp_temp)
                            exp_rates_comp_vec(1,uu)=exp_rates_comp_temp+length(VAR.select_vars);   %external variables after main variables
                        end;
                    end;                    

                    %Find the index of e_matur in term_spreads_cell
                    exp_rates_matur_ind_vec = NaN(1,no_exp_rates);
                    for uu=1:no_exp_rates
                        exp_rates_matur_ind_vec(1,uu)=find(strcmp(exp_rates_comp_cell{1,uu},term_spreads_cell),1);
                    end;
                                        
                    VAR.exp_rates       =   exp_rates_vec;
                    VAR.exp_rates_comp  =   exp_rates_comp_vec;                    
                    VAR.exp_rates_matur = exp_rates_matur;
                    VAR.exp_rates_matur_ind = exp_rates_matur_ind_vec;                    
                    VAR.exp_rates_label_cell = exp_rates_label_cell;
                    VAR.exp_rates_label_short_cell = exp_rates_label_short_cell;                    
                    end;

                    %Run the VAR
                    [VAR,VARChol,VARbs,VARCholbs]=doVAR(VAR,DATASET_VAR,DATASET_FACTORS,nboot,clevel);
                     
                    %Create figures
                    if ii_spreads==no_spreads+1
                            nCol_e    =   2;
                            nRow_e    =   3;
                            nRow    =   3;
                        if length(VAR.select_vars)>6
                            nCol    =   3;
                        else
                            nCol    =   2;
                        end;
                        no_fig=plot_figure_nochol(VAR,VARChol,VARbs,VARCholbs,nRow,nCol,nRow_e,nCol_e,ii{ii_monpol},VAR.switch_extern);
                        ii{ii_monpol} = ii{ii_monpol}+no_fig;
                    else
                        nRow_e    =   3;
                        nCol_e    =   2;                                                    
                        if length(VAR.select_vars)>4
                                nRow    =   3;
                                nCol    =   2;                            
                                no_fig=plot_figure(VAR,VARChol,VARbs,VARCholbs,nRow,nCol,ii{ii_monpol},VAR.switch_extern);
                        elseif length(VAR.select_vars)==4
                                nRow    =   4;
                                nCol    =   2;                            
                                no_fig=plot_figure_sep(VAR,VARChol,VARbs,VARCholbs,nRow,nCol,ii{ii_monpol},VAR.switch_extern);
                        else
                                nRow    =   3;
                                nCol    =   2;                            
                                no_fig=plot_figure_sep(VAR,VARChol,VARbs,VARCholbs,nRow,nCol,ii{ii_monpol},VAR.switch_extern);
                        end;
                        ii{ii_monpol} = ii{ii_monpol}+no_fig;
                    end;
                end;
            end;
        end;
end;
