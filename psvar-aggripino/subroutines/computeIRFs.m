
function res=computeIRFs(dataStructure,modelSpec)
%computes IRFs to shock of interest under different identification schemes
%
%inputs:
%dataStructure  =dataset structure
%modelSpec      =structure with model specification
%
% output: structure (responses to the chosen structural shock only)
% irfs            :[(nH+1)xn] matrix of responses to structural shock
% irfs_u & irfs_l :[(nH+1)xn] matrices of error bands
% relevance       : scalar (nan if Cholesky)
%
% miranda 2015 smirandaagrippino@london.edu

%--------------------------------------------------------------------------

%unpack input structures
select=modelSpec.dataSelection;

%load relevant data series
y =dataStructure.data(:,select);  


%model size
n =size(y,2); 
modelSpec.modelSize =n;


%model basics
nL     =modelSpec.nLags; 
nH     =modelSpec.nHorizons; 
shockS =modelSpec.shockSize;


%identification scheme
iScheme=modelSpec.identification;


%estimate VAR
VAR =estimateVAR(y,nL);


%identification
switch iScheme
    case 'CHOL'
        
        Bzero =chol(cov(VAR.resid)); 
        Bzero =bsxfun(@rdivide,Bzero,diag(Bzero));
        
    case 'PSVAR'
        
        instrumentSet           =loadInstrument(dataStructure,modelSpec);
        modelSpec.instrumentSet =instrumentSet;
        
        Bzero =psvar(VAR.resid,modelSpec,instrumentSet)';

end

%structural IRFs
irfs =SirfBuild(VAR.beta,Bzero,n,nL,nH,shockS);

%bands (wild bootstrap)
bands =computeSirfBands(y,modelSpec);


%load final structure
res.irfs   =irfs;
res.irfs_u =bands.upper;
res.irfs_l =bands.lower;





%--child functions--------------------------------------------------------%

function res=estimateVAR(data,nLags)
%estimate vector autoregression with constant

[T,n] =size(data); 
nL    =nLags;


%build matrix of relevant lagged endogenous
Ylag=nan(T-nL,n*nL); %[y_{t-1},...,y_{t-p}]';

for j=1:nL
    
    Ylag(:,n*(j-1)+1:n*j)=data(nL-j+1:end-j,:);
end

nT =size(Ylag,1); 
Y  =data(nL+1:end,:);

%VAR coefficients
beta =[ones(nT,1) Ylag]\Y; 

%VAR innovations
innovations=Y-[ones(nT,1) Ylag]*beta;


%load final structure
res.Y       =Y; 
res.Ylag    =Ylag;
res.beta    =beta;
res.resid   =innovations;
res.data    =data;
res.nlags   =nLags;
%-------------------------------------------------------------------------%



function irf=SirfBuild(B,Bzero,n,nL,nH,shockS)
%compute structural impulse responses Bzero: u_{t}=Bzero e_{t}
%expects B=[n*nL+1 x n]

%companion form
A=zeros(n*nL,n*nL); 

A(1:n,:)=B(2:end,:)'; 
A(n+1:end,1:n*(nL-1))=eye(n*(nL-1));

%irfs
irf=nan(n,nH+1); irf(:,1)=eye(n)*Bzero'*shockS;
for h=1:nH
    
    Ah=A^h; irf(:,h+1)=Ah(1:n,1:n)*Bzero'*shockS;
    
end
%-------------------------------------------------------------------------%


function res=computeSirfBands(y,modelSpec)
%calculates bands around IRFs using wild bootstrap


%basics
n  =modelSpec.modelSize;
nL =modelSpec.nLags;
nH =modelSpec.nHorizons;
nC =modelSpec.cLevel;
nB =modelSpec.bootSize;


%identification and shock size
iScheme =modelSpec.identification; 
shockS  =modelSpec.shockSize;


%estimate VAR
VAR =estimateVAR(y,nL);

%bootstrap impulse response functions
bootIRFs =nan(n,nH+1,nB);

for b=1:nB
    
    bootSet =computeBootstrapSamples(VAR,modelSpec);
    
    Yboot   =bootSet.Y;
    
    VARboot =estimateVAR(Yboot,nL);
    
    switch iScheme
        
        case 'CHOL'

            Bzero =chol(cov(VARboot.resid)); 
            Bzero =bsxfun(@rdivide,Bzero,diag(Bzero));

        case 'PSVAR'
            
            bootProxy       =modelSpec.instrumentSet; 
            bootProxy.data  =bootSet.Z;
            
            Bzero =psvar(VARboot.resid,modelSpec,bootProxy)';

    end
    
    %
    bootIRFs(:,:,b) =SirfBuild(VARboot.beta,Bzero,n,nL,nH,shockS);
    
    if mod(b,100)==0
        
        fprintf(1,'bootstrap irf: %i of %i \n',b,nB);
    end
        
end

%get relevant quantiles
bootIRFs =sort(bootIRFs,3);

uBound =nC+(100-nC)/2;   uBound =uBound/100;
lBound =(100-nC)/2;      lBound =lBound/100;


%load final structure
res.upper =bootIRFs(:,:,round(uBound*nB));
res.lower =bootIRFs(:,:,round(lBound*nB));
%-------------------------------------------------------------------------%


function res=computeBootstrapSamples(VAR,modelSpec)
%returns wild bootstrap set of endogenous in main VAR (Y); endogenous in
%VAR with fed funds used for long-run expectations (X); instrument (Z)


Yinnovations =VAR.resid;    
betaY        =VAR.beta;

%
[nT,nY] =size(Yinnovations); 
nL      =VAR.nlags;

%bootstrap algorithm: randomly switch sign of the reduced form innovations
bSwitch =1-2*(rand(nT,1)>.5);

%bootstrapped residuals
bootUy  =Yinnovations.*repmat(bSwitch,1,nY);

%build bootsrap sample
bootY   =VAR.data;       bootY(nL+1:end,:)=nan; 
for t=1:nT
    
    Ylag =flipud(bootY(t:t+nL-1,:))'; 
    Ylag =[1; Ylag(:)];
    
    bootY(t+nL,:) =Ylag'*betaY + bootUy(t,:);

end

%identification
iScheme =modelSpec.identification;
switch iScheme
    
    case 'CHOL'
        
        %load final structure
        res.Y =bootY;
        
    case 'PSVAR'
        
        instrumentSet =modelSpec.instrumentSet;
        
        %bootstrapped instruments
        bootZ =instrumentSet.data;       nZ =size(bootZ,2);
        bootZ =bootZ.*repmat(bSwitch(instrumentSet.commonSample),1,nZ);

        %load final structure
        res.Y =bootY;
        res.Z =bootZ;
end
%-------------------------------------------------------------------------%


function Bzero =psvar(VARinnovations,modelSpec,instrumentSet)
%identification using external instruments

n =size(VARinnovations,2);


%instrument
proxy       =instrumentSet.data;

%trim VAR innovations to match the instruments sample
innovations =VARinnovations(instrumentSet.commonSample,:)';

keepRows    =~isnan(proxy);
proxy       =proxy(keepRows,:);
innovations =innovations(:,keepRows);

%identification
res =ProxySVARidentification(innovations,find(modelSpec.shockVar),proxy);
res

%Bzero
Bzero=zeros(n,n); Bzero(:,modelSpec.shockVar)=res.B/res.B(modelSpec.shockVar);
%--------------------------------------------------------------------------


function res =loadInstrument(dataStructure,modelSpec)

%load IV variables
load instruments
    
externalInstrument =IV;

%find relevant proxy variable
selectN  =ismember(lower(externalInstrument.labels),lower(modelSpec.selectedInstrument));
selectT  =~isnan(IV.data(:,selectN));


IVdata  =IV.data(selectT,selectN);
IVdates =IV.dates(selectT);
    

%define common time bounds to align data to instruments
lowerT =max(dataStructure.dates(1+modelSpec.nLags),IVdates(1));
upperT =min(dataStructure.dates(end),IVdates(end));

commonTimeLine =lowerT;
while commonTimeLine(end)<upperT
    
    commonTimeLine=[commonTimeLine; addtodate(commonTimeLine(end),1,'month')];
end

commonSample =ismember(dataStructure.dates,commonTimeLine);


%trim non overlapping instrument observations & update structure
selectTime                      =ismember(IVdates,commonTimeLine);
externalInstrument.data         =IVdata(selectTime,:);
externalInstrument.dates        =IVdates(selectTime,:);
externalInstrument.commonSample =commonSample(modelSpec.nLags+1:end);


%load in final structure
res =externalInstrument;

