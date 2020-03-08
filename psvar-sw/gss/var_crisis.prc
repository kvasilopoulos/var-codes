/* Var Procs for Crisis Paper */
proc(8) = var_crisis(s,lag,tfirst,tlast,iconst,icomp,dfactor,tnow);

/* Computes VAR and covariance matrix of estimated parameters
   
   Input:
   		    s = txn matrix of series
   		    lag = number of lags
   		    tfirst = index of first obs to use
   		    tlast = index of last obs to use
   		    iconst = 0 Do Not Include Constant
   		             1 Include Constant
   		    icomp = 0 Do not compute Companion Matrices
   		            1 Compute companion matrix (ignoring contsant term if iconst = 1)
   		            2 Compute Companion matrix (adding constant as last element of state if iconst = 1)
   		    dfactor = discount factor for discouted least squares
   		    tnow = time period of interest for discounted least squares (irrelevant if dfactor = 1);
   		            
   Output:
          beta = estimated VAR coefficients (including constant)
                 each column has coefficients for a single equation (at time tnow)
          seps = innovation covariance matrix (at time tnow)
          eps =  VAR residuals
          Companion Form Parameter, computed in icomp = 1;
          Q, M, G for model written in SS form
               y(t) = Q*z(t)
               z(t) = M*z(t-1) + G*u(t)
               var(u(t)) = I
         icomp = 2 ... constant term is added to companion state vector
*/

local ns,dtrend,wght,yv,x,i,tmp,beta,eps,seps,m,q,g,b,comp,nobs,ndf,trnd,trndpack,n1,n2, const_coef;

@ Set Up VAR @

if tfirst .< lag+1; "VAREST: Increasing tfirst to tfirst=lag+1"; tfirst=lag+1; endif;
if tlast .> rows(s); "VAREST: Decreasing tfirst to tlast=Last Obs of Data Matrix"; tlast=rows(s); endif;
ns=cols(s);

@ Construct vector of weights using dfactor @
dtrend = seqm(1,dfactor,rows(s));
wght = miss(zeros(rows(s),1),0);
wght[1:tnow] = rev(dtrend[1:tnow]);
if tnow .< rows(s);
 wght[tnow+1:rows(s)] = dtrend[2:rows(s)-tnow+1];
endif;

wght = wght[tfirst:tlast];
yv=s[tfirst:tlast,.];
x=ones(rows(yv),1);
i=1; do while i<=lag;
  x=x~s[tfirst-i:tlast-i,.];
i=i+1; endo;

if iconst .== 0;     @ Eliminate Constant @
 x = x[.,2:cols(x)];
endif;

@ Eliminate Missing Values @
tmp=yv~x~wght;
trnd=seqa(1,1,rows(tmp));
tmp=trnd~tmp;
tmp=packr(tmp);
trndpack=tmp[.,1];
tmp=tmp[.,2:cols(tmp)];
n1=trndpack[1]-1;                             @ Number of periods lost in beginning @
n2=trnd[rows(trnd)]-trndpack[rows(trndpack)]; @ Number of periods lost at end @


format /rd 3,0;
if n1 ./= 0; 
 "VAREST: Missing Values at beginning of period, losing ";;n1;; " data points"; 
 tfirst = tfirst + n1;
endif;

if n2 ./= 0; 
 "VAREST: Missing Values at End of period, losing ";;n2;; " data points"; 
 tlast = tlast - n2;
endif;

yv=tmp[.,1:cols(yv)];
x=tmp[.,cols(yv)+1:cols(yv)+cols(x)];
wght = tmp[.,cols(tmp)];

@ Set Average Wght to equal unity @
wght = wght/meanc(wght);
wght = sqrt(wght);
yv = yv.*wght;
x = x.*wght;
/*format /rd 5,4;
wght;
*/
@ Estimate VAR @
beta=invpd(x'x)*x'yv;
eps=yv-x*beta;
nobs=rows(x);
ndf=rows(x)-cols(x);
seps=(eps'eps)/ndf;


m=0; q=0; g=0;
@
   Transform the VAR so that it is written in Standard form as:
   s(t)=P1*s(t-1) + P2*s(t-2) + ... + Pvarlag*s(t-varlag) + eps(t)
@
if icomp .> 0;
if iconst .== 0;
 b = beta';
 const_coef = zeros(cols(s),1);
else;
 b=beta[2:rows(beta),.]';
 const_coef = beta[1,.]';      @ Coefficients on Constant Term @
endif; 
 comp=zeros(cols(b),cols(b));
 comp[1:rows(b),.]=b;
 if cols(b) .> rows(b);
  comp[rows(b)+1:rows(comp),1:cols(comp)-rows(b)]=eye(cols(comp)-rows(b));
 endif;
 @ -- Write Model in SS Form --
     y(t) = Q*z(t)
     z(t) = M*z(t-1) + G*u(t)
     var(u(t)) = I
 @
 m=comp;
 q=zeros(ns,rows(m));
 q[1:ns,1:ns]=eye(ns);
 g=zeros(cols(m),ns);
 g[1:ns,1:ns]=(chol(seps))';
endif;

if icomp .== 2;  @ Add Constant Term @
 g=g|zeros(1,cols(g));
 q=q~zeros(rows(q),1);
 m=m~zeros(rows(m),1);
 m=m|zeros(1,cols(m));
 m[rows(m),cols(m)]=1;
 m[1:ns,cols(m)]=const_coef;
endif;

@ Pad eps with missing values @
tmp = miss(zeros(rows(s),cols(s)),0);
tmp[tfirst:tlast,.] = eps;
eps = tmp;

retp(beta, seps, eps, nobs, ndf, q, m, g);

endp;
