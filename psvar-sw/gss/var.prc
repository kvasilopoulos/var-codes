/* ... Miscellaneous VAR Procedures ... */

proc(9) = varest(s,lag,tfirst,tlast,ivarmat,icomp,i_const);

/* Computes VAR and covariance matrix of estimated parameters
   
   Input:
   		    s = txn matrix of series
   		    lag = number of lags
   		    tfirst = index of first obs to use
   		    tlast = index of last obs to use
   		    ivarmat = 0 Do not compute Covariance Matrix of VAR parameters
   		              1 Compute covariance matrix using classical formula (with block diagonality between beta and sigma)
   		              2 Compute covariance matrix using GMM formula (hetero-robust and allowing covariance between beta and sigma)
   		    icomp = 0 Do not compute Companion Matrices
   		            1 Compute companion matrix ignoring contsant term
   		            2 Compute Companion matrix adding constant as last element of state -- loadings are zero if iconst = 0
   		    i_const = 0 Do not include constant
   		              1 Include constant
   		            
   Output:
          beta = estimated VAR coefficients (including constant if i_const == 1)
                 each column has coefficients for a single equation
          seps = innovation covariance matrix
          eps =  VAR residuals
          nobs = number of observations used for estimation
          ndf =  number of degrees of freedom in each equation 
          v_bs = covariance matrix of vec(beta)|vech(seps) ... 
                 computed if ivarmat = 1
          Companion Form Parameter, computed in icomp = 1;
          Q, M, G for model written in SS form
               y(t) = Q*z(t)
               z(t) = M*z(t-1) + G*u(t)
               var(u(t)) = I
         icomp = 2 ... constant term is added to companion state vector (if i_const == 1)
*/

local ns, yv, x, i, tmp, beta, eps, seps, v_bs, v_beta, z, t, v_seps, cv_bs, m, q, g, b, comp, nobs, ndf, vhseps;
local trnd, trndpack, n1, n2, xxi_x, j, w, const_coef;

@ Set Up VAR @

if tfirst .< lag+1; "VAREST: Increasing tfirst to tfirst=lag+1"; tfirst=lag+1; endif;
if tlast .> rows(s); "VAREST: Decreasing tfirst to tlast=Last Obs of Data Matrix"; tlast=rows(s); endif;
ns=cols(s);

yv=s[tfirst:tlast,.];
x=ones(rows(yv),1);
i=1; do while i<=lag;
  x=x~s[tfirst-i:tlast-i,.];
i=i+1; endo;
if i_const .== 0;
 x=x[.,2:cols(x)];  @ Eliminate Constant if i_const == 0 @
endif;

@ Eliminate Missing Values @
tmp=yv~x;
trnd=seqa(1,1,rows(tmp));
tmp=trnd~tmp;
tmp=packr(tmp);
trndpack=tmp[.,1];
tmp=tmp[.,2:cols(tmp)];
n1=trndpack[1]-1;  @ Number of periods lost in beginning @
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

@ Estimate VAR @
beta=invpd(x'x)*x'yv;
eps=yv-x*beta;
nobs=rows(x);
ndf=rows(x)-cols(x);
seps=(eps'eps)/ndf;


v_bs=0;
if ivarmat .== 1;
@ Calculate Covariance Matrix of Estimated VAR Parameters @
   @ --- Covariance Matrix of Vec(beta) --- @
     v_beta=seps.*.(invpd(x'x));

   @ --- Covariance Matrix of Vech(seps) --- @
     vhseps=vech(seps);
     z=zeros(rows(vhseps),rows(eps));
     t=1; do while t <= rows(eps);
      z[.,t] = vech(eps[t,.]'eps[t,.])-vhseps;
     t=t+1; endo;
     v_seps=z*z'/(cols(z)*cols(z));

   @ --- Covariance Between Vec(beta) and Vech(seps) --- @
    cv_bs = zeros(rows(v_beta),rows(v_seps));
     
  @ -- Put Results Together, Cov matrix of [Vec(beta)|Vech(seps)] --- @
   v_bs = (v_beta~cv_bs)|(cv_bs'~v_seps);
endif;

if ivarmat .== 2;
 @ Step 1 -- Compute vec(beta_hat - beta) -- @
  xxi_x=invpd(x'x)*x';
  w=miss(zeros(rows(beta)*cols(beta),rows(eps)),0);
  i=1; do while i <= cols(beta);
 	 j=(i-1)*rows(beta);
 	 w[j+1:j+rows(beta),.]=xxi_x.*eps[.,i]';
  i=i+1; endo;
  @ Step 2 -- Compute vech(sigma_hat-sigma) @
     vhseps=vech(seps);
     z=zeros(rows(vhseps),rows(eps));
     t=1; do while t <= rows(eps);
      z[.,t] = vech(eps[t,.]'eps[t,.])-vhseps;
     t=t+1; endo;
     z=z/cols(z);         @ Divide by T @
  @ Compute Covariance Matrix @
  w=w|z;
  v_bs=w*w';
endif;

m=0; q=0; g=0;
if icomp .> 0;
 @
   Transform the VAR so that it is written in Standard form as:
   s(t)=P1*s(t-1) + P2*s(t-2) + ... + Pvarlag*s(t-varlag) + eps(t)
 @

 @ ---- Calculate Companion Matrix ---- @
 if i_const .== 0;
  b = beta';
  const_coef = zeros(cols(s),1);
 elseif i_const .== 1;
  b=beta[2:rows(beta),.]';
  const_coef = beta[1,.]';      @ Coefficients on Constant Term @
 else;
  "Invalid value of i_const -- stopping";stop;
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

retp(beta, seps, eps, nobs, ndf, v_bs, q, m, g);

endp;

proc(3)=var_to_comp(beta,seps,icomp);

@ Convert VAR Parameters into Companion Form @
/*   Input:   Beta   In format from Varest
              Seps:  In format from Varest
              icomp = 2 ... constant term added to state vector

     Output: Q, M, G for model written in SS form
               y(t) = Q*z(t)
               z(t) = M*z(t-1) + G*u(t)
               var(u(t)) = I
*/
local ns,b,comp,m,q,g;

ns=rows(seps);

 @
   Transform the VAR so that it is written in Standard form as:
   s(t)=P1*s(t-1) + P2*s(t-2) + ... + Pvarlag*s(t-varlag) + eps(t)
 @

 @ ---- Calculate Companion Matrix ---- @
 b=beta[2:rows(beta),.]';
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

if icomp .== 2;  @ Add Constant Term @
 g=g|zeros(1,cols(g));
 q=q~zeros(rows(q),1);
 m=m~zeros(rows(m),1);
 m=m|zeros(1,cols(m));
 m[rows(m),cols(m)]=1;
 m[1:ns,cols(m)]=beta[1,.]';
endif;

retp(q, m, g);
endp;


/* ... Miscellaneous VAR Procedures ... */

proc(8) = varest_i1(s,lag,tfirst,tlast,icomp,i1_ind);

/* Computes VAR and covariance matrix of estimated parameters
   
   Input:
   		    s = txn matrix of series
   		    lag = number of lags
   		    tfirst = index of first obs to use
   		    tlast = index of last obs to use
   		    icomp = 0 Do not compute Companion Matrices
   		            1 Compute companion matrix ignoring contsant term
   		            2 Compute Companion matrix adding constant as last element of state
   		    i1_ind = nx1 vector, 0 = I(0), 1 = I(1)
   		             note: levels of data are in S, for series with I_ind=1, the I(1)
   		                   constraint is imposed on the VAR (using restricted least squares formulae)


   Output:
          beta = estimated VAR coefficients (including constant)
                 each column has coefficients for a single equation
          seps = innovation covariance matrix
          eps =  VAR residuals
          nobs = number of observations used for estimation
          ndf =  number of degrees of freedom in each equation 
          Companion Form Parameter, computed in icomp = 1;
          Q, M, G for model written in SS form
               y(t) = Q*z(t)
               z(t) = M*z(t-1) + G*u(t)
               var(u(t)) = I
*/

local ns, yv, x, i, tmp, eps, seps, v_bs, v_beta, z, t, v_seps, cv_bs, m, q, g, b, comp, nobs, ndf, vhseps;
local trnd, trndpack, n1, n2;
local beta_ols,xxi,j,rsum,r,beta;

@ Set Up VAR @

if tfirst .< lag+1; "VAREST: Increasing tfirst to tfirst=lag+1"; tfirst=lag+1; endif;
if tlast .> rows(s); "VAREST: Decreasing tfirst to tlast=Last Obs of Data Matrix"; tlast=rows(s); endif;
ns=cols(s);

yv=s[tfirst:tlast,.];
x=ones(rows(yv),1);
i=1; do while i<=lag;
  x=x~s[tfirst-i:tlast-i,.];
i=i+1; endo;

@ Eliminate Missing Values @
tmp=yv~x;
trnd=seqa(1,1,rows(tmp));
tmp=trnd~tmp;
tmp=packr(tmp);
trndpack=tmp[.,1];
tmp=tmp[.,2:cols(tmp)];
n1=trndpack[1]-1;  @ Number of periods lost in beginning @
n2=trnd[rows(trnd)]-trndpack[rows(trndpack)]; @ Number of periods lost at end @

format /rd 3,0;
if n1 ./= 0; 
 "VAREST: Missing Values at beginning of period, losing ";;n1;; " data points"; 
endif;

if n2 ./= 0; 
 "VAREST: Missing Values at End of period, losing ";;n2;; " data points"; 
endif;

yv=tmp[.,1:cols(yv)];
x=tmp[.,cols(yv)+1:cols(yv)+cols(x)];

@ Estimate VAR -- Unrestricted @
xxi=invpd(x'x);
beta_ols=xxi*(x'yv);
@ Set up Matrix that sums coefficients @
Rsum=zeros(cols(s),rows(beta_ols));
i=1; do while i <= lag;
	j=cols(s)*(i-1)+1;    @ 1 for constant term @
	rsum[.,j+1:j+cols(s)]=eye(cols(s));
i=i+1; endo;
r=eye(cols(s));

@ Save only those rows of Rsum and r for integrated variables @
rsum=selif(rsum,i1_ind);
r=selif(r,i1_ind);

if sumc(i1_ind) .== 0;
 beta=beta_ols;
else;
 beta=beta_ols-xxi*rsum'(invpd(rsum*xxi*rsum'))*(rsum*beta_ols-r);
endif;

eps=yv-x*beta;
nobs=rows(x);
ndf=rows(x)-cols(x);
seps=(eps'eps)/ndf;

m=0; q=0; g=0;
if icomp .> 0;
 @
   Transform the VAR so that it is written in Standard form as:
   s(t)=P1*s(t-1) + P2*s(t-2) + ... + Pvarlag*s(t-varlag) + eps(t)
 @

 @ ---- Calculate Companion Matrix ---- @
 b=beta[2:rows(beta),.]';
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
 m[1:ns,cols(m)]=beta[1,.]';
endif;


@ Pad eps with missing values @
tmp = miss(zeros(rows(s),cols(s)),0);
tmp[tfirst:tlast,.] = eps;
eps = tmp;

retp(beta, seps, eps, nobs, ndf, q, m, g);

endp;

/* ----------------------------------------- */

proc(2) = varspec(s,lag,nord,tfirst,tlast);

/* Computes Spectrum from VAR
   
   Input:
   		    s = txn matrix of series
   		    lag = number of lags
   		    nord = number of ordinates
   		    tfirst = index of first obs to use
   		    tlast = index of last obs to use

   Output:
          specmat = nord x (n*n) spectrum (normalized by 2pi)
          wvec = vector of ordinates          
*/

local step, wvec, ns, specmat, ivarmat, icomp, beta, seps, nobs, ndf, v_bs, q, m, g, i, w, ss, eps;

step=pi/nord;
wvec=seqa(step,step,nord); 
ns=cols(s);
specmat=miss(zeros(nord,ns*ns),0);

 ivarmat=0;
 icomp=1;
 {beta, seps, eps, nobs, ndf, v_bs, q, m, g} =varest(s,lag,tfirst,tlast,ivarmat,icomp);
 i=1; do while i <= nord;
  w=wvec[i];
  ss=spmodh(q,m,g,w);
  specmat[i,.]=(vec(ss))';
 i=i+1; endo;
 specmat=specmat/(2*pi);
 retp(specmat,wvec);
endp;

/* ----------------------------------------------------------------- */

proc(1) = spmodh(h,f,g,w);
/* Procedure for Calculating Spectral Density Matrix for a VAR of the
   form:
     x(t)=h*s(t)
     s(t)=f*s(t-1) + g*e(t)

     where var(e(t))=I

   Output:  ss = spectrum of x at frequency w.
   Note: Spectrum is not divided by 2*pi

*/
local im, z, sm, smi, ss;


 let im = 0+1i;
 z=exp(-w*im);
 sm=eye(rows(f));
 sm=sm-z*f;
 smi=inv(sm);
 ss=h*smi*g*g'*smi'*h';

retp(ss);
endp;

/* ----------------------------------------------------------------- */