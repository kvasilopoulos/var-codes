proc(1) = andrews_eos(y,x,m);

@ Andrews End of sample stability test @
/*  y = xb + u
    check for stability in last m periods
     
    output 
    pvalue = pvalue of test statistic
*/
local n, bhat, d, uhat, sigmahat, t, sinv, xj, uj, aj, vj, sj, s, j, ii, jj, bj, ujhat, pvalue;
local svec, yj;
n = rows(y)-m;
d = cols(x);
 
 @ Estimate Full Sample regression @
bhat = y/x;
uhat = y-x*bhat;

@ Estimate Covariance Matrix over m periods @
sigmahat = zeros(m,m);
for t (1,n+1,1);
	sigmahat = sigmahat + uhat[t:t+m-1]*uhat[t:t+m-1]';
endfor;
sigmahat = sigmahat/(n+1);
sinv = invpd(sigmahat);

xj = x[n+1:n+m,.];
uj = uhat[n+1:n+m];
if m .>= d;
  aj = xj'sinv*uj;
  vj = xj'sinv*xj;
  sj = aj'invpd(vj)*aj;
else;
  sj = uj'sinv*uj;
endif;

s = sj;  @ Test Statistic @

@ Compute Simulations @
svec = zeros(n-m+1,1);
for j (1,n-m+1,1);
	ii = ones(n,1);
	jj = trunc(m/2);
	ii[j:j+jj-1]=zeros(jj,1);
	yj = selif(y[1:n],ii);
	xj = selif(x[1:n,.],ii);
	bj = yj/xj;
	ujhat = y-x*bj;
	xj = x[j:j+m-1,.];
  uj = ujhat[j:j+m-1];
  if m .>= d;
   aj = xj'sinv*uj;
   vj = xj'sinv*xj;
   sj = aj'invpd(vj)*aj;
  else;
   sj = uj'sinv*uj;
  endif;
   svec[j] = sj;
 endfor;
 
 pvalue = meanc(svec .> s);
 
retp(pvalue);
endp;


