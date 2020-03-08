proc(3) = est_eta(z,eps,neta1);

/*  The structural error is give by eta_1 = b'eps. Z is a vector of instruments that is correlated with eta_1 but uncorrelated with the other eta's. 
*/
local tmp, zp, ep, alpha, beta, gam, lmbd, eta_1, tstat, df, pvalue;

format /ro 11,4;
pvalue = missmat(1,1);
if neta1 .> cols(z);
 "Too Few Instruments";
 retp(pvalue,pvalue,pvalue);
endif;
tmp=packr(z~eps);
zp = tmp[.,1:cols(z)];
ep = tmp[.,cols(z)+1:cols(z)+cols(eps)];
tmp = 1;
{alpha,beta,gam,lmbd}=canon(zp,ep,tmp,neta1);
eta_1 = eps*alpha;
if cols(z) .> 1;
 @ lmbd'; @
 tmp = ones(rows(lmbd)-neta1,1)-lmbd[neta1+1:rows(lmbd)];
 if minc(tmp) .> .0000001;
   tstat = -rows(ep)*(sumc(ln(tmp)));
   df = (cols(zp)-neta1)*(cols(ep)-neta1);
   pvalue = cdfchic(tstat,df);
   /*
   format /rd 6,2;
   "Test for RR Restriction, tstat, df, pvalue:";;tstat;;df;;
   format /rd 10,4;
   pvalue;
   */
  endif;
endif;

retp(eta_1,lmbd[1:neta1],pvalue);
endp;
 
 