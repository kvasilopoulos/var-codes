proc(2) = bw_trend_se(x,bw,nma);
 @ Compute BW Estimate of Trend and SEs @
 local xmean, z, tmp, mtmp, nobs,t,dt,bw_weight,trend,acv,i,var_tmp,se;
  
 xmean = missmat(rows(x),1);
 z = seqa(1,1,rows(x));
 tmp = packr(x~z);
 x = tmp[.,1];
 z = tmp[.,2];
 nobs = rows(x);
 mtmp = zeros(nobs,1);
 trend = seqa(1,1,nobs);
 for t (1,nobs,1);
 	dt = (trend-t)/bw_bw;
  bw_weight = (15/16)*((1-dt.^2).^2);   @ Bi-Weight @
  bw_weight = bw_weight.*(abs(dt).< 1);
  bw_weight = bw_weight/sumc(bw_weight);
  mtmp[t]=bw_weight'x;
 endfor;
 xmean[z[1]:z[rows(z)]]=mtmp;
 @ Compute SEs @
 x = x - mtmp;
 @ Compute ACVs @
 acv = missmat(nma+1,1);
 i = 1; do while i <= nma+1;
 	acv[i] = x[i:nobs]'x[1:nobs+1-i]/(nobs+1-i);
 i = i+1; endo;
 var_tmp = missmat(nobs,1);
 se = missmat(rows(xmean),1);
 for t (1,nobs,1);
  dt = (trend-t)/bw_bw;
  bw_weight = (15/16)*((1-dt.^2).^2);   @ Bi-Weight @
  bw_weight = bw_weight.*(abs(dt).< 1);
  bw_weight = bw_weight/sumc(bw_weight);
  var_tmp[t] = acv[1]*(bw_weight'bw_weight);
  i = 2; do while i <= nma+1;
  	var_tmp[t] = var_tmp[t] + 2*acv[i]*(bw_weight[i:rows(bw_weight)]'bw_weight[1:rows(bw_weight)+1-i]);
  i = i+1; endo;
 endfor;
 se[z[1]:z[rows(z)]]=sqrt(var_tmp); 
retp(xmean,se);
endp;