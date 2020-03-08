proc(4) = canon(y,x,z,r);

/* 
   mww, 6/7/2004
   This procedure computes the r largest canonical correlations
   between y and x after partialling out z
   
   The regression model has the form
   
   Y = X*Gam + e
   
   where Y is Nxp
         X is Nxk
         Gam is kxp
         
         Rank(Gam) = r
         
         Gam is factored as Gam=Alpha*Beta'
         where Alhpa is kxr and Beta is pxr
         
         Note, an alternative way to write this is 
         
         y(i) = Beta*alpha'x(i) + e(i)
         
         Alpha is normalized so that Alpha'Sxx*Alpha = I(r)
         
         where Sxx = X'X/N
         
    Inputs:
         Y
         X
         Z = (arbitrary scalar if not used)
         r
     Outputs
         Alpha
         Beta
         Gam
         lamda = column vector of squared canonical correlations
         
     Note: This program would receive a C- in a computing course, but it is transparent.
   
*/
local b, n, sxxi, syyi, cxxi, sxy, tmp, vax, vex, lmbd, alpha, beta, gam;

@ Step 1; partial out Z @
n=rows(x);
if rows(z) .== n;
 b=y/z;
 y=y-z*b;
 b=x/z;
 x=x-z*b;
endif;

@ Step 2; construct Moment Matrices @
sxxi=invpd(x'x/n);
syyi=invpd(y'y/n);
cxxi=chol(sxxi);
sxy=(x'y/n);
tmp=cxxi*sxy;
tmp=tmp*syyi*tmp';
{vax,vex}=eighv(tmp);
vax=rev(vax);
vex=(rev(vex'))';
lmbd=vax;
alpha=vex[.,1:r];
alpha=cxxi'alpha;
beta=sxy'alpha;
gam=alpha*beta';
retp(alpha,beta,gam,lmbd);
endp;


