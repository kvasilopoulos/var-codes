proc(3) = y4plot(y,yf,tcode);

  local y4,y4f,y4r, dnobs, tmp1;
    
  dnobs = rows(y);
  y4=miss(zeros(dnobs,1),0);
  y4f=miss(zeros(dnobs,1),0);
  y4r=miss(zeros(dnobs,1),0);  
 	
 	if tcode .== 1;
 	 y4[5:dnobs]=y[5:dnobs]-y[1:dnobs-4];
   y4f[5:dnobs]=yf[5:dnobs]-yf[1:dnobs-4];
     
  elseif tcode .== 2;
   y4[4:dnobs]=y[4:dnobs]+y[3:dnobs-1]+y[2:dnobs-2]+y[1:dnobs-3];
   y4f[4:dnobs]=yf[4:dnobs]+yf[3:dnobs-1]+yf[2:dnobs-2]+yf[1:dnobs-3];
  
  elseif tcode .== 3;
   tmp1 = miss(zeros(dnobs,1),0);
   tmp1[4:dnobs]=y[4:dnobs]+y[3:dnobs-1]+y[2:dnobs-2]+y[1:dnobs-3];
   y4[4:dnobs]=tmp1[4:dnobs]+tmp1[3:dnobs-1]+tmp1[2:dnobs-2]+tmp1[1:dnobs-3];
   
   tmp1 = miss(zeros(dnobs,1),0);
   tmp1[4:dnobs]=yf[4:dnobs]+yf[3:dnobs-1]+yf[2:dnobs-2]+yf[1:dnobs-3];
   y4f[4:dnobs]=tmp1[4:dnobs]+tmp1[3:dnobs-1]+tmp1[2:dnobs-2]+tmp1[1:dnobs-3]; 

  elseif tcode .== 4;
 	 y4[5:dnobs]=y[5:dnobs]-y[1:dnobs-4];
   y4f[5:dnobs]=yf[5:dnobs]-yf[1:dnobs-4]; 

  elseif tcode .== 5;
   y4[4:dnobs]=y[4:dnobs]+y[3:dnobs-1]+y[2:dnobs-2]+y[1:dnobs-3];
   y4f[4:dnobs]=yf[4:dnobs]+yf[3:dnobs-1]+yf[2:dnobs-2]+yf[1:dnobs-3];
      
  elseif tcode .== 6;
   tmp1 = miss(zeros(dnobs,1),0);
   tmp1[4:dnobs]=y[4:dnobs]+y[3:dnobs-1]+y[2:dnobs-2]+y[1:dnobs-3];
   y4[4:dnobs]=tmp1[4:dnobs]+tmp1[3:dnobs-1]+tmp1[2:dnobs-2]+tmp1[1:dnobs-3];
   
   tmp1 = miss(zeros(dnobs,1),0);
   tmp1[4:dnobs]=yf[4:dnobs]+yf[3:dnobs-1]+yf[2:dnobs-2]+yf[1:dnobs-3];
   y4f[4:dnobs]=tmp1[4:dnobs]+tmp1[3:dnobs-1]+tmp1[2:dnobs-2]+tmp1[1:dnobs-3]; 
 
  endif;
  
  y4r=y4-y4f;
  
 retp(y4,y4f,y4r);
endp;
     