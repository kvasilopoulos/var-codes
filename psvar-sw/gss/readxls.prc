proc(5)=readxls(xlsname,sheet,ns,nd,ndesc,ncodes);

/* Read Data in Excel File 
   Assumed file structure:
   Col. 1 = Date (string format)
   Cols 2-end = series
   
   Row 1 = Series Names (String)
   Row 2:1+Ndesc = Series Description Matrix (String)
   Row 2+Ndesc:1+ndesc+ncodes = Transformation Code Matrix
   Rows 2+Ndesc+Ncodes:end = data
   
   Input:
    xlsname = name of excel file (must be in directory of calling program )
    sheet = 1,2,3, ... (Sheet in Excel File Where data are found)
    ns = number of series 
    nd = number of time periods 
    ndesc = number of descriptors for each series (>= 0)
    ncodes = number of codes for each series (>= 0)
    
   Output:
    namevec = vector of names (string), ns x 1
    descmat = Matrix of descriptions (string), ns x ndesc
    tcodemat = Matrix of transformation codes, ns x ncodes 
    datevec = vector of dates (string)
    datamat = matrix of data (ns x nd)
*/
local tmp,endcol,range,namevec,descmat,datevec,tcodemat,datamat;
local endrow_desc, firstrow_codes, endrow_codes, firstrow_data, endrow_data;

@ Get ending Column Label @
tmp=ns+1;
endcol=fixed_to_aa(tmp);
endrow_desc=1+ndesc;
firstrow_codes=1+endrow_desc;
endrow_codes=firstrow_codes+ncodes-1;
firstrow_data=1+endrow_codes;
endrow_data=firstrow_data+nd-1;

@ Get Series Names @
range="b1:" $+ endcol $+ "1";
namevec=spreadsheetreadsa(xlsname,range,sheet);
namevec=namevec';

@ Get Series descriptors @
descmat=miss(0,0);
if ndesc .> 0;
 range="b2:" $+ endcol $+ ftocv(endrow_desc,1,0);
 descmat=spreadsheetreadsa(xlsname,range,sheet);
 descmat=descmat';
endif;

@ Get Transformation Codes @
tcodemat=miss(0,0);
if ncodes .> 0;
 range="b" $+ ftocv(firstrow_codes,1,0) $+ ":" $+ endcol $+ ftocv(endrow_codes,1,0);
 tcodemat=spreadsheetreadm(xlsname,range,sheet);
 tcodemat=tcodemat';
endif;

@ Get Date String @
range="a" $+ ftocv(firstrow_data,1,0) $+ ":a" $+ ftocv(endrow_data,1,0);
datevec=spreadsheetreadsa(xlsname,range,sheet);

@ Get Data @
range="b" $+ ftocv(firstrow_data,1,0) $+ ":" $+ endcol $+ ftocv(endrow_data,1,0);
datamat=spreadsheetreadm(xlsname,range,sheet);

retp(namevec,descmat,tcodemat,datevec,datamat);
endp;

proc(1) = fixed_to_aa(x);
 @ Convert fixed value to Excel Column Heading, A-Z, AA-AZ, BA-BZ, etc. @
 local small, strvec, strx, y, y1, y2;

 small = 1.0e-04;
 if x .> 26*26;
  "x too large, stop";
  stop;
 endif;
 let strvec[26,1]=
 A
 B
 C
 D 
 E 
 F 
 G 
 H 
 I 
 J
 K 
 L
 M
 N
 O
 P
 Q
 R
 S
 T
 U
 V
 W
 X
 Y
 Z;
 strx="";
 y=(x-1)/26;
 y1=trunc(y);
 y2=(y-y1+small)*26;
 y2=floor(y2+1);
 if y1 .== 0;
  strx="" $+ strvec[y2];
 else;
  strx="" $+ strvec[y1] $+ strvec[y2];
 endif;
retp(strx);
endp;
