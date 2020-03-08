function no_fig=plot_figure_nochol(VAR,VARChol,VARbs,VARCholbs,nCol,nRow,nCol_e,nRow_e,M,switch_extern)

xmin    =   1;
xmax    =   VAR.irhor;

if VAR.switch_extern==1
    string_cell     =   {'','_e','_e_matur','_ts','_rr','_er','_bkeven'};
    if VAR.switch_exp==1
        string_cell = [string_cell, {'_exp'}];
    end;
    for gg=1:VAR.n_e
        string_cell = [string_cell, {['_e_all{1,' num2str(gg) '}']}];
    end;                
else
    string_cell     =   {''};
end;
no_string_cell  =   length(string_cell);
for ii=1:no_string_cell
    eval(['irs' string_cell{1,ii} '    =   VAR.irs'  string_cell{1,ii} '(:,:);']);
    eval(['irsH' string_cell{1,ii} '    =   VARbs.irsH'  string_cell{1,ii} '(:,:);']);
    eval(['irsL' string_cell{1,ii} '    =   VARbs.irsL'  string_cell{1,ii} '(:,:);']);
end;

fig=figure(M);
set(gcf,'DefaultAxesFontSize',VAR.fontsize);
set(gcf,'DefaultTextFontSize',VAR.fontsize);
for nvar=1:VAR.n
    subplot(nCol,nRow,nvar); p1=plot(1:VAR.irhor,irs(:,nvar),'LineWidth',2,'Color','k'); title(VAR.select_vars_label(1,nvar)); set(gcf, 'Color', 'w');  ylabel('%');
    hold on
    subplot(nCol,nRow,nvar); plot(1:VAR.irhor,irsH(:,nvar),'LineWidth',1,'Color','k','LineStyle','--');
    hold on
    subplot(nCol,nRow,nvar); plot(irsL(:,nvar),'LineWidth',1,'Color','k','LineStyle','--');
    hold on
    
    set(gca,'YGrid','off','XGrid','on');
    grid
end;

M_new = M;
if switch_extern==1
[M_new,fig,yy]=new_figure(VAR,M_new,fig);

%The same axes 
ymin_temp    =   100;
ymax_temp    =   -100;
for nvar=1:nCol_e*nRow_e
    %Create axes
    ymin_temp =  min(min(irsL_e(:,nvar)),ymin_temp);
    ymax_temp =  max(max(irsH_e(:,nvar)),ymax_temp);
end;
ydiff = abs(ymax_temp-ymin_temp);
ymin_e =  floor(10*(ymin_temp-ydiff*0.025))/10;
ymax_e =  ceil(10*(ymax_temp+ydiff*0.025))/10;

for nvar=1:VAR.n_e
    subplot(nCol_e,nRow_e,yy); plot(1:VAR.irhor,irs_e(:,nvar),'LineWidth',2,'Color','k'); title(VAR.extern_vars_label(1,nvar)); set(gcf, 'Color', 'w'); axis([xmin xmax ymin_e ymax_e]); ylabel('%');
    hold on
    subplot(nCol_e,nRow_e,yy); plot(1:VAR.irhor,irsH_e(:,nvar),'LineWidth',1,'Color','k','LineStyle','--');
    hold on
    subplot(nCol_e,nRow_e,yy); plot(irsL_e(:,nvar),'LineWidth',1,'Color','k','LineStyle','--');
    hold on
        
    set(gca,'YGrid','off','XGrid','on');
    grid

    yy=yy+1;
    if and(yy>(nCol_e*nRow_e),nvar<VAR.n_e)    %Starting a new figure if necessary
        [M_new,fig,yy]=new_figure(VAR,M_new,fig);
        %The same axes 
        ymin_temp    =   100;
        ymax_temp    =   -100;
        for zz=1:min(nCol_e*nRow_e,VAR.n_e-nvar);
            %Create axes
                ymin_temp =  min(min(irsL_e(:,nvar+zz)),ymin_temp);
                ymax_temp =  max(max(irsH_e(:,nvar+zz)),ymax_temp);
        end;
        ydiff = abs(ymax_temp-ymin_temp);
        ymin_e =  floor(10*(ymin_temp-ydiff*0.025))/10;
        ymax_e =  ceil(10*(ymax_temp+ydiff*0.025))/10;
    end;
end;

if yy>1    %Starting a new figure if necessary
        [M_new,fig,yy]=new_figure(VAR,M_new,fig);
end;


%The same axes 
ymin_temp    =   100;
ymax_temp    =   -100;
for nvar=1:VAR.n_ts_1
    %Create axes
    ymin_temp =  min(min(irsL_ts(:,nvar)),ymin_temp);
    ymax_temp =  max(max(irsH_ts(:,nvar)),ymax_temp);
end;
for nvar=1:VAR.n_er
    %Create axes
    ymin_temp =  min(min(irsL_er(:,nvar)),ymin_temp);
    ymax_temp =  max(max(irsH_er(:,nvar)),ymax_temp);
end;
ydiff = abs(ymax_temp-ymin_temp);
ymin_ts_1 =  floor(10*(ymin_temp-ydiff*0.025))/10;
ymax_ts_1 =  ceil(10*(ymax_temp+ydiff*0.025))/10;


irs_all     = [irs irs_e];
irs_term    = irs_all(:,VAR.term_spreads);


for nvar=1:VAR.n_ts_1
        
    subplot(nCol_e,nRow_e,yy); p3=plot(1:VAR.irhor,irs_ts(:,nvar),'LineWidth',2,'Color','k'); title(VAR.term_spreads_label_cell(1,nvar)); set(gcf, 'Color', 'w'); axis([xmin xmax ymin_ts_1 ymax_ts_1]); ylabel('%');
    hold on
    subplot(nCol_e,nRow_e,yy); plot(1:VAR.irhor,irsH_ts(:,nvar),'LineWidth',1,'Color','k','LineStyle','--');
    hold on
    subplot(nCol_e,nRow_e,yy); plot(irsL_ts(:,nvar),'LineWidth',1,'Color','k','LineStyle','--');
    hold on

    subplot(nCol_e,nRow_e,yy); p4=plot(1:VAR.irhor,irs_term(:,nvar),'LineWidth',2,'Color','r');
    
     if yy==1
         legend([p3 p4],'Term premium','Nominal rates');
     end;

    set(gca,'YGrid','off','XGrid','on');
    grid
    yy=yy+1;
    if and(yy>(nCol_e*nRow_e),nvar<(VAR.n_ts_1+VAR.n_er))    %Starting a new figure if necessary
        [M_new,fig,yy]=new_figure(VAR,M_new,fig);
    end;    
end;    

irs_all      =   [irs irs_e];
irs_term     = irs_all(:,VAR.term_spreads);
irs_er_tot   = irs_all(:,VAR.excess_return)+irs_term(:,VAR.excess_return_matur);

for nvar=1:VAR.n_er
    subplot(nCol_e,nRow_e,yy); p5=plot(1:VAR.irhor,irs_er(:,nvar),'LineWidth',2,'Color','k'); title(VAR.excess_return_label_cell(1,nvar)); set(gcf, 'Color', 'w'); axis([xmin xmax ymin_ts_1 ymax_ts_1]); ylabel('%');
    hold on
    subplot(nCol_e,nRow_e,yy); plot(1:VAR.irhor,irsH_er(:,nvar),'LineWidth',1,'Color','k','LineStyle','--');
    hold on
    subplot(nCol_e,nRow_e,yy); plot(irsL_er(:,nvar),'LineWidth',1,'Color','k','LineStyle','--');
    hold on

    subplot(nCol_e,nRow_e,yy); p6=plot(1:VAR.irhor,irs_er_tot(:,nvar),'LineWidth',2,'Color','r');
    
     if nvar==1
         legend([p5 p6],'Combined excess premium','Nominal rates');
     end;
    
    
    set(gca,'YGrid','off','XGrid','on');
    grid
    yy=yy+1;
    if and(yy>(nCol_e*nRow_e),nvar<VAR.n_er)    %Starting a new figure if necessary
        [M_new,fig,yy]=new_figure(VAR,M_new,fig);
    end;    
end;    

if yy>1    %Starting a new figure if necessary
        [M_new,fig,yy]=new_figure(VAR,M_new,fig);
end;


%The same axes 
ymin_temp    =   100;
ymax_temp    =   -100;
for nvar=VAR.n_ts_1+1:VAR.n_ts %max(nCol_e*nRow_e,VAR.n_ts-VAR.n_ts_1)
    %Create axes
    ymin_temp =  min(min(irsL_ts(:,nvar)),ymin_temp);
    ymax_temp =  max(max(irsH_ts(:,nvar)),ymax_temp);
end;
ydiff = abs(ymax_temp-ymin_temp);
ymin_ts_2 =  floor(10*(ymin_temp-ydiff*0.025))/10;
ymax_ts_2 =  ceil(10*(ymax_temp+ydiff*0.025))/10;


for nvar=VAR.n_ts_1+1:VAR.n_ts
        
    subplot(nCol_e,nRow_e,yy); p3=plot(1:VAR.irhor,irs_ts(:,nvar),'LineWidth',2,'Color','k'); title(VAR.term_spreads_label_cell(1,nvar)); set(gcf, 'Color', 'w'); axis([xmin xmax ymin_ts_2 ymax_ts_2]); ylabel('%');
    hold on
    subplot(nCol_e,nRow_e,yy); plot(1:VAR.irhor,irsH_ts(:,nvar),'LineWidth',1,'Color','k','LineStyle','--');
    hold on
    subplot(nCol_e,nRow_e,yy); plot(irsL_ts(:,nvar),'LineWidth',1,'Color','k','LineStyle','--');
    hold on

    subplot(nCol_e,nRow_e,yy); p4=plot(1:VAR.irhor,irs_term(:,nvar),'LineWidth',2,'Color','r');
    
     if yy==1
         legend([p3 p4],'Term premium','Nominal rates');
     end;

    set(gca,'YGrid','off','XGrid','on');
    grid
    yy=yy+1;
    if and(yy>(nCol_e*nRow_e),nvar<VAR.n_ts)    %Starting a new figure if necessary
        [M_new,fig,yy]=new_figure(VAR,M_new,fig);
    end;    
end;    

if yy>1    %Starting a new figure if necessary
        [M_new,fig,yy]=new_figure(VAR,M_new,fig);
end;

%The same axes 
ymin_temp    =   100;
ymax_temp    =   -100;
for nvar=1:VAR.n_rr %max(nCol_e*nRow_e,VAR.n_ts-VAR.n_ts_1)
    %Create axes
    ymin_temp =  min(min(irsL_rr(:,nvar)),ymin_temp);
    ymax_temp =  max(max(irsH_rr(:,nvar)),ymax_temp);
end;
ydiff = abs(ymax_temp-ymin_temp);
ymin_rr =  floor(10*(ymin_temp-ydiff*0.025))/10;
ymax_rr =  ceil(10*(ymax_temp+ydiff*0.025))/10;


irs_all     =   [irs irs_e];
irs_nom     =   irs_all(:,VAR.real_rates);

for nvar=1:VAR.n_rr
        
    subplot(nCol_e,nRow_e,2*nvar-1); p1=plot(1:VAR.irhor,irs_rr(:,nvar),'LineWidth',2,'Color','k'); title(VAR.real_rates_label_cell(1,nvar)); set(gcf, 'Color', 'w'); axis([xmin xmax ymin_rr ymax_rr]); ylabel('%');
    hold on
    subplot(nCol_e,nRow_e,2*nvar-1); plot(1:VAR.irhor,irsH_rr(:,nvar),'LineWidth',1,'Color','k','LineStyle','--');
    hold on
    subplot(nCol_e,nRow_e,2*nvar-1); plot(irsL_rr(:,nvar),'LineWidth',1,'Color','k','LineStyle','--');
    hold on

    subplot(nCol_e,nRow_e,2*nvar-1); p2=plot(1:VAR.irhor,irs_nom(:,nvar),'LineWidth',2,'Color','r');
    
     if (2*nvar-1)==1
         legend([p1 p2],'Real','Nominal');
     end;
    
    set(gca,'YGrid','off','XGrid','on');
    grid
    yy=yy+1;
    if and(yy>(nCol_e*nRow_e),nvar<VAR.n_rr)    %Starting a new figure if necessary
        [M_new,fig,yy]=new_figure(VAR,M_new,fig);
    end;    
end;    

%The same axes 
ymin_temp    =   100;
ymax_temp    =   -100;
for nvar=1:VAR.n_rr %max(nCol_e*nRow_e,VAR.n_ts-VAR.n_ts_1)
    %Create axes
    ymin_temp =  min(min(irsL_bkeven(:,nvar)),ymin_temp);
    ymax_temp =  max(max(irsH_bkeven(:,nvar)),ymax_temp);
end;
ydiff = abs(ymax_temp-ymin_temp);
ymin_bkeven =  floor(10*(ymin_temp-ydiff*0.025))/10;
ymax_bkeven =  ceil(10*(ymax_temp+ydiff*0.025))/10;


for nvar=1:VAR.n_rr
        
    subplot(nCol_e,nRow_e,2*nvar); plot(1:VAR.irhor,irs_bkeven(:,nvar),'LineWidth',2,'Color','k'); title(VAR.bkeven_label_cell(1,nvar)); set(gcf, 'Color', 'w'); axis([xmin xmax ymin_bkeven ymax_bkeven]); ylabel('%');
    hold on
    subplot(nCol_e,nRow_e,2*nvar); plot(1:VAR.irhor,irsH_bkeven(:,nvar),'LineWidth',1,'Color','k','LineStyle','--');
    hold on
    subplot(nCol_e,nRow_e,2*nvar); plot(irsL_bkeven(:,nvar),'LineWidth',1,'Color','k','LineStyle','--');
    hold on

    set(gca,'YGrid','off','XGrid','on');
    grid
    yy=yy+1;
    if and(yy>(nCol_e*nRow_e),nvar<VAR.n_rr)    %Starting a new figure if necessary
        [M_new,fig,yy]=new_figure(VAR,M_new,fig);
    end;    
end;    

if VAR.switch_exp==1
%plotting expectations
if yy>1    %Starting a new figure if necessary
        [M_new,fig,yy]=new_figure(VAR,M_new,fig);
end;

irs_all     =   [irs irs_e];
irs_comp    =   irs_all(:,VAR.exp_rates_comp);
irs_exp_matur = irs_e_matur(:,VAR.exp_rates_matur_ind);
nCol_exp    =   2;
nRow_exp    =   2;  %4

%The same axes 
ymin_temp    =   100;
ymax_temp    =   -100;
for nvar=1:VAR.n_exp
    %Create axes
    ymin_temp =  min(min(irsL_exp(:,nvar)),ymin_temp);
    ymax_temp =  max(max(irsH_exp(:,nvar)),ymax_temp);
end;
ydiff = abs(ymax_temp-ymin_temp);
ymin_exp =  floor(10*(ymin_temp-ydiff*0.025))/10;
ymax_exp =  ceil(10*(ymax_temp+ydiff*0.025))/10;


for nvar=1:VAR.n_exp
%    yy = (nvar<=4)*(2*nvar-1)+(nvar>4)*(2*(nvar-4));
    subplot(nRow_exp,nCol_exp,yy); p1=plot(1:VAR.irhor,irs_exp(:,nvar),'LineWidth',2,'Color','k'); title(VAR.exp_rates_label_cell(1,nvar)); set(gcf, 'Color', 'w');  axis([xmin xmax ymin_exp ymax_exp]); ylabel('%');
    hold on
    subplot(nRow_exp,nCol_exp,yy); plot(1:VAR.irhor,irsH_exp(:,nvar),'LineWidth',1,'Color','k','LineStyle','--');
    hold on
    subplot(nRow_exp,nCol_exp,yy); plot(irsL_exp(:,nvar),'LineWidth',1,'Color','k','LineStyle','--');
    hold on

    subplot(nRow_exp,nCol_exp,yy); p2=plot(1:VAR.irhor,irs_comp(:,nvar),'LineWidth',2,'Color','r','LineStyle','-');  %-.

    subplot(nRow_exp,nCol_exp,yy); p3=plot(1:VAR.irhor,irs_exp_matur(:,nvar),'LineWidth',2,'Color','b','LineStyle',':');
    
     if nvar==1
         legend([p1 p2 p3],'Expectations','Nominal','Implied riskless rates');
     end;
    
    set(gca,'YGrid','off','XGrid','on');
    grid
    yy=yy+1;
    if and(yy>(nCol_exp*nRow_exp),nvar<VAR.n_exp)    %Starting a new figure if necessary
        [M_new,fig,yy]=new_figure(VAR,M_new,fig);
    end;    
end;    

end;

for gg=1:VAR.n_e
    if yy>1    %Starting a new figure if necessary
        [M_new,fig,yy]=new_figure(VAR,M_new,fig);
    end;
    nCol_g=nCol;
    nRow_g=nRow;
    for nvar=1:VAR.n
            subplot(nCol_g,nRow_g,nvar); p1=plot(1:VAR.irhor,irs_e_all{1,gg}(:,nvar),'LineWidth',2,'Color','k'); title(VAR.select_vars_label(1,nvar)); set(gcf, 'Color', 'w');  ylabel('%');
            hold on
            subplot(nCol_g,nRow_g,nvar); plot(1:VAR.irhor,irsH_e_all{1,gg}(:,nvar),'LineWidth',1,'Color','k','LineStyle','--');
            hold on
            subplot(nCol_g,nRow_g,nvar); plot(irsL_e_all{1,gg}(:,nvar),'LineWidth',1,'Color','k','LineStyle','--');
            hold on
    
            set(gca,'YGrid','off','XGrid','on');
            grid
            if yy==1
                legend(['XI: ' VAR.extern_vars_label{1,gg} ]);
            end;

            yy=yy+1;
            if and(yy>(nCol_g*nRow_g),nvar<VAR.n_e)    %Starting a new figure if necessary
                [M_new,fig,yy]=new_figure(VAR,M_new,fig);
            end;
    end;
end;



end;

[M_new,fig,yy]=new_figure(VAR,M_new,fig);

close(fig);
M_new=M_new-1;

no_fig=M_new-M+1;