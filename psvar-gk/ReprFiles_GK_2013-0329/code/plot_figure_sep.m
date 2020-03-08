function no_fig=plot_figure_sep(VAR,VARChol,VARbs,VARCholbs,nCol,nRow,M,switch_extern)

fig=figure(M);
set(gcf,'DefaultAxesFontSize',VAR.fontsize);
set(gcf,'DefaultTextFontSize',VAR.fontsize);
M_new=M;
xmin    =   1;
xmax    =   VAR.irhor;
ymin    =   NaN(1,VAR.n);
ymax    =   NaN(1,VAR.n);
for nvar=1:VAR.n
    %Create axes
    ymin_temp =  min(min([VARCholbs.irsL(:,nvar); VARbs.irsL(:,VARChol.chol_order(nvar))]));
    ymax_temp =  max(max([VARCholbs.irsH(:,nvar); VARbs.irsH(:,VARChol.chol_order(nvar))]));
    ydiff = abs(ymax_temp-ymin_temp);
    ymin(1,VARChol.chol_order(nvar)) =  floor(10*(ymin_temp-ydiff*0.025))/10;
    ymax(1,VARChol.chol_order(nvar)) =  ceil(10*(ymax_temp+ydiff*0.025))/10;
end;

for nvar=1:VAR.n
    
    subplot(nCol,nRow,2*nvar-1); plot(1:VAR.irhor,VAR.irs(:,nvar),'LineWidth',2,'Color','k'); title(VAR.select_vars_label(1,nvar));  axis manual; axis([xmin xmax ymin(1,nvar) ymax(1,nvar)]); ylabel('%');
    hold on
    subplot(nCol,nRow,2*nvar-1); plot(1:VAR.irhor,VARbs.irsH(:,nvar),'LineWidth',1,'Color','k','LineStyle','--');
    hold on
    subplot(nCol,nRow,2*nvar-1); plot(VARbs.irsL(:,nvar),'LineWidth',1,'Color','k','LineStyle','--');
    set(gca,'YGrid','off','XGrid','on');
    grid
    if nvar==1
        legend('External Instruments');
    end;
    
end;
for nvar=1:VAR.n
    subplot(nCol,nRow,2*VARChol.chol_order(nvar)); plot(1:VARChol.irhor,VARChol.irs(:,nvar),'LineWidth',2,'Color','r');  title(VAR.select_vars_label(1,VARChol.chol_order(nvar))); axis manual; axis([xmin xmax ymin(1,VARChol.chol_order(nvar)) ymax(1,VARChol.chol_order(nvar))]);  set(gcf, 'Color', 'w'); ylabel('%');
    hold on
    subplot(nCol,nRow,2*VARChol.chol_order(nvar)); plot(1:VARChol.irhor,VARCholbs.irsH(:,nvar),'LineWidth',1,'Color','r','LineStyle','--');
    hold on
    subplot(nCol,nRow,2*VARChol.chol_order(nvar)); plot(VARCholbs.irsL(:,nvar),'LineWidth',1,'Color','r','LineStyle','--');
    hold on

    set(gca,'YGrid','off','XGrid','on');
    grid
    if VARChol.chol_order(nvar)==1
        legend('Cholesky');
    end;
end;

if switch_extern==1
yy=VAR.n+1;
for nvar=1:VAR.n_e
    if yy>(nCol*nRow);     %Starting new figure if necessary 
        [M_new,fig,yy]=new_figure(VAR,M_new,fig);
    end;
    subplot(nCol,nRow,yy); p3=plot(1:VAR.irhor,VAR.irs_e(:,nvar),'LineWidth',2,'Color','k'); title(VAR.extern_vars_label(1,nvar)); set(gcf, 'Color', 'w'); ylabel('%');
    hold on
    subplot(nCol,nRow,yy); plot(1:VAR.irhor,VARbs.irsH_e(:,nvar),'LineWidth',1,'Color','k','LineStyle','--');
    hold on
    subplot(nCol,nRow,yy); plot(VARbs.irsL_e(:,nvar),'LineWidth',1,'Color','k','LineStyle','--');
    hold on

    set(gca,'YGrid','off','XGrid','on');
    grid

    if yy==1
        legend([p3],'External Instruments');
    end;
    yy=yy+1;
end;
end;

[M_new,fig,yy]=new_figure(VAR,M_new,fig);

close(fig);
M_new=M_new-1;

no_fig=M_new-M+1;