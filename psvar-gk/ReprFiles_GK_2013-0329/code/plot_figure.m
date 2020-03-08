function no_fig=plot_figure(VAR,VARChol,VARbs,VARCholbs,nCol,nRow,M,switch_extern)

fig=figure(M);
set(gcf,'DefaultAxesFontSize',VAR.fontsize);
set(gcf,'DefaultTextFontSize',VAR.fontsize);
M_new=M;
for nvar=1:VAR.n
    subplot(nCol,nRow,nvar); p1=plot(1:VAR.irhor,VAR.irs(:,nvar),'LineWidth',2,'Color','k'); title(VAR.select_vars_label(1,nvar)); set(gcf, 'Color', 'w');
    hold on
    subplot(nCol,nRow,nvar); plot(1:VAR.irhor,VARbs.irsH(:,nvar),'LineWidth',1,'Color','k','LineStyle','--');
    hold on
    subplot(nCol,nRow,nvar); plot(VARbs.irsL(:,nvar),'LineWidth',1,'Color','k','LineStyle','--');
    hold on
    
    subplot(nCol,nRow,VARChol.chol_order(nvar)); p2=plot(1:VARChol.irhor,VARChol.irs(:,nvar),'LineWidth',2,'Color','r'); 
    hold on
    subplot(nCol,nRow,VARChol.chol_order(nvar)); plot(1:VARChol.irhor,VARCholbs.irsH(:,nvar),'LineWidth',1,'Color','r','LineStyle','--');
    hold on
    subplot(nCol,nRow,VARChol.chol_order(nvar)); plot(VARCholbs.irsL(:,nvar),'LineWidth',1,'Color','r','LineStyle','--');
    hold on

    set(gca,'YGrid','off','XGrid','on');
    grid
    if nvar==1
        legend([p1 p2],'External Instruments','Cholesky');
    end;
end;

if switch_extern==1
yy=VAR.n+1;
for nvar=1:VAR.n_e
    if yy>(nCol*nRow);     %Starting new figure if necessary 
        [M_new,fig,yy]=new_figure(VAR,M_new,fig);
    end;
    subplot(nCol,nRow,yy); p3=plot(1:VAR.irhor,VAR.irs_e(:,nvar),'LineWidth',2,'Color','k'); title(VAR.extern_vars_label(1,nvar)); set(gcf, 'Color', 'w');
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