function [M_new,fig,yy]=new_figure(VAR,M_new,fig)

fig_pos=get(fig,'Position');
set(fig,'Position',[100 100 900 850]);        

h=axes('Position',[0,0,1,1],'Xlim',[0 1],'Ylim',[0 1]);
set(h,'Visible','off');
string_F_m=['First stage regression: F: ' num2str(VAR.F_m(1),'%2.2f')];
text('Position',[0.1 0.05],'string',string_F_m);
string_F_m_rob=['robust F: ' num2str(VAR.F_m_rob(1),'%2.2f')];
text('Position',[0.41 0.05],'string',string_F_m_rob);
string_R2_m = ['R2: ' num2str(VAR.R2_m(1)*100,'%1.2f') '%'];
text('Position',[0.575 0.05],'string',string_R2_m);
string_R2adj_m = ['Adjusted R2: ' num2str(VAR.R2adj_m(1)*100,'%1.2f') '%'];
text('Position',[0.7 0.05],'string',string_R2adj_m);        

M_new=M_new+1;
fig=figure(M_new); set(gcf,'DefaultAxesFontSize',VAR.fontsize); set(gcf,'DefaultTextFontSize',VAR.fontsize);  set(gcf, 'Color', 'w');
yy=1;

