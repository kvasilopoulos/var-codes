display1= cell2mat(values(VAR.MAP,plotdisplay));
for nvar = 1:length(display1)
                  
        f=figure;    
        box off

            plot(0:VAR.irhor-1,VAR.irs(:,display1(nvar),shock),'-','MarkerSize',4,'LineWidth',2,'Color', [0 0 0.5]);
            hold on
            
            if isempty(VARci)==0
                if isempty(VARci.irsH)==0
                p1=plot(0:VAR.irhor-1,VARci.irsH(:,display1(nvar)),'LineWidth',1,'Color', [0 0 0.5],'LineStyle','--');
                plot(0:VAR.irhor-1,VARci.irsL(:,display1(nvar)),'LineWidth',1,'Color', [0 0 0.5],'LineStyle','--');               
                end

                if isempty(VARci.irsH2)==0
                p2=plot(0:VAR.irhor-1,VARci.irsH2(:,display1(nvar)),'LineWidth',1,'Color', [0.9 0 0],'LineStyle','-.');
                plot(0:VAR.irhor-1,VARci.irsL2(:,display1(nvar)),'LineWidth',1,'Color', [0.9 0 0],'LineStyle','-.');
                end
                
                if isempty(VARci.irsH3)==0
                p3=plot(0:VAR.irhor-1,VARci.irsH3(:,display1(nvar)),'LineWidth',1,'Color', [0 0 0],'LineStyle','-');
                plot(0:VAR.irhor-1,VARci.irsL3(:,display1(nvar)),'LineWidth',1,'Color', [0 0 0],'LineStyle','-');
                end
                
                if isempty(VARci.irsH4)==0
                p4=plot(0:VAR.irhor-1,VARci.irsH4(:,display1(nvar)),'LineWidth',1,'Color', [0 0.5 0],'LineStyle',':');
                plot(0:VAR.irhor-1,VARci.irsL4(:,display1(nvar)),'LineWidth',1,'Color', [0 0.5 0],'LineStyle',':');
                end
            end
            hline(0,'k-')
            ti=title(FIGLABELS{nvar});
            xl=xlabel('horizon (quarters)');
            if DATASET.UNIT(cell2mat(values(DATASET.MAP,{plotdisplay{nvar}})))==1
            yl=ylabel('percent');
            elseif DATASET.UNIT(cell2mat(values(DATASET.MAP,{plotdisplay{nvar}})))==2
                yl=ylabel('percentage points');
            end
            set(gca,'XTick',0:4:VAR.irhor)
            if isempty(FIG.AXIS)==0
            axis([0 VAR.irhor+1 FIG.AXIS(nvar,1) FIG.AXIS(nvar,2)])
            end
        if legendflag == 1 
                    if isempty(VARci.irsH4)==0
                        l = legend([p1,p2,p3,p4],legendlabels);
                    elseif isempty(VARci.irsH3)==0
                        l = legend([p1,p2,p3],legendlabels);
                    elseif isempty(VARci.irsH2)==0
                         l = legend([p1,p2],legendlabels);               
                    end
                     set(l,'FontName', 'AvantGarde','FontSize',10,'Location','NorthEast','EdgeColor','white','box','off')
        end


        box off
            set([xl,yl], 'FontName', 'AvantGarde','FontSize',14);
            set([ti], 'FontName', 'AvantGarde','FontSize',16);           
end