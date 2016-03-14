classdef TemperatureAnalyze3
    %     This type of class analyzes our frameset variable along the
    %     time/frames
    
    properties(GetAccess = 'public', SetAccess = 'private')
        CoordinatesToolTip;
        SlopeRakeFace;
        SlopeClearanceFace;
        HeatCarriedByChip;
        HeatCarriedByTool;
        ExceedingPoints;
        AmountValidFrames;
        CurrentInternalEnergyTool;
        MaximumTemperatureCuttingZone;
        InternalEnergyToolRate;
    end
    
    properties(GetAccess = 'private', SetAccess = 'private')
        BinaryImageTool;
        Frames;
    end
    
    methods
        
        function obj = TemperatureAnalyze3(Frames)
            obj.Frames = Frames;
            b = 20;%Beginning
            e = 38;%End
            delta = e-b+1;
            obj.AmountValidFrames = [b e delta];
            auxCoordinates = zeros(delta,2);
            auxSlopeRF = zeros(delta,1);
            auxSlopeCF = zeros(delta,1);
            auxExceedingpoints = zeros(delta,1);
            auxHeatbyTool = zeros(delta,1);
            auxHeatbyChip = zeros(delta,1);
            auxInternalEnergyTool = zeros(delta,1);
            auxBinary = struct('b',[]);
            auxMaxTempCZ = zeros(delta,1);
            for i = b:e
                h = TemperatureAnalyze2(Frames(i).f);
                auxCoordinates(i-b+1,:) = h.CoordinateToolTip;
                auxSlopeCF(i-b+1) = h.ClearanceAngle;
                auxSlopeRF(i-b+1) = h.RakeAngle;
                auxExceedingpoints(i-b+1) = h.exceedingPoints(320);
                auxHeatbyTool(i-b+1) = h.HeatFluxAwayFromToolTip;
                auxHeatbyChip(i-b+1) = h.HeatCarriedAwayByChip;
                auxInternalEnergyTool(i-b+1) = h.InternalEnergyTool;
                auxBinary(i-b+1).b = h.passBinaryImageTool();
                auxMaxTempCZ(i-b+1) = h.MaximumTemperatureCuttingZone;
            end
            obj.CoordinatesToolTip = auxCoordinates;
            obj.SlopeRakeFace = auxSlopeRF;
            obj.SlopeClearanceFace = auxSlopeCF;
            obj.HeatCarriedByChip = auxHeatbyChip;
            obj.HeatCarriedByTool = auxHeatbyTool;
            obj.ExceedingPoints = auxExceedingpoints;
            obj.CurrentInternalEnergyTool = auxInternalEnergyTool;
            obj.BinaryImageTool = auxBinary;
            obj.MaximumTemperatureCuttingZone = auxMaxTempCZ;
        end
        
        function obj = toolFormatAnalysis(obj)
            trf = obj.SlopeRakeFace;
            tcf = obj.SlopeClearanceFace;
            cTT = obj.CoordinatesToolTip;
            trf(trf == 6) = [];
            tcf(tcf == 3) = [];
            m = size(cTT,1);
            for i = m:-1:1
                if isequal(cTT(i,:),[198 73])
                    cTT(i,:) = [];
                end
            end
            figure
            boxplot(cTT)
            title('Coordinates x and y of the Tool Tip')
            figure
            boxplot(trf)
            title('Rake face slope boxplot')
            figure
            boxplot(tcf)
            title('Clearance face slope boxplot')
        end
        
        function obj = displayExceedingPoints(obj)
            v = obj.AmountValidFrames(1):obj.AmountValidFrames(2);
            y = obj.ExceedingPoints';
            plot(v,obj.ExceedingPoints,'xb')
            hold on
            p = polyfit(v,y,1);
            yp = polyval(p,v);
            plot(v,yp,'--k')
            xlabel('Frame')
            ylabel('Number of points exceeding the defined temperature')
            legend('Data points','Linear fitting','Location','SouthEast')
        end
        
        function obj = displayHeatFlux(obj)
            green = [23 156 125]/255;
            orange = [255 102 0]/255;
            Hc = obj.HeatCarriedByChip;
            Ht = obj.HeatCarriedByTool;
            v = obj.AmountValidFrames(1):obj.AmountValidFrames(2);
            figure
            plot(v,Hc,'x','Color',green)
            ym = mean(Hc);
            H = refline(0,ym);
            H.Color = 'k';
            H.LineStyle = '-.';
            ylim([0 1.5*max(Hc)])
            legend('Heat flux','Mean heat flux','Location','SE')
            ylabel('$q''_{c}  (W/m)$','Interpreter','latex')
            xlabel('Frame number','Interpreter','latex')
            title('Heat conducted away by  the chip','Interpreter','latex')
            figure
            plot(v,Ht,'x','Color',green)
            ym = mean(Ht);
            H = refline(0,ym);
            H.Color = 'k';
            H.LineStyle = '-.';
            ylim([0 1.5*max(Ht)]) 
            legend('Heat flux','Mean heat flux','Location','SE')
            ylabel('$q''_{t}  (W/m)$','Interpreter','latex')
            xlabel('Frame number','Interpreter','latex')
            title('Heat conducted away from the tool tip','Interpreter','latex')
            figure
            plot(v,Hc+Ht,'x','Color',green)
            ym = mean(Hc+Ht);
            H = refline(0,ym);
            H.Color = 'k';
            H.LineStyle = '-.';
            ylim([0 1.5*max(Hc+Ht)])
            yt = 1475*100*10^3/(60*4.4);
            H = refline(0,yt);
            H.Color = orange;
            H.LineStyle = '-.';
            legend('Experimental data','Mean total heat flux','Total heat flux produced by the cutting process','Location','SE')
            ylabel('$q''_{total}  (W/m)$','Interpreter','latex')
            xlabel('Frame number','Interpreter','latex')
            title('Heat conducted away by the chip and from the tool tip','Interpreter','latex')
        end   
        
        function obj = displayHeatAmountTool(obj)
            green = [23 156 125]/255;
%             orange = [255 102 0]/255;
            v = obj.AmountValidFrames(1):obj.AmountValidFrames(2);
            Ha = obj.CurrentInternalEnergyTool';
            figure
            plot(v,Ha,'x','Color',green)
%             lsline
            hold on
            p = polyfit(v,Ha,1);
            yp = polyval(p,v);
            plot(v,yp,'-.k')
            xlabel('Frame')
            ylabel('U'' (J/m)')
            legend('Heat amount','Linear fitting','Location','SouthEast')
            title('Internal energy of the tool of valid pixels')
%             xlabel('Frame','Interpreter','latex')
%             ylabel('U'' $(J/m)$','Interpreter','latex')
%             legend('Heat amount','Linear fitting','Location','SouthEast')
%             title('Internal energy of the tool of valid pixels','Interpreter','latex')
        end
        
        function obj = displayTempCZbehavior(obj)
            green = [23 156 125]/255;
            T = obj.MaximumTemperatureCuttingZone;
            v = obj.AmountValidFrames(1):obj.AmountValidFrames(2);
            plot(v,T,'x','LineWidth',1.5,'Color',green)
            hold on
            p = polyfit(v,T',1);
            yp = polyval(p,v);
            plot(v,yp,'-.k')
            ylabel('Maximum temperature on the cutting zone ($^{o}C$)','Interpreter','latex')
            xlabel('Frame','Interpreter','latex')
            ylim([0.5*min(T) 1.5*max(T)])
            legend('Maximum Temperature','least - squares fitting')
            hold off
        end
        
        function obj = displayIsothermsBehavior(obj)
            for i = obj.AmountValidFrames(1):obj.AmountValidFrames(2)
                h = TemperatureAnalyze2(obj.Frames(i).f);
                h.displayIsotherms();
                waitforbuttonpress
            end
        end
        
        function obj = heatFluxInsideTool(obj)
            Q = obj.CurrentInternalEnergyTool;
            v = (obj.AmountValidFrames(1):obj.AmountValidFrames(2))';
            p = polyfit(v,Q,1);
            q = obj.HeatCarriedByTool;
            y = mean(q);
            qi = q + p(1)/0.003;
            yi = mean(qi);
%             l = length(Q);
%             auxQ = zeros(l-1,1);
%             for i = 1:l-1
%                 auxQ(i) = Q(i+1) - Q(i);
%             end
%             qa = auxQ/0.003;
            obj.InternalEnergyToolRate = p(1)/0.003;
            figure
            subplot(2,1,1)
            plot(v,qi,'rx')
            H = refline(0,y);
            H.Color = 'r';
            H.LineStyle = '-.';
            title('inside')
            subplot(2,1,2)
            plot(v,q,'x')
            H = refline(0,yi);
            H.LineStyle = '-.';
            title('outside')
        end
    end
    
end