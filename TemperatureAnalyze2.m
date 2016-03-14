classdef TemperatureAnalyze2
    % This class was built to analyze the temperature inside the tool
    % shape, the temperature gradient, the isotherms...
    
    properties(GetAccess = 'public', SetAccess = 'private')
        % define the properties of the class here, (like fields of a struct)
        CoordinateToolTip;
        TemperatureToolTip;
        RakeAngle;%Rake face slope
        ClearanceAngle;%Clearance face slope
        ShearAngle;
        FrictionAngle;
        MeanTemperatureTool;
        MaximumTemperatureChip;
        MaximumTemperatureCuttingZone;
        HeatCarriedAwayByChip;
        HeatFluxAwayFromToolTip
        ;
        TotalPowerProduced;
        InternalEnergyTool;
        CuttingForceParallelToolFace;
        CuttingForcePowerDirection;
        CuttingForceUncutChipThicknessDirection;
        CuttingForceParallelShearPlane;
        CuttingForcePerpendicularShearPlane;
        CuttingForcePerpendicularToolFace;
        CoefficientFriction;
        ShearStress;
        NormalStress;
        RatioR;
        ShearEnergyVolume;
        FrictionEnergyVolume;
        CuttingVelocity;
        UnCutChipThickness;
        ContactLength;
    end
    
    properties(GetAccess = 'private', SetAccess = 'private')
        coordRF;
        coordCF;
        BW;
        lines;
        frame;
        pointCF;%auxiliar to plot the cutting edge
        pointRF;
        Tx;%auxiliar to plot the gradients of the frame
        Ty;
        biImageTool;%Binary image of the tool shape
        biImageChip;
        biShearLine;
        xyMaxTemp;%coordinates of the point inside the chip with maximum Temperature
        lineChip;
        lineTool;
        validTemperature;
        heatCapacity;
        nExcPoints;
    end
    
    methods
        % methods, including the constructor are defined in this block
        
        function obj = TemperatureAnalyze2(Frame)%constructor
            Fp = 1500;%Newtons
            Fq = 1000;%Newtons
            width = 4.4*10^-3;%meters
            Vp = 100/60;%meters/second
            tuc = 200*10^-6;%meters
            clength = [];%for the example VP41 clength = 1.5351 mm
            if isequal(clength,[])
                clength = obj.contactLength();
            end
            obj.ContactLength = clength;
            obj.CuttingVelocity = Vp*60;%m/minute
            obj.TotalPowerProduced = Fp*Vp/width;
            obj.UnCutChipThickness = tuc;
            obj.frame = Frame;            
            obj = obj.calculateCoordinates();
            if isempty(obj.coordRF)==0&&isempty(obj.coordCF)==0
                obj = obj.coordinateToolTip();
            else %Default conditions 
                obj.CoordinateToolTip = [198 73];%Average of the TT coordinates
                if isempty(obj.coordRF)
                    obj.RakeAngle = 6;
                end
                if isempty(obj.coordCF)
                    obj.ClearanceAngle = 3;
                end
            end
            obj = obj.pointsRFandCF();
            obj = obj.TempTT();
            obj = obj.calculateGradient();
            obj = obj.maximumTemperature();
            obj = obj.meanTemperatureTool();
            obj = obj.chipContour();
            obj = obj.maxTemperatureChip();
            obj = obj.findLineChip();
            obj = obj.findLineTool();
            obj = obj.heatBalance(tuc,Vp);
            obj = obj.internalEnergyTool();
            obj = obj.shearLine();
            obj = obj.forcesValues(Fp,Fq,width,tuc);
        end
        
        function obj = maximumTemperature(obj)
            obj.MaximumTemperatureCuttingZone = max(max(obj.frame));
            [~,lin] = max(obj.frame);
            [~,col] = max(max(obj.frame));
            lin = lin(col);
            obj.xyMaxTemp = [col lin];
%             imagesc(obj.frame)
%             hold on
%             plot(col,lin,'xr')
%             hold off
        end
        
        function l = contactLength(obj)
            imagesc(obj.frame)
            imdistline
            v = input('What is the value of the contact length for this frame?');
            close all
            l = 15*10^-6*v;
        end
        
        function obj = maxTemperatureTool(obj)
           A = round(obj.CoordinateToolTip);
           xt = A(1);
           yt = A(2);
           y1 = round(yt + xt*tan(obj.RakeAngle*pi/180));
           x2 = round(xt - yt*tan((90 - obj.ClearanceAngle)*pi/180));
           c = [xt 0 0 x2];
           r = [yt y1 0 0];
           B = roipoly(obj.frame,c,r);
           [m,n] = size(B);
           A = ones(m,n);
           C = A - B;
           Frame = C.*obj.frame;
           T = max(max(Frame)); 
           K = sprintf('Maximum temperature in the tool: %3.2f ?C\n',T);
           disp(K)
        end
        
        function obj = maxTemperatureChip(obj)
%             [~,lin] = max(obj.frame);
%             [~,col] = max(max(obj.frame));
%             lin = lin(col);
%             obj.xyMaxTchip = [col lin];
            Frame = obj.biImageChip.*obj.frame;
            obj.MaximumTemperatureChip = max(max(Frame));
        end
        
        function obj = meanTemperatureTool(obj)
           A = round(obj.CoordinateToolTip);
           [m,~] = size(obj.frame);
           xt = A(1);
           yt = A(2);
           y1 = round(yt + (xt - 1)*tan(obj.RakeAngle*pi/180));
           x2 = round(xt - (m - yt)*tan(pi/2 - (90 - obj.ClearanceAngle)*pi/180));
           c = [xt 0 0 x2];
           r = [yt y1 m m];
           B = roipoly(obj.frame,c,r);
           obj.biImageTool = B;
           Frame = B.*obj.frame;
%            T = input('Minimum valid temperature: ');
           obj.validTemperature = 310;
           B = Frame > obj.validTemperature;
%            imshow(B)
           Frame = B.*Frame;
%            imagesc(Frame)
           s = sum(sum(Frame));
           n = sum(sum(B));
           meanT = s/n;
           obj.MeanTemperatureTool = meanT;
        end    
        
        function obj = chipContour(obj)
            aux = round(obj.frame);
            aux(aux < 0) = 0;
            y = label2rgb(aux,'parula');
%             image(y)
            B1 = y(:,:,2) > 180 & y(:,:,3) < 200;
            B2 = obj.biImageTool ==1 & B1 == 1;
            B = B1 - B2;
            obj.biImageChip = B;
%             imshow(B)
%             K = B.*obj.frame;
%             figure
%             imagesc(K)
%             figure
%             imagesc(obj.frame)
        end
        
        function obj = displayBinary(obj)
            imshow(obj.BW);
            hold on
            plot(obj.coordRF(:,1),obj.coordRF(:,2),'bx')
            plot(obj.coordCF(:,1),obj.coordCF(:,2),'yx')
            plot(obj.CoordinateToolTip(1),obj.CoordinateToolTip(2),'xm')
            hold off
        end
        
        function obj = TempTT(obj)
            p1 = round(obj.CoordinateToolTip + 5*[-cos(obj.RakeAngle*pi/180) sin(obj.RakeAngle*pi/180)]);
            p2 = round(obj.CoordinateToolTip + 5*[-cos((90 - obj.ClearanceAngle)*pi/180) sin((90 - obj.ClearanceAngle)*pi/180)]);
            p3 = round(obj.CoordinateToolTip + 5*[-(cos(obj.RakeAngle*pi/180)+cos((90 - obj.ClearanceAngle)*pi/180)) (sin(obj.RakeAngle*pi/180)+sin((90 - obj.ClearanceAngle)*pi/180))]);
            T1 = obj.frame(p1(2),p1(1));
            T2 = obj.frame(p2(2),p2(1));
            T3 = obj.frame(p3(2),p3(1));
            TT = obj.frame(round(obj.CoordinateToolTip(2)),round(obj.CoordinateToolTip(1)));
            T = [T1 T2 T3 TT];
            obj.TemperatureToolTip = mean(T);
%             imagesc(obj.frame)
%             hold on
%             plot([p1(1) p2(1) p3(1) obj.CoordinateToolTip(1)],[p1(2) p2(2) p3(2) obj.CoordinateToolTip(2)],'xm')
        end
        
        function obj = calculateCoordinates(obj)
            obj.BW = edge(obj.frame,'sobel');
%             imshow(obj.BW)
%             figure
%             imagesc(f.*(ones(size(f))-obj.BW))
%-------------------------Finding the clearance face---------------------------- 
            [H, THETA, RHO] = hough(obj.BW,'Theta',2:5);%Hough transformation
            P  = houghpeaks(H, 10);
            obj.lines = houghlines(obj.BW, THETA, RHO, P, 'FillGap', 15,'MinLength',10);%Here we can find the lines of cutting edge and afterwards find the coordinate of the tool tip
            l = length(obj.lines);
            obj.coordCF = [];
            for i=1:l
                Theta=obj.lines(i).theta;
                t1 = obj.lines(i).point1;
                t2 = obj.lines(i).point2;
                rho = obj.lines(i).rho;
                
%                 imagesc(f)
%                 hold on
%                 plot([t1(1) t2(1)]',[t1(2) t2(2)]','m')
%                 hold off
                
                if rho < 204 && rho > 198
                    obj.coordCF = [t1;t2];
                    obj.ClearanceAngle = Theta;
                end
            end           
%-------------------------Finding the rake face----------------------------           
            [H, THETA, RHO] = hough(obj.BW,'Theta',81:85);%Hough transformation
            P  = houghpeaks(H, 10);
            obj.lines = houghlines(obj.BW, THETA, RHO, P, 'FillGap', 15,'MinLength',10);%Here we can find the lines of cutting edge and afterwards find the coordinate of the tool tip
            l = length(obj.lines);
            obj.coordRF = [];
            for i=1:l
                Theta=obj.lines(i).theta;
                t1 = obj.lines(i).point1;
                t2 = obj.lines(i).point2;
                rho = obj.lines(i).rho;
                
%                 imagesc(f)
%                 hold on
%                 plot([t1(1) t2(1)]',[t1(2) t2(2)]','m')
%                 hold off
                
                if rho < 103 && rho > 98
                    obj.coordRF = [t1;t2];
                    obj.RakeAngle = 90 - Theta;
                end
            end
        end
        
        function obj = coordinateToolTip(obj)
            a=(obj.coordRF(1,2)-obj.coordRF(2,2))/(obj.coordRF(1,1)-obj.coordRF(2,1));%The slope of the rake face hardly will be Inf(Infinite) or NaN(Not-a-number),
            %because we took for this face a slope smaller than 45?
            b=obj.coordRF(1,2)-a*obj.coordRF(1,1);
            m=(obj.coordCF(1,2)-obj.coordCF(2,2))/(obj.coordCF(1,1)-obj.coordCF(2,1));%Slope of the cf, in some cases may be Inf(inclination of 90?, for example)
            h=@(x)(a*x+b);%line of the clearance face represented by f
            if m==Inf||m==-Inf%if the slope of the cf is 90? or -90?(Inf or -Inf)
                xi=obj.coordCF(1,1);%xi represents the coordinate x of the intersection(tool tip)
            else
                n=obj.coordCF(1,2)-m*obj.coordCF(1,1);
                xi=(n-b)/(a-m);
            end
            yi=h(xi);
            obj.CoordinateToolTip = [xi yi];
        end
        
        function obj = displayImageAndToolTip(obj)
            figure
            imagesc(obj.frame);
            hold on
            plot(obj.CoordinateToolTip(1),obj.CoordinateToolTip(2),'xm')
            hold off
        end
        
        function obj = pointsRFandCF(obj)
            alpha = (90 - obj.ClearanceAngle)*pi/180;
            gamma = obj.RakeAngle*pi/180;
            obj.pointRF = obj.CoordinateToolTip + 90*[-cos(gamma) sin(gamma)];
            obj.pointCF = obj.CoordinateToolTip + 90*[-cos(alpha) sin(alpha)];
        end
        
        function obj = temperatureRFandCF(obj)
%             green = [23 156 125]/255;
%             orange = [255 102 0]/255;
            pixelpitch = 15/1000;
            extCF = obj.pointCF;
            extRF = obj.pointRF;
            l1 = round(abs(obj.CoordinateToolTip(1)-extRF(1)));
            l2 = round(abs(obj.CoordinateToolTip(2)-extCF(2)));
            vRFx = round(linspace(obj.CoordinateToolTip(1),extRF(1),l1));
            vRFy = round(linspace(obj.CoordinateToolTip(2),extRF(2),l1));
            vCFx = round(linspace(obj.CoordinateToolTip(1),extCF(1),l2));
            vCFy = round(linspace(obj.CoordinateToolTip(2),extCF(2),l2));
            T_RF = zeros(1,l1);
            T_CF = zeros(1,l2);
            for t=1:l1
                T_RF(t)=obj.frame(vRFy(t),vRFx(t));
            end
            for t=1:l2
                T_CF(t)=obj.frame(vCFy(t),vCFx(t));
            end
            d1=zeros(1,l1);
            d2=zeros(1,l2);
            for t=1:l1 - 1
                d1(t+1)=(((vRFx(t+1)-vRFx(1))^2)+((vRFy(t+1)-vRFy(1))^2))^(1/2);
            end
            for t=1:l2 - 1
                d2(t+1)=(((vCFx(t+1)-vCFx(1))^2)+((vCFy(t+1)-vCFy(1))^2))^(1/2);
            end
            d1=d1*pixelpitch;
            d2=d2*pixelpitch;
            figure
            hold on
            plot(d1,T_RF)
            plot(d2,T_CF)
            xlabel('Distance from the tool tip (mm)')
            ylabel('Temperature (?C)')
            legend('Rake face','Clearance face')
            %saveas(fig,file,'jpeg')
            hold off
            figure
            imagesc(obj.frame)
            hold on
            plot(vRFx,vRFy,'k','LineWidth',1)
            plot(vCFx,vCFy,'k','LineWidth',1)
            hold off
        end
        
        function obj = displayIsotherms(obj)
            tRF = obj.RakeAngle*pi/180;
            tCF = (90 - obj.ClearanceAngle)*pi/180;
            vRF = [-cos(tRF) sin(tRF)];
            vCF = [-cos(tCF) sin(tCF)];
            %p1 RF direction
            t = (obj.CoordinateToolTip(1) - 1)/vRF(1);
            p1 = obj.CoordinateToolTip - t*vRF;
            %p2 CF direction
            t = (256 - obj.CoordinateToolTip(2))/vCF(2);
            p2 = obj.CoordinateToolTip + t*vCF;
            
            auxX = [p1(1) obj.CoordinateToolTip(1) p2(1)]';
            auxY = [p1(2) obj.CoordinateToolTip(2) p2(2)]';
           
%             figure
            Tmax = max(max(obj.biImageTool.*obj.frame));
            Tv = obj.validTemperature;
            v=round(linspace(Tmax,Tv,8));
            contour(obj.frame,v);
            hold on
            text(150,150,'Tool')
            text(100,75,'Chip')
            daspect([1,1,1])
            plot(auxX,auxY,'k')
            xlabel('pixel');
            ylabel('pixel');
            title('Isotherms');
            cb = colorbar('vert');
            zlab = get(cb,'ylabel');
            set(zlab,'String','Temperature (?C)');
            axis([50 250 25 200])
            %saveas(fig,t,'jpeg')
            colormap jet
            hold off
%             figure
%             mesh(obj.frame)
%             axis tight
%             colormap jet
%             figure
%             imagesc(obj.frame)
%             hold on
%             plot(auxX,auxY,'k')
%             daspect([1,1,1])
%             colormap jet
%             hold off
        end
        
        function obj = calculateGradient(obj)
            pp = 15*10^-6;
            tx = zeros(size(obj.frame));
            ty = zeros(size(obj.frame));
            k = 0;
            for j = 1:5
                [auxx,auxy]=gradaux_v2(obj.frame,j);
                tx = tx + auxx;
                ty = ty + auxy;
                k = k + 1;
            end
            obj.Tx = tx/(k*pp);
            obj.Ty = ty/(k*pp);
        end
        
        function obj = displayGradient(obj)
            auxx = [obj.pointCF(1) obj.CoordinateToolTip(1) obj.pointRF(1)];
            auxy = [obj.pointCF(2) obj.CoordinateToolTip(2) obj.pointRF(2)];
            k = 67;
            qx = -k*obj.Tx;
            qy = -k*obj.Ty;
            figure
            quiver(qx,qy)
            hold on
            plot(auxx,auxy,'k')
            xmin = obj.CoordinateToolTip(1) - 10;
            xmax = obj.CoordinateToolTip(1) + 5;
            ymin = obj.CoordinateToolTip(2) - 5;
            ymax = obj.CoordinateToolTip(2) + 10;
            axis([xmin xmax ymin ymax])
            title('Tool Tip')
            daspect([1,1,1])
            figure
            quiver(qx,qy)
            hold on
            plot(auxx,auxy,'k')
            xmin = obj.CoordinateToolTip(1) - 30;
            xmax = obj.CoordinateToolTip(1) - 10;
            ymin = obj.CoordinateToolTip(2) - 5;
            ymax = obj.CoordinateToolTip(2) + 15;
            axis([xmin xmax ymin ymax])
            title('Rake Face')
            daspect([1,1,1])
            figure
            quiver(qx,qy)
            hold on
            plot(auxx,auxy,'k')
            xmin = obj.CoordinateToolTip(1) - 10;
            xmax = obj.CoordinateToolTip(1) + 10;
            ymin = obj.CoordinateToolTip(2) + 10;
            ymax = obj.CoordinateToolTip(2) + 20;
            axis([xmin xmax ymin ymax])
            title('Clearance Face')
            daspect([1,1,1])
        end
        
        function obj = displayGradientContour(obj)
            auxx = [obj.pointCF(1) obj.CoordinateToolTip(1) obj.pointRF(1)];
            auxy = [obj.pointCF(2) obj.CoordinateToolTip(2) obj.pointRF(2)];
            k = 67;
            qx = -k*obj.Tx;
            qy = -k*obj.Ty;
%             numLines = 5;
%             vx = round(linspace(obj.CoordinateToolTip(1)+3,obj.CoordinateToolTip(1)+20,numLines));
%             vy = round(linspace(obj.CoordinateToolTip(2)-3,obj.CoordinateToolTip(2)-20,numLines));
%             v = zeros(1,numLines);
%             for i = 1:numLines
%                 v(i) = obj.frame(vy(i),vx(i));
%             end
            figure
            quiver(qx,qy)
            hold on
            plot(auxx,auxy,'k')
            contour(obj.frame,10)
            xmin = obj.CoordinateToolTip(1) - 20;
            xmax = obj.CoordinateToolTip(1) + 5;
            ymin = obj.CoordinateToolTip(2) - 5;
            ymax = obj.CoordinateToolTip(2) + 20;
            axis([xmin xmax ymin ymax])
            daspect([1,1,1])
        end
        
        function obj = findLineChip(obj)
            [m,n] = size(obj.frame);
            xm = obj.xyMaxTemp(1);
            ym = obj.xyMaxTemp(2);
            o = obj.RakeAngle*pi/180;
            x1 = xm - tan(o)*(ym - 1);
            x2 = x1 + tan(o)*(m - 1);
            vx = round(linspace(x1,x2,m));
            vy = linspace(1,m,m);
            B1 = zeros(m,n);
            for i = 1:m
                B1(vy(i),vx(i)) = 1;
            end
            B2 = B1 == 1 & obj.biImageChip == 1;
            obj.lineChip = B2;
%             imshow(B2)
%             figure
%             imshow(obj.biImageChip)
%             aux = obj.biImageChip - B2;
%             figure
%             imshow(aux)
        end
        
        function obj = findLineTool(obj)
             Tv = obj.validTemperature;
%------------------------------------------------------------------------
            % first method - small range of pixels             
%             Frame = obj.biImageTool.*obj.frame;
%             B1 = Frame < (Tv + 16) & Frame > (Tv + 6);
%             imshow(B1)
%             B2 = obj.biImageTool == 1 & B1 == 1;
%             obj.lineTool = B2;
%             figure
%             imshow(B2)
%------------------------------------------------------------------------
            % second method  -  removing a layer       
%             Frame = obj.biImageTool.*obj.frame > obj.validTemperature + 10;
%             aux = Frame;
%             [m,n] = size(Frame);
%             for i = 1:m
%                 for j = 1:n
%                     if Frame(i,j) == 1
%                         Frame(i,j) = 0;
%                         break
%                     end
%                 end
%             end
%             imshow(Frame)
%             B = aux - Frame;
%             imshow(B)
%             obj.lineTool = B;
%------------------------------------------------------------------------            
%             %third method - isotherms method
%             [c,~] = contour(obj.frame, [250 250]);
%             aux = round(c(:,2:end));
%             B = zeros(size(obj.frame));
%             col = size(aux,2);
%             for i = 1:col
%                 B(aux(2,i),aux(1,i)) = 1;
%             end
%             imshow(B)
%             B1 = B.*obj.biImageTool;
%             imshow(B1)
%             obj.lineTool = B1;

%------------------------------------------------------------------------
            %fourth method - two lines method
            tRF = obj.RakeAngle*pi/180;
            tCF = (90 - obj.ClearanceAngle)*pi/180;
            vRF = [-cos(tRF) sin(tRF)];
            vCF = [-cos(tCF) sin(tCF)];
            v = vRF + 2*vCF;%normal vector to the lines
            modv = (v(1)^2 + v(2)^2)^(1/2);
            v = v/modv;
            thetav = abs(atan(v(2)/v(1)));
            dt = thetav-tRF;
            ptos = zeros(2,2);
            k = 1;
            aux = obj.biImageTool.*obj.frame > Tv;
%             imshow(aux)
            for i = 80:-2:10
                p = round(obj.CoordinateToolTip + i*vRF);
                if aux(p(2),p(1)) == 1
                    ptos(k,:) = round(obj.CoordinateToolTip + i*cos(dt)*v);
                    k = k + 1;
                end
                if k == 3
                    break
                end
            end
            
            % first line 
            t = [v(2) -v(1)];%tangent vector
            p1 = ptos(1,:) + t*150;
            p2 = ptos(1,:) - t*150;
            m = round(p1(1) - p2(1) + 1);
            n = round(p1(2) - p2(2) + 1);
            if n < m
                s = m;
            else
                s = n;
            end
            vx = round(linspace(p1(1),p2(1),s));
            vy = round(linspace(p1(2),p2(2),s));
            B1 = zeros(size(obj.frame));
            for i = 1:m
                if vy(i) > 0 && vy(i) < 257 && vx(i) > 0 && vx(i) < 321
                    B1(vy(i),vx(i)) = 1;
                end
            end
            vM = obj.biImageTool.*obj.frame > Tv;
            B1f = B1 == 1 & vM == 1;
            aux(:,:,1) = B1f;
            
            %second line
            p1 = ptos(2,:) + t*150;
            p2 = ptos(2,:) - t*150;
            m = round(p1(1) - p2(1) + 1);
            n = round(p1(2) - p2(2) + 1);
            if n < m
                s = m;
            else
                s = n;
            end
            vx = round(linspace(p1(1),p2(1),s));
            vy = round(linspace(p1(2),p2(2),s));
            B2 = zeros(size(obj.frame));
            for i = 1:m
                if  vy(i) > 0 && vy(i) < 257 && vx(i) > 0 && vx(i) < 321
                    B2(vy(i),vx(i)) = 1;
                end
            end
            B2f = B2 == 1 & vM == 1;
            aux(:,:,2) = B2f;
            
%             figure
%             b = obj.biImageTool - B2f - B1f;
%             imshow(b)
%             
%             figure
%             imagesc(b.*obj.frame)
            
            obj.lineTool = aux;
        end
        
        function obj = heatBalance(obj,tuc,Vc)
            k = 67;
            pp = 15*10^-6;
            %First part - Heat carried away by the chip
            cp = [1.395833333333336e-06 -8.374999999999993e-04 0.531666666666665 4.350000000000002e+02];
            obj.heatCapacity = cp;
            M = obj.lineChip.*obj.frame;
            MH = polyval(cp,M);
            MH(MH == cp(4)) = 0;
            Ht = MH.*obj.frame;%J/kg
            Ht = sum(sum(Ht));
            n = sum(sum(obj.lineChip));
            Hc = Ht/n; %mean entalpy on the line chip
%             Vchip = 100*200/(60*n*15);
            p = 7874; %kg/m^3
            Qc = Hc*Vc*tuc*p;
            obj.HeatCarriedAwayByChip = Qc;
%-------------------------------------------------------------------------            
            %Second part - Heat carried away by the tool
            
            
            tRF = obj.RakeAngle*pi/180;
            tCF = (90 - obj.ClearanceAngle)*pi/180;
            v = [-(cos(tRF)+cos(tCF)) (sin(tCF)+sin(tRF))];%normal vector to the lines
            modv = ((v(1))^2 + (v(2))^2)^(1/2);
            v = v/modv;
            %first line
            aux = obj.lineTool(:,:,1);
            L = sum(sum(aux));
            Gx = sum(sum((aux.*obj.Tx)*v(1)));
            Gy = sum(sum((aux.*obj.Ty)*v(2)));
            G = (Gx + Gy)/L;
            Qt1 = -k*L*pp*G;
            
            %second line
            aux = obj.lineTool(:,:,2);
            L = sum(sum(aux));
            Gx = sum(sum((aux.*obj.Tx)*v(1)));
            Gy = sum(sum((aux.*obj.Ty)*v(2)));
            G = (Gx + Gy)/L;
            Qt2 = -k*L*pp*G;
            
            Qm = (Qt1 + Qt2)/2;
            obj.HeatFluxAwayFromToolTip = Qm;
            
%             % Second part - heat (by tool) method isotherms
%             ax = obj.lineTool.*obj.Tx;
%             ay = obj.lineTool.*obj.Ty;
%             B = obj.frame.*(ones(size(obj.frame)) - obj.lineTool);
%             close all
%             imagesc(B)
        end
        
        function n = exceedingPoints(obj,Temperature)
            B = obj.frame.*obj.biImageTool > Temperature;
            n = sum(sum(B));
%             K = sprintf('Number of points exceeding the defined temperature is: %3d \n',n);
%             disp(K)
        end
        
        function obj = internalEnergyTool(obj)
            pp = 15*10^-6;
            p = 7874;%density kg/m^3
            Te = 22;
            B = obj.frame.*obj.biImageTool > obj.validTemperature;
            n = sum(sum(B));
            B1 = obj.frame.*B;
            B2 = polyval(obj.heatCapacity,B1);
            B2(B2 == obj.heatCapacity(4)) = 0;
            H = B2.*(obj.frame - Te)*p*pp^2;%Heat Amount (J/m^3)
            Ha = sum(sum(H))/n;
            obj.InternalEnergyTool = Ha;
        end
        
        function B = passBinaryImageTool(obj)
            B = obj.biImageTool;
        end
       
        function B = passBinaryImageChip(obj)
            B = obj.biImageChip;
        end
        
        function obj = shearLine(obj)
            B = obj.biImageChip;
            v1 = sum(B);
            v1(v1 == 0) = [];
            l1 = length(v1);
            C = imcrop(B,[20 20 l1 100]);
            [m,n] = size(C);
            pto = zeros(1000,2);
            count = 1;
            for i = 2:m-1
                for j = 2:n-1
                    if C(i,j+1) == 1 && C(i,j-1) == 1 && C(i+1,j) == 1 && C(i-1,j) == 1
                        pto(count,:) = [i j];
                        count = count + 1;
                    end
                end
            end
            for i =1000:-1:1
                if isequal(pto(i,:),[0 0]) == 1
                    pto(i,:) = [];
                end
            end
            l = size(pto,1);
            for i = 1:l
                C(pto(i,1),pto(i,2)) = 0;
            end
            
            [H, THETA, RHO] = hough(C,'Theta',-40:-30);%Hough transformation
            P  = houghpeaks(H, 5);
            lin = houghlines(C, THETA, RHO, P, 'FillGap', 15,'MinLength',10);
            l=length(lin);
            p1 = [];
            p2 = [];
            for i=1:l
                Theta=lin(i).theta;
                t1 = lin(i).point1;
                t2 = lin(i).point2;
                y = abs(t1(2)-t2(2));
%                 imshow(C)
%                 hold on
%                 plot([t1(1) t2(1)]',[t1(2) t2(2)]','xm')
%                 hold off
                if isempty(p1) && isempty(p2) && abs(Theta + 34) < 5
                    p1 = t1 + [19 19];
                    p2 = t2 + [19 19];
                    ym = y;
                    obj.ShearAngle = abs(Theta);
                end
                if abs(Theta + 34) < 5 && y > ym
                    p1 = t1 + [19 19];
                    p2 = t2 + [19 19];
                    obj.ShearAngle = abs(Theta);
                end
            end
            
        end
        
        function obj = forcesValues(obj,Fp,Fq,w,tuc)
            phi = obj.ShearAngle*pi/180;%shear angle
            gamma = obj.RakeAngle*pi/180;%Rake angle
            Fs = Fp*cos(phi) - Fq*sin(phi);%Cutting force component parallel to shear plane
            Ns = Fq*cos(phi) + Fp*sin(phi);%Cutting force component perpendicular to shear plane
            Fc = Fp*sin(gamma) + Fq*cos(gamma);%Cutting force component parallel to tool face
            Nc = Fp*cos(gamma) - Fq*sin(gamma);%Cutting force component perpendicular to tool face
            mu = Fc/Nc; % coefficient of friction
            As = w*tuc/sin(phi);%Area shear plane
            tau = Fs/As;%shear stress
            sigma = Ns/As;%Normal stress
            r = sin(phi)/cos(phi - gamma); %ratio r = t/tc = lc/l
            ss = cos(gamma)/(sin(phi)*cos(phi-gamma));%shear strain
            us = tau*ss;%shear energy per volume
            uf = Fc*r/(tuc*w);%friction energy per volume
            beta = atan(Fc/Nc);%friction angle on tool face
            obj.CuttingForceParallelToolFace = Fc;
            obj.CuttingForcePowerDirection = Fp;
            obj.CuttingForceUncutChipThicknessDirection = Fq;
            obj.CuttingForceParallelShearPlane = Fs;
            obj.CuttingForcePerpendicularShearPlane = Ns;
            obj.CuttingForcePerpendicularToolFace = Nc;
            obj.CoefficientFriction = mu;
            obj.ShearStress = tau;
            obj.NormalStress = sigma;
            obj.RatioR = r;
            obj.ShearEnergyVolume = us;
            obj.FrictionEnergyVolume = uf;
            obj.FrictionAngle = beta*180/pi;
        end
        
        function s = resultsStruct(obj)
            f1 = 'CoordinateToolTip';
            value1 = {'-',obj.CoordinateToolTip,'pixel'};%eval(['obj.' f1])
            f2 = 'TemperatureToolTip';
            value2 = {'-',eval(['obj.' f2]),'?C'};
            f3 = 'RakeAngle';
            value3 = {'gamma',eval(['obj.' f3]),'?'};
            f4 = 'ClearanceAngle';
            value4 = {'alpha',eval(['obj.' f4]),'?'};
            f5 = 'ShearAngle';
            value5 = {'phi',eval(['obj.' f5]),'?'};
            f6 = 'FrictionAngle';
            value6 = {'beta',eval(['obj.' f6]),'?'};
            f7 = 'MeanTemperatureTool';
            value7 = {'-',eval(['obj.' f7]),'?C'};
            f8 = 'MaximumTemperatureChip';
            value8 = {'-',eval(['obj.' f8]),'?C'};
            f9 = 'MaximumTemperatureCuttingZone';
            value9 = {'-',eval(['obj.' f9]),'?C'};
            f10 = 'HeatCarriedAwayByChip';
            value10 = {'q_c',eval(['obj.' f10]),'W/m'};
            f11 = 'HeatFluxAwayFromToolTip';
            value11 = {'q_t',eval(['obj.' f11]),'W/m'};
            f12 = 'InternalEnergyTool';
            value12 = {'U''',eval(['obj.' f12]),'J/m'};
            f13 = 'CuttingForceParallelToolFace';
            value13 = {'Fc',eval(['obj.' f13]),'N'};
            f14 = 'CuttingForcePowerDirection';
            value14 = {'Fp',eval(['obj.' f14]),'N'};
            f15 = 'CuttingForceUncutChipThicknessDirection';
            value15 = {'Fq',eval(['obj.' f15]),'N'};
            f16 = 'CuttingForceParallelShearPlane';
            value16 = {'Fs',eval(['obj.' f16]),'N'};
            f17 = 'CuttingForcePerpendicularShearPlane';
            value17 = {'Ns',eval(['obj.' f17]),'N'};
            f18 = 'CuttingForcePerpendicularToolFace';
            value18 = {'Nc',eval(['obj.' f18]),'N'};
            f19 = 'CoefficientFriction';
            value19 = {'mu',eval(['obj.' f19]),'-'};
            f20 = 'ShearStress';
            value20 = {'tau',eval(['obj.' f20]),'N/m^2'};
            f21 = 'NormalStress';
            value21 = {'sigma',eval(['obj.' f21]),'N/m^2'};
            f22 = 'RatioR';
            value22 = {'r',eval(['obj.' f22]),'-'};
            f23 = 'ShearEnergyVolume';
            value23 = {'us',eval(['obj.' f23]),'J/m^3'};
            f24 = 'FrictionEnergyVolume';
            value24 = {'uf',eval(['obj.' f24]),'J/m^3'};
            f25 = 'CuttingVelocity';
            value25 = {'Vp',eval(['obj.' f25]),'m/s'};
            f26 = 'TotalPowerProduced';
            value26 = {'Wt',eval(['obj.' f26]),'W/m'};
            f27 = 'UnCutChipThickness';
            value27 = {'tuc',eval(['obj.' f27]),'m'};
            s = struct(f1,value1,f2,value2,f3,value3,f4,value4,...
                f5,value5,f6,value6,f7,value7,f8,value8,f9,value9,...
                f10,value10,f11,value11,f12,value12,f13,value13,f14,value14,...
                f15,value15,f16,value16,f17,value17,f18,value18,f19,value19,...
                f20,value20,f21,value21,f22,value22,f23,value23,f24,value24,...
                f25,value25,f26,value26,f27,value27);
        end
    end
end