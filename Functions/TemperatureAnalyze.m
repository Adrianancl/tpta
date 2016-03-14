   classdef TemperatureAnalyze
   % This class was built to analyze the temperature,  

       properties(GetAccess = 'public', SetAccess = 'private')
       % define the properties of the class here, (like fields of a struct)
           coordRF;
           coordCF;
           coordTT;
           TemperatureTT;
           ThetaRF;
           ThetaCF;
       end
       
       properties(GetAccess = 'private', SetAccess = 'private')
            cropped;
            BW;
            lines;
            frame;
            pointCF;%auxiliar to plot the cutting edge
            pointRF;
            Tx;%auxiliar to plot the gradients of the frame
            Ty;
       end

       methods
       % methods, including the constructor are defined in this block

           function obj = TemperatureAnalyze(Frame)%constructor
              obj.frame = Frame;  
              obj.cropped = imcrop(Frame,[275 1 306 287]);%here we cut the image to make it clearer
              obj.BW = edge(obj.cropped,'sobel');
              [H, THETA, RHO] = hough(obj.BW);%Hough transformation
              P  = houghpeaks(H, 3);
              obj.lines = houghlines(obj.BW, THETA, RHO, P, 'FillGap', 25,'MinLength',30);%Here we can find the lines of cutting edge and afterwards find the coordinate of the tool tip
              obj = obj.calculateCoordinates();
              if isempty(obj.coordRF)==0&&isempty(obj.coordCF)==0
                  obj = obj.coordinateToolTip();
              else %Default conditions
                  obj.coordTT = [310 239];%Average of the TT coordinates
                  obj.ThetaCF = 1.5242;%mean slope
              end
              obj = obj.pointsRFandCF();
              aux = round(obj.coordTT);
              obj.TemperatureTT = Frame(aux(2),aux(1));
              obj = obj.calculateGradient();
           end

           function obj = displayBinary(obj)
               imshow(obj.BW);
               hold on
               plot(obj.coordRF(:,1)-275,obj.coordRF(:,2),'bx')
               plot(obj.coordCF(:,1)-275,obj.coordCF(:,2),'yx')
               hold off
           end
           
           
           function obj = calculateCoordinates(obj)

                l=length(obj.lines);%length of the struct lines that we gonna use to take the coordinates
                %of the rake face and the clearance face
                %The following loop was written to get the coordinates of the
                %clearance/rake face, because the function houghlines may give to us more coordinates
                %than we need, so here we make sure that we won't take the same lines for
                %the cf and rf(difference between the slope of these lines)
                for i=1:l
                    t1=obj.lines(i).point1;
                    t2=obj.lines(i).point2;
                    tg=(abs(t1(2)-t2(2)))/(abs(t1(1)-t2(1)));%slope
                    if tg > 1%For the cf the slope must be bigger than 45°(tg45=1)
                        if isempty(obj.coordCF)
                            obj.coordCF=[t1;t2]+[275 0;275 0];%shift the coordinates because of the crop
                            obj.ThetaCF = atan(tg);%Here we have to assure that we are taking the best coordinate for the CF(others conditions)
                        end                        
                    elseif isempty(obj.coordRF)&&tg>0&&tg<1%If the slope is smaller than 45°,then it should be the rf
                        obj.coordRF=[t1;t2]+[275 0;275 0];
                        obj.ThetaRF = atan(tg);
                    end
%                     if isempty(obj.coordCF)==0&&isempty(obj.coordRF)==0%If the coordinates are filled, break the loop
%                         break
%                     end
                end
           end
           
           function obj = coordinateToolTip(obj)
                a=(obj.coordRF(1,2)-obj.coordRF(2,2))/(obj.coordRF(1,1)-obj.coordRF(2,1));%The slope of the rake face hardly will be Inf(Infinite) or NaN(Not-a-number),
                %because we took for this face a slope smaller than 45°
                b=obj.coordRF(1,2)-a*obj.coordRF(1,1);
                m=(obj.coordCF(1,2)-obj.coordCF(2,2))/(obj.coordCF(1,1)-obj.coordCF(2,1));%Slope of the cf, in some cases may be Inf(inclination of 90°, for example)
                h=@(x)(a*x+b);%line of the clearance face represented by f
                if m==Inf||m==-Inf%if the slope of the cf is 90° or -90°(Inf or -Inf)
                    xi=obj.coordCF(1,1);%xi represents the coordinate x of the intersection(tool tip)
                    %of the cropped image(x is shiftted by 275)   
                else
                    n=obj.coordCF(1,2)-m*obj.coordCF(1,1);
                    xi=(n-b)/(a-m);
                end
                yi=h(xi);
                obj.coordTT = [xi yi];
           end
           
           function obj = displayImageAndToolTip(obj)
                figure
                imagesc(obj.frame);
                hold on
                plot(obj.coordTT(1),obj.coordTT(2),'xm')
                hold off
           end
           
           function obj = pointsRFandCF(obj)
                obj.pointCF = obj.coordTT + 170*[cos(obj.ThetaRF) -sin(obj.ThetaRF)];
                obj.pointRF = obj.coordTT - 170*[-cos(obj.ThetaCF) sin(obj.ThetaCF)];           
           end
           
           function obj = temperatureRFandCF(obj)
                pixelpitch = 15/1000;
                extCF = obj.pointCF;
                extRF = obj.pointRF;
                l1 = round(abs(obj.coordTT(1)-extRF(1)));
                l2 = round(abs(obj.coordTT(2)-extCF(2)));
                vRFx = round(linspace(obj.coordTT(1),extRF(1),l1));
                vRFy = round(linspace(obj.coordTT(2),extRF(2),l1));
                vCFx = round(linspace(obj.coordTT(1),extCF(1),l2));
                vCFy = round(linspace(obj.coordTT(2),extCF(2),l2));
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
                plot(d2,T_CF,'r')
                xlabel('Distance from the tool tip (mm)')
                ylabel('Temperature (°C)')
                legend('Rake face','Clearance face')
                %saveas(fig,file,'jpeg')
                hold off
                figure
                imagesc(obj.frame)
                hold on
                plot(vRFx,vRFy,'LineWidth',1)
                plot(vCFx,vCFy,'LineWidth',1)
                hold off
           end
           
           function obj = displayIsotherms(obj)
                auxX = [obj.pointCF(1) obj.coordTT(1) obj.pointRF(1)]';
                auxY = [obj.pointCF(2) obj.coordTT(2) obj.pointRF(2)]';
                auxX2 = auxX-275;
                figure
                Region = imcrop(obj.frame,[275 1 225 287]);
                v=round(linspace(Region(214,65),Region(132,165),8));
                contour(Region,v);
                hold on
                daspect([1,1,1])
                plot(auxX2,auxY,'k')
                xlabel('pixel');
                ylabel('pixel');
                title('Isotherms');
                cb = colorbar('vert');
                zlab = get(cb,'ylabel');
                set(zlab,'String','Temperature (°C)');
                %saveas(fig,t,'jpeg')
                colormap jet
                hold off
                figure
                mesh(Region)
                axis tight
                colormap jet
                figure
                imagesc(obj.frame)
                hold on
                plot(auxX,auxY,'k')
                daspect([1,1,1])
                colormap jet
                hold off
           end
           
           function obj = calculateGradient(obj)
                tx = zeros(size(obj.frame));
                ty = zeros(size(obj.frame));
                for j = 10:30
                    [auxx,auxy]=gradient(obj.frame,j);
                    tx = tx + auxx;
                    ty = ty + auxy;
                end
                obj.Tx = tx;
                obj.Ty = ty;
           end
           
           function obj = displayGradient(obj)
                auxx = [obj.pointCF(1) obj.coordTT(1) obj.pointRF(1)];
                auxy = [obj.pointCF(2) obj.coordTT(2) obj.pointRF(2)];
                k = 88.4;
                qx = -k*obj.Tx;
                qy = -k*obj.Ty;
                figure
                quiver(qx,qy)
                hold on
                plot(auxx,auxy,'k')
                xmin = obj.coordTT(1) - 5;
                xmax = obj.coordTT(1) + 20;
                ymin = obj.coordTT(2) - 20;
                ymax = obj.coordTT(2) + 5;
                axis([xmin xmax ymin ymax])
                title('Tool Tip')
                daspect([1,1,1])
                figure
                quiver(qx,qy)
                hold on
                plot(auxx,auxy,'k')
                xmin = obj.coordTT(1) + 95;
                xmax = obj.coordTT(1) + 125;
                ymin = obj.coordTT(2) - 20;
                ymax = obj.coordTT(2) + 5;
                axis([xmin xmax ymin ymax])
                title('Rake Face')
                daspect([1,1,1])
                figure
                quiver(qx,qy)
                hold on
                plot(auxx,auxy,'k')
                xmin = obj.coordTT(1) - 5;
                xmax = obj.coordTT(1) + 20;
                ymin = obj.coordTT(2) - 80;
                ymax = obj.coordTT(2) - 60;
                axis([xmin xmax ymin ymax])
                title('Clearance Face')
                daspect([1,1,1])
           end
           
           function obj = displayGradientContour(obj)
                auxx = [obj.pointCF(1) obj.coordTT(1) obj.pointRF(1)];
                auxy = [obj.pointCF(2) obj.coordTT(2) obj.pointRF(2)];
                k = 88.4;
                qx = -k*obj.Tx;
                qy = -k*obj.Ty;
                numLines = 5;
                vx = round(linspace(obj.coordTT(1)+3,obj.coordTT(1)+20,numLines));
                vy = round(linspace(obj.coordTT(2)-3,obj.coordTT(2)-20,numLines));
                v = zeros(1,numLines);
                for i = 1:numLines
                    v(i) = obj.frame(vy(i),vx(i));
                end
                figure
                quiver(qx,qy)
                hold on
                plot(auxx,auxy,'k')
                contour(obj.frame,v)
                xmin = obj.coordTT(1) - 5;
                xmax = obj.coordTT(1) + 20;
                ymin = obj.coordTT(2) - 20;
                ymax = obj.coordTT(2) + 5;
                axis([xmin xmax ymin ymax])
                daspect([1,1,1])               
           end
       end
   end