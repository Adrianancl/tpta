t=0.0057; %time interval(seconds) 353,1/2 Hz  (2 = Number of presets)
k=88.4; %material conductivity W/m.?C or W/m.K
pixelpitch = 15/10^6;%in meters
x = zeros(3,2);
x(2,:) = round([301 235]);%mean of the  TT's points during the cut
x(3,:) = round([x(2,1) x(2,2)]+150*[cos(h.ThetaRF) -sin(h.ThetaRF)]);%Here
%we assume that the RF has a 150 pixels length 
x(1,:) = round([x(2,1) x(2,2)]-5*[-cos(h.ThetaCF) sin(h.ThetaCF)]);
l1 = abs(x(3,1)-x(2,1)+1);% the coordinate x of the RF is greater than TT
l2 = abs(-x(1,2)+x(2,2)+1);% the coordinate y of the TT is greater than CF
HeatRF = zeros(1,l1);
HeatCF = zeros(1,l2);
for i=18:2:72  %these are the frames that the cutting is happening 
    aux = num2str(i);
    Frame = eval(['Frame' aux]);
    Tx = zeros(size(Frame));
    Ty = zeros(size(Frame));
    count = 0;
    for j = 10
        [tx,ty]=gradaux(Frame,j);
        Tx = Tx + tx;
        Ty = Ty + ty;
        count = count + 1;
    end
    Tx = Tx/(count*pixelpitch);
    Ty = Ty/(count*pixelpitch);
        
    qxt=-k*Tx;%x direction
    qyt=-k*Ty;%y direction
    
    imagesc(Frame)
    hold on
    quiver(qxt,qyt)
    plot(x(:,1),x(:,2),'k')
    hold off
    
    % Heat flux through the rake face
    qrt = -(qxt*sin(h.ThetaRF)+qyt*cos(h.ThetaRF));%heat flux resultant (W/m^2)
    %normal to the face
    rangex = round(linspace(x(2,1),x(3,1),l1));
    rangey = round(linspace(x(2,2),x(3,2),l1));
    %equally distant from the tool tip in that position
    n = 10;%number of pixels into the tool(heat flux)
    for j = 1:l1
        auxH = 0;
        for k = 1:n
            px = round(rangex(j) + (k - 1)*(-sin(h.ThetaRF)));
            py = round(rangey(j) + (k - 1)*(-cos(h.ThetaRF)));
            auxH = auxH + qrt(py,px);
        end
        HeatRF(j) = HeatRF(j) + auxH/n;
    end
    
    % Heat flux through the clearance face
    % What is the effective length of the clearance face? Assume = 5
    % pixels (75um or 0,075mm)
    
    qrt = qxt*sin(h.ThetaCF)+qyt*cos(h.ThetaCF);%heat flux resultant (W/m2)    
    rangex = round(linspace(x(2,1),x(1,1),l2));
    rangey = round(linspace(x(2,2),x(1,2),l2));
    for j = 1:l2
        auxH = 0;
        for k = 1:n
            px = round(rangex(j) + (k - 1)*(sin(h.ThetaCF)));
            py = round(rangey(j) + (k - 1)*(cos(h.ThetaCF)));
            auxH = auxH + qrt(py,px);
        end
        HeatCF(j) = HeatCF(j) + auxH/n;
    end
end
meanHeatRF = HeatRF/28; %W/m2 average of heat flux in each pixel of the RF
meanHeatCF = HeatCF/28;
TotHeatRF = sum(meanHeatRF)/l1%Average of heat flux on the RF's line
TotHeatCF = sum(meanHeatCF)/l2
TotHeatFlux = (sum(meanHeatRF) + sum(meanHeatCF))/(l1 + l2)%Average of heat flux on the cutting edge
