%% Loading the file and initializing the variables

% clear all
% close all

pixelpitch=15/1000;%in mm
%The following points represent the position(the pixel) of the tool tip and of
%a point on the rake face's line
%Coordinates of the tool tip, rake face and clearance face

imagesc(Frame)
hold on
axis tight
daspect([1,1,1])
%We will plot the line on the rake face's profile afterwords

T1=[];%Vector of temperature - rake face
T2=[];%Vector of temperature - clearance face

%% Findind the RF(rake face) and the CF(clearance face)
%The variables m,n,a,b and(xir,yir) we already got from the script
%(CuttingEdge)
%For CF
if m==Inf||m==-Inf
    auxCF=[xir yir-170];%auxCR is an auxiliar variable to get the clearance face
else
    v=[1+m^2 -2*(xir+m*yir-m*n) n^2-2*yir*n+yir^2+xir^2-170^2];%The number 170 means the length of the line on the faces
    r1=roots(v);
    g=@(x)(m*x+n);
    r2=g(r1);
    if r2(1)<yir
        auxCF=[r1(1) r2(1)];
    else
        auxCF=[r1(2) r2(2)];
    end
end
%For RF
v=[1+a^2 -2*(xir+a*yir-a*b) b^2-2*yir*b+yir^2+xir^2-170^2];
r1=roots(v);
r2=h(r1);
if r1(1)>xir
    auxRF=[r1(1) r2(1)];%auxRF is an auxiliar variable to get the rake face
else
    auxRF=[r1(2) r2(2)];
end
%The numerical order is first line:Coordinates of the clearance face
%second line: Coordinates of tool tip
%Third line: Coordinates of the rake face
Coord=[];
Coord(1,:)=auxCF;
Coord(2,:)=[xir yir];
Coord(3,:)=auxRF;
x=round(Coord(:,1));%x coordinates (first column)
y=round(Coord(:,2));%y coordinates (second column)
%The coordinates were rounded because we are working with pixels
%then, the numbers must be integers

%% Code to plot the rake temperature of the rake face

sizex=abs(x(2)-x(3)+1);
sizey=abs(y(2)-y(3)+1);
if sizex < sizey
    sizef=sizey;
else
    sizef=sizex;
end
%Pixel of the rake face's line
rangex=round(linspace(x(2),x(3),sizef));
rangey=round(linspace(y(2),y(3),sizef));
plot(rangex,rangey,'w','LineWidth',1)


for t=1:sizef
    T1(t)=Frame(rangey(t),rangex(t));
end
d1=zeros(1,sizef);
for t=1:sizef-1
    d1(t+1)=(((rangex(t+1)-rangex(1))^2)+((rangey(t+1)-rangey(1))^2))^(1/2);
end
%% Code to plot the temperature of the clearance face

sizex=abs(x(1)-x(2)+1);
sizey=abs(y(1)-y(2)+1);
if sizex < sizey
    sizef=sizey;
else
    sizef=sizex;
end
%Pixels of the clearance face's line
rangex=round(linspace(x(2),x(1),sizef));
rangey=round(linspace(y(2),y(1),sizef));
plot(rangex,rangey,'w','LineWidth',1)
hold off

for t=1:sizef
    T2(t)=Frame(rangey(t),rangex(t));
end
d2=zeros(1,sizef);
for t=1:sizef-1
    d2(t+1)=(((rangex(t+1)-rangex(1))^2)+((rangey(t+1)-rangey(1))^2))^(1/2);
end
%% Ploting

d1=d1*pixelpitch;
d2=d2*pixelpitch;
fig=figure;
hold on
plot(d1,T1)
plot(d2,T2,'r')
xlabel('Distance from the tool tip (mm)')
ylabel('Temperature (°C)')
legend('Rake face','Clearance face')
%saveas(fig,file,'jpeg')
hold off

%% Draw the isoterms

fig=figure;
Region=imcrop(Frame,[275 1 306 287]);
v=round(linspace(Region(214,65),Region(132,165),8));
imagesc(Region)
xp=x-275;
yp=y;
contour(Region,v);
hold on
daspect([1,1,1])
plot(xp,yp,'k')
xlabel('pixel');
ylabel('pixel');
t=['Isotherms - Frame' f];
title(t);
cb = colorbar('vert');
zlab = get(cb,'ylabel');
set(zlab,'String','Temperature (°C)');
%saveas(fig,t,'jpeg')
hold off
