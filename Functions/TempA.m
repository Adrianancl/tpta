%% Loading the file and initializing the variables
clear all
close all

file=input('What is the file''s name? ','s');
V=load(file,file);
aux=V.(file);
[n,m,l]=size(aux);
nf=num2str(l);
frame=input(['What is the frame? (Frames avaliable:' nf ') '],'s');
f=str2num(frame);
Frame=aux(:,:,f);

pixelpitch=15/1000;%in mm
%The following points represent the position(the pixel) of the tool tip and of
%a point on the rake face's line
%Coordinates of the tool tip, rake face and clearance face
if file(1:5)=='AYA01'
    x=[246;258;329];
    y=[67;152;82];
elseif file(1:5)=='AYA03'
    x=[265;278;327];
    y=[53;131;58];
elseif file(1:5)=='AYA05'
    x=[258;268;325];
    y=[59;137;69];
elseif file(1:5)=='AYA06'
    x=[257;271;362];
    y=[50;149;70];
elseif file(1:5)=='AYB01'
    x=[260;273;344];
    y=[63;150;72];
elseif file(1:5)=='AYB03'
    x=[260;276;362];
    y=[60;144;82];
elseif file(1:5)=='AYB05'
    x=[255;275;349];
    y=[56;137;77];
else
    x=[259;273;349];
    y=[50;137;70];
end

imagesc(Frame)
hold on
axis tight
daspect([1,1,1])
%We will plot the line on the rake face's profile afterwords

T1=[];%Vector of temperature - rake face
T2=[];%Vector of temperature - clearance face
%% Code to plot the rake temperature of the rake face

sizex=abs(x(1)-x(2)+1);
sizey=abs(y(1)-y(2)+1);
if sizex < sizey
    sizef=sizey;
else
    sizef=sizex;
end
%Pixel of the rake face's line
rangex=round(linspace(x(1),x(2),sizef));
rangey=round(linspace(y(1),y(2),sizef));
plot(rangex,rangey,'w','LineWidth',1)


for t=1:sizef
    T1(t)=Frame(rangey(t),rangex(t));
end
d1=zeros(1,sizef);
for t=1:sizef-1
    d1(t+1)=(((rangex(t+1)-rangex(1))^2)+((rangey(t+1)-rangey(1))^2))^(1/2);
end
%% Code to plot the temperature of the clearance face

sizex=abs(x(1)-x(3)+1);
sizey=abs(y(1)-y(3)+1);
if sizex < sizey
    sizef=sizey;
else
    sizef=sizex;
end
%Pixels of the clearance face's line
rangex=round(linspace(x(1),x(3),sizef));
rangey=round(linspace(y(1),y(3),sizef));
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
fig=figure(2);
hold on
plot(d1,T1)
plot(d2,T2,'r')
xlabel('Distance from the tool tip (mm)')
ylabel('Temperature (°C)')
legend('Rake face','Clearance face')
saveas(fig,file,'jpeg')
hold off
%% Image of the tool
Region=[];
for t=x(1):x(3)
    for i=y(1):y(2)
        Region(i-y(1)+1,t-x(1)+1)=Frame(i,t);
    end
end
figure(3)
imagesc(Region)
hold on
daspect([1,1,1])
xn=(x-(x(1)-1)*ones(3,1));
yn=(y-(y(1)-1)*ones(3,1));
xp=[xn(2);xn(1);xn(3)];
yp=[yn(2);yn(1);yn(3)];
plot(xp,yp,'k')
xlabel('pixel');
ylabel('pixel');
disp(' ');
T=Region(yn(1),xn(1));
K=sprintf('Temperature on the tool tip: %3.2f °C\n',T );
disp(K)
hold off


%% Draw the isoterms

%[TX,TY] = gradient(Region,5);
%quiver(TX,TY)
f=figure(4);
[C,h]=contour(Region,20);
%clabel(C)
hold on
daspect([1,1,1])
plot(xp,yp,'k')
xlabel('pixel');
ylabel('pixel');
t=['Isotherms - ' file ' - Frame ' frame];
title(t);
cb = colorbar('vert');
zlab = get(cb,'ylabel');
set(zlab,'String','Temperature (°C)');
saveas(f,t,'jpeg')
hold off
