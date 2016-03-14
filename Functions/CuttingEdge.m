%% Initializing the frame that we wanna analyze
f = input('Which frame do you want to analyze? ');%string input(we must choose a number of frame that exists in the workspace)
f=num2str(f);
Frame = eval(['Frame' f]);%Here we build a string with the same name of the variable that we wanna analyze, then with the eval function
%we attribute to the Frame(variable) the matrix values
%% Finding the cutting edge
%In this section we gonna find the edges of the tool. We have to use the
%functions 'hough','houghpeaks' and 'houghlines' combined to get what we
%want. The last function houghlines will give to us the initial and final
%points of the lines that the function got(in our case, the lines of the
%clearance and rake faces)
cropped=imcrop(Frame,[275 1 306 287]);%We will have to shift the x axis afterwards
%because we are cutting the image (with this function crop), what makes the axis be redefined
BW = edge(cropped,'sobel');%Here we make a binary image that makes us able to see the cutting edge
imshow(BW)%Displays the image
[H, THETA, RHO] = hough(BW);%Hough transformation
P  = houghpeaks(H,3);%Identify peaks in the hough transformation (Helps to find the lines of the cutting edge)
%2 is the maximum number of peaks that we want identify
lines = houghlines(BW,THETA,RHO,P,'FillGap',25,'MinLength',30);%Here we can find the lines of cutting edge and afterwards find the coordinate of the tool tip
%% Extracting the coordinates
l=length(lines);%length of the struct lines that we gonna use to take the coordinates
%of the rake face and the clearance face

rf=[];%empty vectors (rf=rake face/cf=clearance face)
cf=[];
%The following loop was written to get the coordinates of the
%clearance/rake face, because the function houghlines may give to us more coordinates
%than we need, so here we make sure that we won't take the same lines for
%the cf and rf(difference between the slope of these lines)
for i=1:l
    t1=lines(i).point1;
    t2=lines(i).point2;
    tg=(abs(t1(2)-t2(2)))/(abs(t1(1)-t2(1)));%slope
    if tg > 1%For the cf the slope must be bigger than 45°(tg45=1)
        cf=[t1;t2];
        thetaCF=atan(tg);
    elseif isempty(rf)&&tg>0&&tg<1%If the slope is smaller than 45°, it should be the rf
        rf=[t1;t2];
        thetaRF=atan(tg);
    end
    if isempty(rf)==0&&isempty(cf)==0
        break
    end
end

%% Finding the coordinate of the tool tip
%Here we gonna build two lines as functions of the first degree, one for the
%rake face and other for the clearance face, then we can find the
%intersetion(tool tip)
%Rake face: f(x)=a*x+b
%Clearance face: f(x)=m*x+n
%The intersection is a*x+b=m*x+n->x=(n-b)/(a-m)
a=(rf(1,2)-rf(2,2))/(rf(1,1)-rf(2,1));%The slope of the rake face hardly will be Inf(Infinite) or NaN(Not-a-number),
%because we took for this face a slope smaller than 45°
b=rf(1,2)-a*rf(1,1);
m=(cf(1,2)-cf(2,2))/(cf(1,1)-cf(2,1));%Slope of the cf, in some cases may be Inf(inclination of 90°, for example)
h=@(x)(a*x+b);%line of the clearance face represented by f
if m==Inf||m==-Inf%if the slope of the cf is 90° or -90°(Inf or -Inf)
    xi=cf(1,1);%xi represents the coordinate x of the intersection(tool tip)
    %of the cropped image(x is shiftted by 275)   
else
    n=cf(1,2)-m*cf(1,1);
    xi=(n-b)/(a-m);
end
yi=h(xi);%Here we define the function to find the y coordinate of the tool tip
%Tool tip = (xi,yi)
%% Plotting in the binary image
imshow(BW)%This comand displays binary images and others formats, like .tif
hold on %Comand to make the following plots in the same screen
plot(rf(:,1),rf(:,2),'xg',cf(:,1),cf(:,2),'xy',xi,yi,'xm')
hold off
%% Plotting in the real image
figure%Open another screen
imagesc(Frame)%This comand displays matrices(And makes a scale with the values inside)
hold on
xir=xi+275;%Here we shift the x axis because we cropped the image before(xir=xi real)
%and the coordinates that we found do not correspond to the initial image
rfr=rf+[275 0;275 0];
cfr=cf+[275 0;275 0];
plot(rfr(:,1),rfr(:,2),'xg','LineWidth',2)
plot(cfr(:,1),cfr(:,2),'xy','LineWidth',2)
plot(xir,yi,'xm','LineWidth',2)
hold off
%Temperature of the tool tip
%back to the real scale before the crop

v=round([xir yi]);%Coordinates of the tool tip, must be integers (the image have its scale in pixels)
T=Frame(v(2),v(1));%Frame is a matrix, then the coordinate x is the number of the column
%and the coordinate y is the number of the line
K=sprintf('Temperature on the tool tip: %3.2f °C\n',T );
disp(K)
%% Finding the functions(lines in the real image)
%Here we have the same code of the section 'Finding the coordinate of the tool tip'
%But now for the scale of the real image (not the cropped one)
a=(rfr(1,2)-rfr(2,2))/(rfr(1,1)-rfr(2,1));%The slope of the rake face hardly will be Inf(Infinite) or NaN(Not-a-number),
%because we took for this face a slope smaller than 45°
b=rfr(1,2)-a*rfr(1,1);
m=(cfr(1,2)-cfr(2,2))/(cfr(1,1)-cfr(2,1));%Slope of the cf, in some cases may be Inf(inclination of 90°, for example)
h=@(x)(a*x+b);%line of the clearance face represented by f
if m==Inf||m==-Inf%if the slope of the cf is 90° or -90°(Inf or -Inf)
    xir=cfr(1,1);  
else
    n=cfr(1,2)-m*cfr(1,1);
    xir=(n-b)/(a-m);
end
yir=h(xir);
%% Calling the temperature analyzes function
TempA_v2


