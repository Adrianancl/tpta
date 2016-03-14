Region=zeros(300,351);
for k=18:2:72  
    i=num2str(k);
    Frame=eval(['Frame' i]);
    Region=Frame(1:300,250:600);
%     v=round(linspace(Frame(250,300),Frame(157,396),8));
%     contour(Region,v);
    mesh(Region)
    daspect([1,1,1])
%     pause(0.3);
    waitforbuttonpress;
end