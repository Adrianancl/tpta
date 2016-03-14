x=[255;275;349];
y=[56;137;77];


Region=[];
for k=1:39   
    for t=x(1):x(3)
            for i=y(1):y(2)
                Region(i-y(1)+1,t-x(1)+1)=AYB05(i,t,k);
            end
    end
    xn=(x-(x(1)-1)*ones(3,1));
    yn=(y-(y(1)-1)*ones(3,1));
    xp=[xn(2);xn(1);xn(3)];
    yp=[yn(2);yn(1);yn(3)];
    contour(Region,30);
    hold on
    daspect([1,1,1])
    plot(xp,yp,'k')
    hold off
    pause(0.2);
end