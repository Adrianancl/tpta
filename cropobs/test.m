t = imread('digitalizar0080.jpg');
B = edge(t(:,:,1),'sobel');
figure
imshow(B);
v = sum(B);
l = length(v);
ini = [];
for i = 1:l
    if v(i) > 70
        if isempty(ini)
            ini = i;
        end
        fin = i;
    end
end
px = [ini fin];

v = sum(B,2);
l = length(v);
ini = [];
count = 0;
for i = 1:l
    if v(i) > 70
        if isempty(ini)
            ini = i;
        end
        fin = i;
    end
end
py = [ini fin];

I = imcrop(t,[px(1) py(1) px(2)-px(1) py(2)-py(1)]);
figure
image(I)
daspect([1,1,1])
% 
% [H, THETA, RHO] = hough(B);%Hough transformation
% P  = houghpeaks(H, 10);
% lines = houghlines(B, THETA, RHO, P, 'FillGap', 15,'MinLength',50);
% for i=1:l
%         Theta=lines(i).theta;
%         t1 = lines(i).point1;
%         t2 = lines(i).point2;
%         imshow(B)
%         hold on
%         plot([t1(1) t2(1)]',[t1(2) t2(2)]','xm')
%         hold off
% end
