t = imread('digitalizar0080.jpg');
B = edge(t(:,:,1),'sobel');
[m,n] = size(B);
m2 = round(m/2);
n2 = round(n/2);
C1 = imcrop(B,[1 1 n2 m2]);
C2 = imcrop(B,[n2+1 1 n-n2 m2]);
C3 = imcrop(B,[1 m2+1 n2 m-m2]);
C4 = imcrop(B,[n2+1 m2+1 n-n2 m-m2]);
imshow(C1)
figure
imshow(C2)
figure
imshow(C3)
figure
imshow(C4)
[H, THETA, RHO] = hough(C1);%Hough transformation
P  = houghpeaks(H, 10);
lines = houghlines(C1, THETA, RHO, P, 'FillGap', 15,'MinLength',50);
l = length(lines);
for i=1:l
        Theta=lines(i).theta;
        t1 = lines(i).point1;
        t2 = lines(i).point2;
        imshow(C1)
        hold on
        plot([t1(1) t2(1)]',[t1(2) t2(2)]','xm')
        hold off
end