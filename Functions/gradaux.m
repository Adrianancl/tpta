function [tx,ty] = gradaux(frame,step)
[m,n] = size(frame);%m=line and n=columns
gradx = zeros(size(frame));
grady = zeros(size(frame));
%--------------------------------------------------------------------------
%gradient in x direction
for i = 1:m
    for j = 1:n
        if abs(j - 1) < step || abs(j - n) < step
            if abs(j - 1) < step
                if abs(j + step) < n
                    gradx(i,j) = (frame(i,j+step) - frame(i,1))/(2*step);
                else
                    gradx(i,j) = (frame(i,n) - frame(i,1))/(2*step);
                end
            else
                gradx(i,j) = (frame(i,n) - frame(i,j-step))/(2*step);
            end
        else
            gradx(i,j) = (frame(i,j+step) - frame(i,j-step))/(2*step);
        end
    end
end
%--------------------------------------------------------------------------
%gradient in y direction
for j = 1:n
    for i = 1:m
        if abs(i - 1) < step || abs(i - m) < step
            if abs(i - 1) < step
                if abs(i + step) < m
                    grady(i,j) = (frame(i+step,j) - frame(1,j))/(2*step);
                else
                    grady(i,j) = (frame(m,j) - frame(1,j))/(2*step);
                end
            else
                grady(i,j) = (frame(m,j) - frame(i-step,j))/(2*step);
            end
        else
            grady(i,j) = (frame(i+step,j) - frame(i-step,j))/(2*step);
        end
    end
end
%--------------------------------------------------------------------------
tx = gradx;
ty = grady;
end