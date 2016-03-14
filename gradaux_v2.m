function [tx,ty] = gradaux_v2(frame,step)
[m,n] = size(frame);%m=line and n=columns
gradx = zeros(size(frame));
grady = zeros(size(frame));
%--------------------------------------------------------------------------
%gradient in x direction

for j = 1:n
    if abs(j - 1) < step || abs(j - n) < step
        if abs(j - 1) < step
            if abs(j + step) < n
                gradx(:,j) = (frame(:,j+step) - frame(:,1))/(2*step);
            else
                gradx(:,j) = (frame(:,n) - frame(:,1))/(2*step);
            end
        else
            gradx(:,j) = (frame(:,n) - frame(:,j-step))/(2*step);
        end
    else
        gradx(:,j) = (frame(:,j+step) - frame(:,j-step))/(2*step);
    end
end

%--------------------------------------------------------------------------
%gradient in y direction

for i = 1:m
    if abs(i - 1) < step || abs(i - m) < step
        if abs(i - 1) < step
            if abs(i + step) < m
                grady(i,:) = (frame(i+step,:) - frame(1,:))/(2*step);
            else
                grady(i,:) = (frame(m,:) - frame(1,:))/(2*step);
            end
        else
            grady(i,:) = (frame(m,:) - frame(i-step,:))/(2*step);
        end
    else
        grady(i,:) = (frame(i+step,:) - frame(i-step,:))/(2*step);
    end
end

%--------------------------------------------------------------------------
tx = gradx;
ty = grady;
end