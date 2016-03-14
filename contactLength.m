function l = contactLength(Frames,b,e)
    v = zeros(e-b+1,1);
    for i = b:e
        imagesc(Frames(i).f)
        imdistline
        v(i-b+1) = input('What is the value of the contact length for this frame?');
        close all
    end
    boxplot(v)
    l = mean(v);
end