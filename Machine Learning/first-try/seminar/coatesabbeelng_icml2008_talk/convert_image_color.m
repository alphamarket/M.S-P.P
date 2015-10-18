I = imread('C:\users\dariush\desktop\Untitled.png');
%c = [1 1 1].*255;
c = [0.9333,0.9529,0.9765].* 255;
for i=1:size(I,1)
    for j=1:size(I,2)
        if mean(I(i,j)./255) >= 0.9
            I(i,j,:) = c;
        end
    end
end
imshow(I);
print -dmeta
disp('Copied to clipboard!');