function [newImg] = PadPictureWithZeros(img, padSize, color)
%PADPICTUREWITHZEROS Summary of this function goes here
%   Detailed explanation goes here
    if size(padSize, 2) == 1
        padX = padSize;
        padY = padSize;
    else 
        padX = padSize(1);
        padY = padSize(2);
    end

    newImg = ones( size(img, 1) + 2 * padX, ...
                    size(img, 2) + 2 * padY);

    if isequal(color, 'black')
        newImg = newImg * 0;
    end

    newImg(padX+1:end-padX, padY+1:end-padY) = img;

end

