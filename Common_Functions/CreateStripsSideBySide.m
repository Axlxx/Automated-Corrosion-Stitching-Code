function final = CreateStripsSideBySide(strips)
%CREATESTRIPSSIDEBYSIDE Summary of this function goes here
%   Detailed explanation goes here
    noStrips = size(strips,1);
    noRows = size(strips, 2);
    noCols = size(strips, 3);

    final = double(zeros(noRows, noStrips*noCols));

    for i = 1:noStrips
        final(:, (i-1)*noCols+1:i*noCols) = squeeze(strips(i,:,:));
    end
end

