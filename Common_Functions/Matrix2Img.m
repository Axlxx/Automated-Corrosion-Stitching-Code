function [img] = Matrix2Img(matrix, mins, maxs)
    if mins < 0
        mins = min(min(min(matrix)));
    end
    if maxs < 0
        maxs = max(max(max(matrix)));
    end

    % when grayscale - scale based on max and min
    img = (matrix - mins) ./ (maxs-mins);

    % remove values outside [0, 1)
    img(img < 0) = 0;
    img(img > 1) = 1;
end

