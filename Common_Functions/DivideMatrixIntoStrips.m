% Produce columns of width = ***width*** that have an overlap = ***overlap*** between them
function [strips, step] = DivideMatrixIntoStrips(Matrix, width, overlap)
    % Length of matrix
    rows = length(Matrix(:,1,1));
    cols = length(Matrix(1,:,1));

    % Step size
    step = (width * (1-overlap));

    % Number of strips
    nostrips = double(cols / step);
    % nostrips = int32(floor(nostrips))+1;
    nostrips = int32(floor(nostrips));
    strips = zeros(nostrips, rows, width);

    for i = 1:nostrips
        init_pos = floor(double((i-1))*step);
        if (init_pos+width > rows)
            % copyS = strips(1:i-1, :, :);
            % % delete strips;
            % % strips = zeros(nostrips-1, cols, width);
            % strips = copyS(1:end-1, :, :);
            strips(end, :, :) = Matrix(:, end-width+1:end);
            break;
        end
        strips(i, :, :) = Matrix(:, init_pos+1:init_pos+width);
    end
    
end

