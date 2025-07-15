function [strip_cropped] = RemoveNonDataStrip(strip)
    % For top and bottom
    rows = size(strip, 1);
    black_row_top = 0;
    black_row_bottom = 0;
    for r = 1:rows
        % If more than half of values are 0
        if nnz(strip(r,:)==0) >= size(strip, 2)/2
            black_row_bottom = black_row_bottom + 1;
        else 
            break;
        end
    end

    for r = rows:-1:1
        if sum(strip(r,:)) == 0
            black_row_top = black_row_top + 1;
        else 
            break;
        end
    end

    strip_cropped = zeros(size(strip,1) - black_row_bottom - black_row_top, size(strip, 2));
    strip_cropped(:,:) = strip(black_row_bottom+1:end-black_row_top, :);



end

