function strip = BlackoutPartOfPicture(strip, removePercentage, color, flag)
    n_x = size(strip, 1);
    n_y = size(strip, 2);
    edge_size = floor(n_y*removePercentage);

    col = 0;
    if isequal(color, 'white')
        col = 1;
    end

    if flag == "left"
        strip(1:end, n_y-edge_size+1:end) = col;
        % strip = strip(:, 1:edge_size);
    else 
        strip(1:end, 1:edge_size) = col;
        % strip = strip(:, n_y-edge_size:n_y);
    end
end

