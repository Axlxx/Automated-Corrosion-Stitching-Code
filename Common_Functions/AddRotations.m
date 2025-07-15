function strips = AddRotations(strips,angles, method)
    noStrips = size(strips, 1);

    for i = 2:noStrips
        if method == 2
            max_str = max(strips(i, :,:), [], "all");

            S = squeeze(max_str-strips(i, :, :));
            S = imrotate(S, angles(i), 'bicubic', 'crop');
            strips(i, :, :) = max_str - S;
        else % method == 1
            strips(i, :, :) = imrotate(squeeze(strips(i, :, :)), angles(i), 'bicubic', 'crop');
        end
    end
end

