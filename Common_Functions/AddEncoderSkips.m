function strips = AddEncoderSkips(strips, encoder_skips)
%ADDENCODERSKIPS Summary of this function goes here
%   Detailed explanation goes here
    noStrips = size(strips, 1);
    noRows = size(strips, 2);
    noColumns = size(strips, 3);

    verbose = 0;

    for i =  1:noStrips
        %%% Generate skipped rows
        
        
        skips = randperm(noRows,encoder_skips);
        how_many = int32(rand() * (encoder_skips+1)) - 1;
        if how_many < 0
            how_many = 0;
        end

        skips = skips(1:how_many);
        skips = sort(skips,2, "ascend");

        %%% Remove rows
        for skip_index = skips
            strips(i, skip_index,:) = zeros(1, noColumns);
        end
        if verbose
            figure(50);
            imshow(squeeze(strips(i,:,:)))
        end
        

        %%% Reorganize rows
        for j = how_many:-1:1
            strips(i, skips(j):end-1, :) = strips(i, skips(j)+1:end, :);
        end
        if verbose
            imshow(squeeze(strips(i,:,:)))
        end

        %%% Remove last rows
        strips(i, end-how_many+1:end, :) = ones(how_many, noColumns) * max(strips(i, :, :), [], 'all');
        if verbose
            imshow(squeeze(strips(i,:,:)))
        end


    end

end

