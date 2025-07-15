function ShowStrips(strips)
    %SHOWSTRIPS Show strips

    min_depth = min(strips, [],  'all');
    max_depth = max(strips, [],  'all');
    gca_range = [min_depth, max_depth];
    % gca_range = [0, 11];


    for i = 1:size(strips, 1)
        subplot(1, size(strips, 1), i)

        imagesc(squeeze(strips(i, :, end:-1:1))) %ImgBack2OutputMatrix(img))
        set(gca,'XColor', 'none','YColor','none')
        axis image
        colorbar;
        clim(gca, gca_range);
        title(sprintf("Strip %d", i))


    end
end

