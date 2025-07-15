function panoramas = CreatePanorama(transfMatrices,strips)
    numImages = numel(transfMatrices);
    imageSize = zeros(numImages,2);

% Compute Tranf Matrix
    for n = 2:numImages
        transfMatrices(n).A = transfMatrices(n-1).A * transfMatrices(n).A; 
    end

% Compute Image Sizes & panorama sizes
    for n = 1:numImages
        imageSize(n,:) = size(squeeze(strips(n,:,:)));
        [xlim(n,:), ylim(n,:)] = outputLimits(transfMatrices(n), [1 imageSize(n,2)], [1 imageSize(n,1)]);
    end
    maxImageSize = max(imageSize);
    if numImages == 1
        xMin = min([1; xlim(:)]); xMax = max([maxImageSize(1); xlim(:)]);
        yMin = min([1; ylim(:)]); yMax = max([maxImageSize(1); ylim(:)]);
    else 
        xMin = min([1; xlim(:)]); xMax = max([maxImageSize(2); xlim(:)]);
        yMin = min([1; ylim(:)]); yMax = max([maxImageSize(1); ylim(:)]);
    end
    width  = round(xMax - xMin);
    height = round(yMax - yMin);
    xLimits = [xMin xMax];
    yLimits = [yMin yMax];
    
% Create the panorama.
    % blender = vision.AlphaBlender('Operation', 'Binary mask', ...
    %     'MaskSource', 'Input port');  
    % blender = vision.AlphaBlender('Operation', 'Blend', ...
    %     'MaskSource', 'Input port'); 
    panoramaMax = zeros([height width 1], 'like', strips(1,:,:)); % use with max
    panorama = ones([height width 1], 'like', strips(1,:,:)); % use with min
    panoramas = ones([10, size(panorama)]);
    panoramas(4, :, :) = panoramaMax;
    panoramaView = imref2d([height width], xLimits, yLimits);
    
    for i = 1:numImages
        % If tforms is I then break;
        if(isequal(transfMatrices(i), transfMatrices(1))  && i ~= 1)
            break;
        end

        % Compute Warped Image
        if length(size(strips)) == 3
            I = squeeze(strips(i, :, :)); 
        else 
            I = squeeze(strips(:, :)); 
        end
       
        % Transform I into the panorama.
        warpedImage = imwarp(I, transfMatrices(i), 'OutputView', panoramaView);
                      
        % Generate a binary mask.    
        mask = imwarp(true(size(I,1),size(I,2)), transfMatrices(i), 'OutputView', panoramaView);
        
        % Overlay the warpedImage onto the panorama.
        % Old
        % panorama = step(blender, panorama, warpedImage, mask);

        % show differences in strips
        panoramas(1, :, :) = imblend(warpedImage,squeeze(panoramas(1, :, :)),mask,mode="Guided");

        % show regions of overlap
        panoramas(2, :, :)  = imblend(warpedImage,squeeze(panoramas(2, :, :)),mask,mode="Alpha");

        % Choose Average pixel (shadows rest of picture)
        panoramas(3, :, :)  = imblend(warpedImage,squeeze(panoramas(3, :, :)),mask,mode="Average");

        % Choose Max pixel
        panoramas(4, :, :)  = imblend(warpedImage,squeeze(panoramas(4, :, :)),mask,mode="Max");

        % Choose Min pixel
        panoramas(5, :, :)  = imblend(warpedImage,squeeze(panoramas(5, :, :)),mask,mode="Min");

        % imshow(panorama);
    end
end