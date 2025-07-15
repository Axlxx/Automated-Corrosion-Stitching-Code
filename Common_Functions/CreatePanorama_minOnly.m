function panorama_min = CreatePanorama_minOnly(transfMatrices,strips)
    numImages = numel(transfMatrices);
    imageSize = zeros(numImages,2);

    % strips(strips == 0) = 0.00001;
    % strips(strips == 1) = 0.99999;
    
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
    panorama_min = ones([height width 1], 'like', strips(1,:,:)) ; % use with min
    panoramaView = imref2d([height width], xLimits, yLimits);
    
    % figure(99);
    for i = 1:numImages
        % If tforms is I then break;
        if(isequal(transfMatrices(i), transfMatrices(1))  && i ~= 1)
            break;
        end

        % Compute Warped Image
        I = squeeze(strips(i, :, :)); 
       
        % Transform I into the panorama.
        warpedImage = imwarp(I, transfMatrices(i), 'OutputView', panoramaView);
        % warpedImage(warpedImage == 0) = 1;
        % imshow(warpedImage)
                      
        % Generate a binary mask.    
        mask = imwarp(true(size(I,1),size(I,2)), transfMatrices(i), 'OutputView', panoramaView);
        
        % Overlay the warpedImage onto the panorama.
        % Choose Min pixel
        panorama_min  = imblend(warpedImage,panorama_min,mask,mode="Min");
        % imshow(panorama_min);

        mask;
        % imshow(panorama);
    end
    % panorama_min(panorama_min == 1) = 0;

    % strips(strips == 0.00001) = 0;
    % strips(strips == 0.99999) = 1;


end