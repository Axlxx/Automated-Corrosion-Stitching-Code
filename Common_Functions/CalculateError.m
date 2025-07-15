% Calculate error
function [SummedError] = CalculateError(ground_truth, panorama, intensityScale, min_value, verbose)
    % Take img 
    FeatureDetector = @(X) detectKAZEFeatures(X);
    FeatureExtractor = @(X, Y) extractFeatures(X, Y);
    FeatureMatcher = @(X, Y) matchFeatures(X, Y, ...
            'Unique', true, ...
            'Method', 'Exhaustive');    
    
    area_ground_truth = size(ground_truth, 1) * size(ground_truth, 2);

    %%% Modify sizes
    sizes = [max(size(ground_truth,1), size(panorama, 1)), max(size(ground_truth,2), size(panorama, 2))];
    new_ground_truth = ones(sizes);
    new_ground_truth(1:size(ground_truth,1), 1:size(ground_truth,2)) = ground_truth;
    ground_truth = new_ground_truth;

    new_panorama = ones(sizes);
    new_panorama(1:size(panorama,1), 1:size(panorama,2)) = panorama;
    panorama = new_panorama;
    
    %%% Get TM of matching panorama onto image
    [TM, PL, PR] = TransfMatrixOf2Strips(ground_truth, panorama, FeatureDetector, FeatureExtractor, FeatureMatcher, transltform2d);
    
    % Maybe
    % TM.Translation = round(TM.Translation);

    %%% Get panorama matched onto ground truth
    transfMatrices = [transltform2d, TM];
    images = zeros([2, sizes]);
    images(1, :, :) = 1;
    images(2, :, :) = panorama;

    % panorama_modified = CreatePanorama_minOnly(transfMatrices,images); %% don't stitch the 2, just transform panorama based on TM
    panorama_modified = panorama;
    %%% Calculate error and get difference image
    difference = double(imfuse(ground_truth, panorama_modified, "diff")) ./ 255;
    difference(difference == 1) = 0;
    errorRate = sum(difference, 'all');
    SummedError = intensityScale * errorRate/area_ground_truth; % 5 mm for 1 pixel * average pixel intensity
    
    %%% Calculate error metrics too
    % ssimval = ssim(panorama_modified , ground_truth);
    % peaksnr = psnr(panorama_modified ,ground_truth);
    % err = immse(panorama_modified ,ground_truth);

    if verbose
        figure(21)
        imagesc(panorama*intensityScale + min_value);
        title('Ground Truth', 'FontSize',15)
        set(gca,'XColor', 'none','YColor','none')
        colorbar;

        figure(22)
        subplot(2, 2, 1)
        imagesc(ground_truth*intensityScale + min_value);
        title('Ground Truth', 'FontSize',15)
        set(gca,'XColor', 'none','YColor','none')
        colorbar;
        % clim(gca, gca_range);
        
        subplot(2, 2, 2)
        imagesc(panorama_modified*intensityScale + min_value)
        title('Stitched Panorama', 'FontSize',15)
        set(gca,'XColor', 'none','YColor','none')
        colorbar;
        % clim(gca, gca_range);

        subplot(2, 2, 3)
        showMatchedFeatures(ground_truth, panorama_modified, PL, PR);
        set(gca, 'XColor', 'none', 'YColor', 'none')
        colorbar;
        title('Panorama matching', 'FontSize',15)
        
        subplot(2, 2, 4)
        imagesc(difference*intensityScale);
        set(gca, 'XColor', 'none', 'YColor', 'none')
        colorbar;
        % clim(0, 1);
        title('Absolute pixel difference', 'FontSize',15)
        sgtitle(sprintf("Error per pixel = %0.3f (mm) with colour range = %0.2f mm", SummedError, intensityScale), 'FontSize',20)
    end

    

end

function panorama = makePano(transfMatrices, reconstr, img)
    numImages = 2;
    imageSize = zeros(numImages,2);

% Compute Image Sizes & panorama sizes
    imageSize(1,:) = size(img);
    [xlim(1,:), ylim(1,:)] = outputLimits(transfMatrices(1), [1 imageSize(1,2)], [1 imageSize(1,1)]);
    imageSize(2,:) = size(reconstr);
    [xlim(2,:), ylim(2,:)] = outputLimits(transfMatrices(2), [1 imageSize(2,2)], [1 imageSize(2,1)]);
    

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


    panorama = ones([height width 1], 'like', reconstr); % use with min
    panoramaView = imref2d([height width], xLimits, yLimits);
    
    for i = 1:numImages
        % If tforms is I then break;
        if(isequal(transfMatrices(i), transfMatrices(1))  && i ~= 1)
            break;
        end

        % Compute Warped Image
        if i == 1
            I = img;
        else 
            I = reconstr;
        end
       
        % Transform I into the panorama.
        warpedImage = imwarp(I, transfMatrices(i), 'OutputView', panoramaView);
                      
        % Generate a binary mask.    
        mask = imwarp(true(size(I,1),size(I,2)), transfMatrices(i), 'OutputView', panoramaView);
        
        % Overlay the warpedImage onto the panorama.
        % Choose Min pixel
        panorama  = imblend(warpedImage,panorama,mask,mode="Min");

        % imshow(panorama);
    end

end