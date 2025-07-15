    
%%%%% Function for calculating difference between 2 pictures
%%% always takes the lowest by matching the left sides first
    function [threshold, histo] = CalculateError_percentile(ground_truth, panorama, intensityScale, min_value, verbose)
        
    %%%%% Overlap ground truth with panorama in the best way possible
    %%% Set up matching parameters
        FeatureDetector = @(X) detectKAZEFeatures(X);
        FeatureExtractor = @(X, Y) extractHOGFeatures(X, Y);
        FeatureMatcher = @(X, Y) matchFeatures(X, Y, ...
                'Unique', true, ...
                'Method', 'Exhaustive');    
        
    %%% Modify ground truth and pano sizes
        sizes = [max(size(ground_truth,1), size(panorama, 1)), max(size(ground_truth,2), size(panorama, 2))];
        new_ground_truth = ones(sizes);
        new_ground_truth(1:size(ground_truth,1), 1:size(ground_truth,2)) = ground_truth;
        ground_truth = new_ground_truth;
        new_panorama = ones(sizes);
        new_panorama(1:size(panorama,1), 1:size(panorama,2)) = panorama;
        panorama = new_panorama;
    
    %%% Blackout right parts of images
        blackout_ground_truth = BlackoutPartOfPicture(ground_truth, 0.6, 'black', "left") ;
        blackout_panorama = BlackoutPartOfPicture(panorama, 0.6, 'black', "left") ;
                
    %%% Get TM of matching panorama onto image
        [TM, PL, PR] = TransfMatrixOf2Strips(blackout_ground_truth, blackout_panorama, ...
                                            FeatureDetector, FeatureExtractor, ...
                                            FeatureMatcher, transltform2d);
        
    %%% Get panorama matched onto ground truth
        transfMatrices = [transltform2d, TM];
        images = zeros([2, sizes]);
        images(1, :, :) = 1;
        images(2, :, :) = panorama;
        panorama_modified = CreatePanorama_minOnly(transfMatrices,images); %% don't stitch the 2, just transform panorama based on TM
    
    %%% Graphs if verbose
        if verbose
        %%% Full panorama
            figure(21)
            subplot(1, 3, 1);
            imagesc(blackout_ground_truth*intensityScale + min_value);
            title('Ground Truth', 'FontSize',15)
            set(gca,'XColor', 'none','YColor','none')
            colorbar;
    
            subplot(1, 3, 2);
            imagesc(blackout_panorama*intensityScale + min_value);
            title('Reconstructed Panorama', 'FontSize',15)
            set(gca,'XColor', 'none','YColor','none')
            colorbar;
    
            subplot(1, 3, 3);
            imagesc(panorama_modified*intensityScale + min_value);
            title('Overlain panorama', 'FontSize',15)
            set(gca,'XColor', 'none','YColor','none')
            colorbar;
        end
    
    %%%%% Get difference between images
    %%% Calculate error and get difference image
        difference = double(imfuse(ground_truth, panorama_modified, "diff")) ./ 255;
        difference(difference == 1) = 0; % making sure

        bins = linspace(0,intensityScale, 101);
        histo = histcounts(difference, bins);
        if verbose
            figure(22);
            histogram(difference, bins);
            
            figure(23);
            ecdf(reshape(difference, 1, []));
            xlim([0, intensityScale]);
        end

        [f, x] = ecdf(reshape(difference, 1, []));
        index = min(find(f>0.95 == 1, 'FIRST'));
        threshold = x(index);

            %%% This is other method
            % errorRate = sum(difference, 'all');
            % SummedError = intensityScale * errorRate/area_ground_truth; % 5 mm for 1 pixel * average pixel intensity
    
    %%% Graphs if verbose
        if verbose
        %%% Ground truth + panorama + matched features + difference
            figure(24)
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
            sgtitle(sprintf("95%% of pixels lower than %0.3f error for colour range = %0.2f mm", threshold, intensityScale), 'FontSize',20)
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

