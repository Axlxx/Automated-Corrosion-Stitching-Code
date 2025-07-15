%% Simulation program to find the performance of each Image Registration algorithm
% Will load a data set
% Divide each map into two strips and add inspection errors
% Try to stitch them back using the specified algorithms
% Outputs the results for each method 
%
%
% Author: Alexandru Nichita
% Lab: Non-Destructive Evaluation, Department of Mechaincal Engineering 
% University: Imperial College London
% Last Edit Date: July 2025


% Flag to not clear previous data (e.g. for comparison)
newRun = 1;
if newRun
    clear;
    newRun = 1;
    ClearFigures();
end

%%% Get pictures
% Will get all depth maps from the folder specified that has the given
% token inside. 
% names = GetAllWildcardNamesFromFolder("fodler_name\", 'part_of_name', 'format', number_of_maps);
names = sort(GetAllWildcardNamesFromFolder("DatabaseCombined\", 'np_', '.mat', -1), 'ascend')';

%%% Range in which to look
irange = [1:1:length(names)]; % Which depth maps to include in dataset
jrange = [1, 3, 5:5:50];    % Which overlap levels to analyse

%%% Detector parameters and Save locations 
%%% Flags 
% Type of noise used
noiseTypes = ["WHITE", "GAUSSIAN"]; 

% Type of feature detectors used, to add more, modify the getDetExtr() 
% function at the bottom of the page and add a short name for the detector here
detectors = ["KAZEHoG", "KAZE", "SIFT", "SURF", "SIFTHoG", "SURFHoG"];
dets_used = [1, 2, 3, 4, 5, 6]; % Use all detectors above or just some (remove as needed)

% Name of current run
folder_suffix = "NewNoise_Oneshots";
folders = [];

%%% Create folders
if newRun
    for index_det = 1:length(detectors)
        % check the "Results" folder for the runs that have the "folder_suffix" in their 
        % name that you have defined above to find the results of the latest run
        folder_name = strcat("Runs_", detectors(index_det), "_" ,folder_suffix);
        folders = [folders, folder_name];
    
        if ~exist("folder_name", "dir")
            mkdir(folder_name);
            addpath(folder_name);
        end
    
        save(sprintf('%s/General.mat',folder_name), "irange", "jrange");
    end
end

%%% Noise levels
% Selects thresholds for the errors tested 
% To add a new level simply add another column to each of the vectors below
noise_intensities = [0, 0.1, 0.5, 1, 2]; % Extend for noise values
encoder_skips_arr = [0, 5, 10, 15, 20];  % Extend for % of encoder skip
degs =              [0, 5, 10, 15, 20];  % Extend for maximum rotation (in degrees) of second strip

% You can also select a subset of thresholds (e.g. [1, 3, 4]) 
% or add your new thresholds at the end [1:6]
threshold_levels = [1, 2, 3, 4, 5];    

%%% Simulator parameters
strip_size = 100;
defaultRIGID = rigidtform2d();

%%% Results parammeters
if newRun
    %%%                    det              , intensity,              , i             , j
    % If stitch was a success for each detector/ each noise threshold/ each depth map/ each overlapping level
    success =       zeros(length(detectors), length(noise_intensities), length(irange), 50);

    % Associated error when comparing stitched panorama to ground truth
    error =         zeros(length(detectors), length(noise_intensities), length(irange), 50);

    % What SSIM thought was the overlap and the associated value within [-1, 1]
    perceived_j =   zeros(length(noise_intensities), length(irange), 50);
    metric_ssim =   zeros(length(noise_intensities), length(irange), 50);

    % Depth range of each map 
    intensity =     zeros(length(irange));
end

%%% verbose flag
verbose = true;

%%% Load next image 
tic; % timing
for i=irange % For each depth map
    fprintf("\n\n            Loading new image sample no =       %d\n\n", i);
    loader = load(names(i));
    output_matrix = loader.output_matrix;

    %%% Get intensity range
    min_depth = min(output_matrix, [],  'all');
    max_depth = max(output_matrix, [],  'all');
    gca_range = [min_depth, max_depth];
    intensity(i) = max_depth - min_depth;

    % Useful function for later
    ImgBack2OutputMatrix = @(img) img*(max_depth-min_depth) + min_depth;
    
    %%% Show ground truth
    if verbose
        figure(1)
        imagesc(output_matrix) %ImgBack2OutputMatrix(img))
        set(gca,'XColor', 'none','YColor','none')
        axis image
        colorbar;
        clim(gca, gca_range);
        title('Simulated Depth Corrosion Map (cm)', 'FontSize',20);
    end

    %%% Set error intensity
    for threshold = threshold_levels % for each noise threshold
        fprintf("\nAnalysing noise threshold %d\n", threshold);

        noise_intensity = noise_intensities(threshold);
        encoder_skips = encoder_skips_arr(threshold);
        maxDeg = degs(threshold);
    
        %%% Set overlap percentage
        for j = jrange % For each overlapping size

            percentage = double(j/100);
            ground_truth = output_matrix(:, 1:int32(floor(double(strip_size*(2-percentage)))));
            fprintf("Analysing percentage = %d\n", j);
            
            %%% Divide into strips with overlap percentage
            [strips_noisefree_more, step] = DivideMatrixIntoStrips(output_matrix, strip_size, percentage);
            strips_noisefree = strips_noisefree_more(1:2,:,:);

            %%% Modify strips (add inspection errors)
            % Add Noise 
            noiseType = noiseTypes(1);
            strips = AddNoiseToMatrix(strips_noisefree,noise_intensity, noiseType);
            % Add Encoder Skips
            strips = AddEncoderSkips(strips, encoder_skips);
            % Add Rotation
            angle = normrnd(maxDeg, maxDeg/3*1) - maxDeg; % idk what this wizardry i wrote is either
            strips = AddRotations(strips,[0, angle], 2);
            strips_angled = AddRotations(strips_noisefree,[0, angle], 2);
            strips_angled(1,:,:) = Matrix2Img(squeeze(strips_angled(1, :, :)), -1, -1);
            strips_angled(2,:,:) = Matrix2Img(squeeze(strips_angled(2, :, :)), -1, -1);

           
            %%% SSIM find largest overlap
            jSSIMvals = zeros(1, 50);
            for t = 5:50
                jSSIMvals(t) = ssim(squeeze(strips(1, :, end-t+1:end)), squeeze(strips(2, :, 1:t)));
            end

            [SSIM_val, j_used] = max(jSSIMvals(:), [], "all");
            perceived_j(threshold, i, j) = j_used;
            metric_ssim(threshold, i, j) = SSIM_val;

            % img = Matrix2Img(output_matrix, -1, -1); % 5, 10
            %%% Strip As-is
            StripLeft =     Matrix2Img(squeeze(strips(1, :, :)), -1, -1);
            StripRight =    Matrix2Img(squeeze(strips(2, :, :)), -1, -1);

            if verbose
                % Show the 2 strips before transformations
                figure(2)
                subplot(1,2, 1)
                imagesc(ImgBack2OutputMatrix(StripLeft)) %, [L, R]);
                title('Example Strip 1', 'FontSize',20);
                set(gca,'XColor', 'none','YColor','none')
                axis image
                colorbar;
                clim(gca, gca_range);
                subplot(1,2, 2)
                imagesc(ImgBack2OutputMatrix(StripRight))
                axis image
                title('Example Strip 2', 'FontSize',20);
                set(gca,'XColor', 'none','YColor','none')
                colorbar;
                clim(gca, gca_range);

                figure(7)
                subplot(1,2, 1)
                imagesc(squeeze(strips_angled(1, :, :))); %, [L, R]);
                title('Example Strip 1', 'FontSize',20);
                set(gca,'XColor', 'none','YColor','none')
                axis image
                colorbar;
                clim(gca, gca_range);
                subplot(1,2, 2)
                imagesc(squeeze(strips_angled(2, :, :)));
                axis image
                title('Example Strip 2', 'FontSize',20);
                set(gca,'XColor', 'none','YColor','none')
                colorbar;
                clim(gca, gca_range);
            end

            %%% Prepare Strip
            percentage = double(j_used/100);
            StripLeft_padded =     prepareStrip(StripLeft, percentage, "right") ;
            StripRight_padded =    prepareStrip(StripRight, percentage, "left") ;

            %%% Check for all detector types
            for det = dets_used %5:length(detectors)
                
                %%% Set detector
                detector = detectors(det);
                [FeatureDetector, FeatureExtractor] = getDetExtr(detector);

                %%% Prepare for matrix loop
                % TransfMatrix = defaultRIGID;
                ratio = 0.05;
                found_flag = 0;
                
                %%% Do this while you haven't found enough matches or
                %%% the matches found give a bad solution
                while(~found_flag)
                    %%% Update matcher
                    FeatureMatcher = @(X, Y) matchFeatures(X, Y, 'Unique', true, 'Method', 'Approximate', 'MaxRatio', ratio);    
                        
                    %%% Compute transformation
                    [TransfMatrix, matchedPointsLeft, matchedPointsRight, inlierIndex] = TransfMatrixOf2Strips( ...
                        StripLeft_padded, StripRight_padded, ...
                        FeatureDetector, FeatureExtractor, ...
                        FeatureMatcher, defaultRIGID );

                    if verbose
                        figure(3)
                        subplot(1, 3, 1)
                        imagesc(ImgBack2OutputMatrix(StripLeft)) %, [L, R]);
                        axis image
                        features = FeatureDetector(StripLeft);
                        hold on;
                        plot(features.selectStrongest(30))
                        hold off;
                        title('Example Strip 1', 'FontSize',10);
                        set(gca,'XColor', 'none','YColor','none')
                        colorbar;
                        clim(gca, gca_range);
                        
                        subplot(1, 3, 2)
                        imagesc(ImgBack2OutputMatrix(StripRight))
                        axis image
                        features = FeatureDetector(StripRight);
                        hold on;
                        plot(features.selectStrongest(30))
                        hold off;
                        title('Example Strip 2', 'FontSize',10);
                        set(gca,'XColor', 'none','YColor','none')
                        colorbar;
                        clim(gca, gca_range);

                        subplot(1, 3, 3)
                        showMatchedFeatures(    StripLeft, StripRight, ...
                                                matchedPointsLeft, matchedPointsRight, ...
                                                PlotOptions={'bo', 'ro', 'y--'});
                        

                        if (sum(inlierIndex) > 0)
                            hold on;
                            % linespecs = {'y-', 'LineWidth', 3};
                            % showMatchedFeatures(    StripLeft, StripRight, ...
                            %                         matchedPointsLeft(inlierIndex, :) , matchedPointsRight(inlierIndex, :), ...
                            %                         PlotOptions={'bo', 'ro', 'y-'});
                            pointsShowLeft = parsePoints(matchedPointsLeft(inlierIndex, :), 1);
                            pointsShowRight = parsePoints(matchedPointsRight(inlierIndex, :), 2);
                            lineX = [pointsShowLeft(:,1)'; pointsShowRight(:, 1)'];
                            numPts = numel(lineX);
                            lineX = [lineX; NaN(1,numPts/2)];
                            
                            lineY = [pointsShowLeft(:, 2)'; pointsShowRight(:, 2)'];
                            lineY = [lineY; NaN(1,numPts/2)];
                            
                            plot(lineX(:), lineY(:), 'y-', 'LineWidth', 2); % line
                            hold off;
                        end
                        
                        title(sprintf("Identified Matches = %d", length(matchedPointsRight)), 'FontSize', 10)
                        
                    end

                    %%% Check for success
                    if (sum(inlierIndex) >= 2 && ... %%% At least 2 matches found
                        abs(TransfMatrix.RotationAngle - angle) < 15) %%% Angle not too large
                        found_flag = 1;
                        success(det, threshold, i, j)=1;
                        sgtitle('Succeded to stitch', 'FontSize', 15)
                        break; % If it found the neccesary matches, it exits the loop
                    end
    
                    %%% Relax required match strength (loop condition) each failed loop
                    ratio = ratio + 0.2;

                    %%% Check for fail
                    if ratio > 1
                        success(det, threshold, i, j)=0;
                        fprintf("Not found for detector = %s, angle = %d, matched points = %d\n", detector, TransfMatrix.RotationAngle, length(matchedPointsLeft));
                        sgtitle('Failed to stitch', 'FontSize', 15)
                        
                        break;
                    end
                    
                end % End while
            
                %%% It either found something or nothing
                %%% Generate panorama
                % TransfMatrix.Translation = round(TransfMatrix.Translation);
                % If nothing was found the second strip is not stitched so
                % it will return an image with only the first strip
                reconstruction_min = CreatePanorama_minOnly([defaultRIGID, TransfMatrix],strips_angled);

                %%% Get error metric
                % This is for analysis using histograms
                [thr, histo] =  CalculateError_percentile(Matrix2Img(ground_truth, -1, -1), reconstruction_min, intensity(i), min_depth, verbose);
                % This was used in the paper
                [SummedError] = CalculateError(Matrix2Img(ground_truth, -1, -1), reconstruction_min, intensity(i), min_depth, verbose);

                error(det, threshold, i, j) = SummedError;

                if(verbose)
                    % Show the Image Reconstruction compared to ground truth
                    figure(4)
                    sgtitle('Image Reconstruction compared to ground truth', 'FontSize', 20);
        
                    subplot(2, 2, 1)
                    hold on;
                    imagesc(ground_truth)
                    plot(ones(size(ground_truth,2))*60, 'k-', 'LineWidth',2)
                    hold off;
                    set(gca,'XColor', 'none','YColor','none')
                    axis image
                    colorbar;
                    title('Simulated corrosion map', 'FontSize', 15)
                    % clim(gca, [5,10]);
        
                    subplot(2, 2, 2)
                    hold on;
                    imagesc(ImgBack2OutputMatrix(reconstruction_min))
                    plot(ones(size(reconstruction_min,2))*60, 'r-','LineWidth',2)
                    hold off;
                    set(gca,'XColor', 'none','YColor','none')
                    axis image
                    colorbar;
                    title('Reconstructed Panorama', 'FontSize', 15)
                    % clim(gca, [5,10]);
                    
                    subplot(2, 2, 3:4)
                    hold on;
                    plot(ground_truth(60,1:end), 'k-', 'LineWidth',2);
                    plot(ImgBack2OutputMatrix(reconstruction_min(60,:)), 'r-','LineWidth',2);
                    hold off;
                    ylabel('Depth (cm)');
                    title('Cross-section comparison of ground truth and reconstruction',  'FontSize', 15)
                    legend('Ground Truth', 'Reconstructed Panorama', 'Location','southeast',  'FontSize', 15);
                end
          
            end
        end
    end
end

stop = toc % Show time taken for everything to run

%%% Save results
%    det              , intensity,              , i             ,j
% success error metric_ssim intensity perceived_j 
for det = dets_used %5:length(detectors)
    for threshold  = 1:length(noise_intensities)
        
        success_det = squeeze(success(det, threshold, :, :));
        error_det = squeeze(error(det, threshold, :, :));
        metric_ssim_det = squeeze(metric_ssim(threshold, :, :));
        perceived_j_det = squeeze(perceived_j(threshold, :, :));
        intensity_det = intensity(:);

        save(sprintf('%s/Detector_%s__Type_%s__NoiseLvl_%d.mat', ...
            folders(det), detectors(det), "WHITE", noise_intensities(threshold)), ...
            'success_det', 'error_det', ...
            'metric_ssim_det', 'perceived_j_det', ...
            'intensity_det');
    end
end


% Done
figure(1);
title("Done");


%% Functions

% Selects the correct detector and extraction method at the start based on
% the name
% Mofify this if you want to add more feature detectors
function [FeatureDetector, FeatureExtractor] = getDetExtr(detector)
    if isequal(detector,'KAZE') | isequal(detector,'KAZEHoG')
        FeatureDetector = @(X) detectKAZEFeatures(X);
    elseif isequal(detector,'SIFT') | isequal(detector,'SIFTHoG')
        FeatureDetector = @(X) detectSIFTFeatures(X);
    elseif isequal(detector,'ORB')
        FeatureDetector = @(X) detectORBFeatures(X);
    elseif isequal(detector,'FAST')
        FeatureDetector = @(X) detectFASTFeatures(X);
    elseif isequal(detector,'SURF') | isequal(detector,'SURFHoG')
        FeatureDetector = @(X) detectSURFFeatures(X);
  % elseif isequal(detector,'HarrisDet') | isequal(detector,'anythingAtAll') % Example of addition
  %     FeatureDetector = @(X) detectHarrisFeatures(X);
    end
        
    % or more extraction methods
    if isequal(detector,'KAZEHoG') | isequal(detector,'SIFTHoG') | isequal(detector,'SURFHoG')
        FeatureExtractor = @(X, Y) extractHOGFeatures(X, Y);
    else
        FeatureExtractor = @(X, Y) extractFeatures(X, Y);
    end
end

% Function that helps blackout the area not identified by SSIM and
% pads the strips with blank content so HOG does not reject features just
% because it has no space to compute the values at the edges
function modifiedStrip = prepareStrip(strip, percentage, orientation)
    %%% Remove non-data part of strip
    modifiedStrip = RemoveNonDataStrip(strip);

    %%% Remove non overlapping part
    modifiedStrip = BlackoutPartOfPicture(modifiedStrip, 1 - percentage, 'black', orientation) ;
            
    %%% Add padding to picture 
    modifiedStrip = PadPictureWithZeros(modifiedStrip, [0, 10], 'black');
end

%%% Thing to get point locations
% Just for visualisation, don't touch this 
% Why would you ever need to touch this
% Just don't 
function points=parsePoints(points, ptsInputNumber)
    
    fcnInputVarNumber = 2 + ptsInputNumber; 
    varName = ['matchedPoints', num2str(ptsInputNumber)];
    
    if ~isa(points, 'vision.internal.FeaturePoints') && ~isa(points, 'MSERRegions')
        validateattributes(points,{'int16', 'uint16', 'int32', 'uint32', ...
            'single', 'double'}, {'2d', 'nonsparse', 'real', 'size', [NaN 2]},...
            mfilename, varName, fcnInputVarNumber);
    else
        points = points.Location;
    end
    
    points = double(points);
end