%% Program used to stitch data obtained in the NDE lab using a Verasonics probe 
% (more details about acquistion method in the paper)
% Will load a data set of A-scans
% Will convert to depth maps
% Takes the first two depth maps (can be extended for more with a for loop)
% Tries to stitch them back using the specified algorithm
% Outputs the result and the workings of the algorithm (especially if it
% cannot find enough matches)
%
%
% Author: Alexandru Nichita
% Lab: Non-Destructive Evaluation, Department of Mechaincal Engineering 
% University: Imperial College London
% Last Edit Date: July 2025

ClearFigures();

%%% General parameters
detectors = ["KAZEHoG", "KAZE", "SIFT", "SURF", "SIFTHoG", "SURFHoG"];
defaultRIGID = rigidtform2d();
angle = 0;

%%% Run parameters
% Which detector and extractor to use to try and stitch. 
dets_used = [1];


% Select aperture of PAUT probe
aperture = 6;
if ~exist("last_aperture", "var")
    last_aperture = aperture+1;
end

%%% verbose flag
verbose_signal = false; % If you want to see how the depth maps are generated
verbose = true; % How the image registration works (and everything else)

%%% Results parammeters
success =       0; % For debug / looking in workspace or console

%%% get strips image 
if ~exist("strips", "var") || last_aperture ~= aperture
    tic;
    % Set speed of sound based on your calibration 
    % I preloaded the ones used for the datasets in the paper here
    % Uncomment the file loading and speed you are using
    load("Steel_Database_Steel_2.mat");
    speed = 5.9809e+06; 

    % load("Steel_Database_Aluminium_3.mat");
    % speed = 5.7003e+07;

    % Open or replace to make any modifications to A-scan processing
    strips = stripsFromRealData(steel_db, aperture,  speed, 62.5 * 10^6, verbose_signal);
    toc;
end

% If running the program multiple times with the same apertur, there is no
% need to redo the depth maps
% If you switch the loading file, clear the entire workspace 
last_aperture = aperture;

%%% Get intensity range
% For display in figures
min_depth = 7.5; %min(strips, [],  'all');
max_depth = 10.5; %max(strips, [],  'all');
gca_range = [min_depth, max_depth];
intensity = max_depth - min_depth; % not used anywhere

% Useful function for later
ImgBack2OutputMatrix = @(img) img*(max_depth-min_depth) + min_depth;

% if you want to find the ECDF of each data set
if (false)
    % ecdfs_of_apertures_S1 = zeros(4000, 1);
    % ecdfs_of_tics_S1 = zeros(4000, 1);
    [ecdf_tics1,ecdf_aperture1]=ecdf(reshape(squeeze(strips(1, :, :)), 1, []));
    ecdfs_of_tics_S1(aperture, 1:length(ecdf_tics1)) = ecdf_tics1;
    ecdfs_of_apertures_S1(aperture, 1:length(ecdf_aperture1)) = ecdf_aperture1;
    
    % ecdfs_of_apertures_S2 = zeros(4000, 1);
    % ecdfs_of_tics_S2 = zeros(4000, 1);
    [ecdf_tics2,ecdf_aperture2]=ecdf(reshape(squeeze(strips(2, :, end/2-20:end/2+20)), 1, []));
    ecdfs_of_tics_S2(aperture, 1:length(ecdf_tics2)) = ecdf_tics2;
    ecdfs_of_apertures_S2(aperture, 1:length(ecdf_aperture2)) = ecdf_aperture2;
end

%%% Show strips
if verbose
    index_arbitrary = 1;
    for i = 1:size(strips, 1)
        index_arbitrary = index_arbitrary + 1;
        figure(100+index_arbitrary)
        x_data = [0.0:1.0:size(strips, 3)-1]*0.6;
        y_data = [0.0:1.0:size(strips, 2)-1]*0.2;
        z_data = squeeze(strips(i, :, :));
        % subplot(1, size(strips, 1), i)
        % imagesc() %end:-1:1))) %ImgBack2OutputMatrix(img))
        % set(gca,'XColor', 'none','YColor','none')
        pcolor(x_data, y_data, z_data);
        xticks(gca, x_data(1:10:end));
        yticks(gca, y_data(1:50:end));
        xlabel(gca,'Width (mm) (B-scan)')
        ylabel(gca,'Length (mm) (number of scans)');
        
        % axis image
        daspect([1 3 1])
        shading flat;
        colorbar;
        clim(gca, gca_range);
        fontsize(15,"points")
        title(sprintf("Strip %d", i), 'FontSize',18)
        
    end
    
    clf(figure(101),"reset");
    axis image
    figure(101);
    hold on;
    iauighnsk = find(sum(ecdfs_of_apertures_S1, 2)'>0);
    % cols(iauighnsk) = ['r', 'g', 'b'];
    iauighnsk = 18
    for indeces = iauighnsk
        indy = ecdfs_of_apertures_S1(indeces, :) > 0;
        if indeces == iauighnsk(1)
            hold off;
        end
        semilogy(squeeze(ecdfs_of_apertures_S1(indeces, indy)), squeeze(ecdfs_of_tics_S1(indeces, indy)), cols(indeces), 'LineWidth',2);
        if indeces == iauighnsk(1)
            hold on;
        end
    end
    grid on;
    title('ECDF comparison (for first strip)', 'FontSize',20)
    legend('Original image', 'Strips with no blending', 'Stitched panorama', Location='northwest');
    % xlabel('Depth (mm)');
    % ylabel('ECDF (%)')
    yticks([0.00001, 0.0001, 0.001, 0.01, 0.1, 1])
    fontsize(15,'points');
    legend('Aperture 6', 'Aperture 12', 'Aperture 18');


end

%%%%%%% Algorithm Start

tic;
%%% SSIM find largest overlap
jSSIMvals = zeros(1, 50);
% Starting at 15 as feature-less space in scans can make it so a few blank columns
% look the same (e.g. two white sheets of paper) 
for t = 15:size(strips, 3)
    jSSIMvals(t) = ssim(squeeze(strips(1, :, end-t+1:end)), squeeze(strips(2, :, 1:t)));
end

[SSIM_val, j_used] = max(jSSIMvals(:), [], "all");
perceived_j = j_used;
perceived_SSIM = SSIM_val;

% img = Matrix2Img(output_matrix, -1, -1); % 5, 10
%%% Strip As-is
StripLeft =     Matrix2Img(squeeze(strips(1, :, :)) , min_depth, max_depth); %end:-1:1))
StripRight =    Matrix2Img(squeeze(strips(2, :, :)) , min_depth, max_depth); %end:-1:1))

%%% Prepare Strip
percentage = double(j_used/100);
StripLeft_padded =     imgaussfilt( prepareStrip(StripLeft, percentage, "right") , 0.1);
StripRight_padded =    imgaussfilt( prepareStrip(StripRight, percentage, "left") , 0.1);

%%% Check for all detector types
for det = dets_used %5:length(detectors)
                
    %%% Set detector
    detector = detectors(det);
    [FeatureDetector, FeatureExtractor] = getDetExtr(detector);

    %%% Prepare for matrix loop
    % TransfMatrix = defaultRIGID;
    ratio = 0.25;
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
            fontsize(15,"points")
            subplot(1, 3, 1)
            
            % imagesc(ImgBack2OutputMatrix(StripLeft_padded)) %, [L, R]);
            imagesc(StripLeft_padded) 
            axis image
            features = FeatureDetector(StripLeft_padded);
            hold on;
            plot(features.selectStrongest(30))
            hold off;
            set(gca,'XColor', 'none','YColor','none')
            colorbar;
            title('Left Strip (cropped)', 'FontSize',18);
            
            % clim(gca, gca_range);
            subplot(1, 3, 2)
            % imagesc(ImgBack2OutputMatrix(StripRight_padded))
            imagesc(StripRight_padded) 
            axis image
            features = FeatureDetector(StripRight_padded);
            hold on;
            plot(features.selectStrongest(30))
            hold off;
            set(gca,'XColor', 'none','YColor','none')
            colorbar;
            title('Right Strip (cropped)', 'FontSize',18);
            
            % clim(gca, gca_range);

            subplot(1, 3, 3)
            showMatchedFeatures(    StripLeft_padded, StripRight_padded, ...
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
            
            title(sprintf("Identified Matches = %d", length(matchedPointsRight)), 'FontSize', 18)
            
        end

        %%% Check for success
        if (sum(inlierIndex) >= 2 && ... %%% At least 2 matches found
            abs(TransfMatrix.RotationAngle - angle) < 50) %%% Angle not too large
            found_flag = 1;
            success=1;
            sgtitle('Succeded to stitch', 'FontSize', 20)
            break;
        end

        %%% Relax required match strength (loop condition) each failed loop
        ratio = ratio + 0.2;

        %%% Check for fail
        if ratio > 1
            success=0;
            fprintf("Not found for detector = %s, angle = %d, matched points = %d\n", detector, TransfMatrix.RotationAngle, length(matchedPointsLeft));
            sgtitle('Failed to stitch', 'FontSize', 20)
            
            break;
        end

                    
    end % End while
            
    %%% It either found something or nothing
    TransfMatrix.Translation = round(TransfMatrix.Translation); %e.g. [41, 0]; 
    TransfMatrix.RotationAngle = round(TransfMatrix.RotationAngle); % e.g. [0];

    %%% Use this if you want to generate a ground truth matrix and you know
    %%% the transformation already.
    % TransfMatrix.Translation = [41, 0]; %
    % TransfMatrix.RotationAngle = [0];

    %%% Generate panorama
    reconstruction_min = CreatePanorama_minOnly([defaultRIGID, TransfMatrix], Matrix2Img(squeeze(strips(1:2, :, :)) , min_depth, max_depth));

    if(verbose)
        % Show the Image Reconstruction 
        figure(4)
        % sgtitle('Image Reconstruction', 'FontSize', 20);

        imagesc(ImgBack2OutputMatrix(reconstruction_min))
        set(gca,'XColor', 'none','YColor','none')
        axis image
        colorbar;
        fontsize(15,"points")
        title('Reconstructed Panorama', 'FontSize', 18)
        % clim(gca, [5,10]);
    end
end
stop = toc;



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
    end
        
    if isequal(detector,'KAZEHoG') | isequal(detector,'SIFTHoG') | isequal(detector,'SURFHoG')
        FeatureExtractor = @(X, Y) extractHOGFeatures(X, Y, "CellSize", [4, 4], "BlockSize", [2, 2]);
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