function [TransfMatrix, matchedPointsLeft, matchedPointsRight, inlierIndex] = TransfMatrixOf2Strips(StripLeft, StripRight, FeatureDetector, FeatureExtractor, FeatureMatcher, MatrixType)
    % Detect Features Left Strip
    pointsLeft = FeatureDetector(StripLeft);
    [featuresLeft, pointsLeft] = FeatureExtractor(StripLeft,pointsLeft);

    % Detect Features Right Strip 
    pointsRight = FeatureDetector(StripRight);
    [featuresRight, pointsRight] = FeatureExtractor(StripRight,pointsRight);

    % Detect Matching features
    indexPairs = FeatureMatcher(featuresLeft, featuresRight); % matchFeatures(featuresLeft, featuresRight, 'Unique', true);   
    matchedPointsLeft   = pointsLeft(indexPairs(:,1), :);
    matchedPointsRight  = pointsRight(indexPairs(:,2), :);  

    % Get Transformation Matrix

    if(isequal(MatrixType, affinetform2d))
        type = 'affine';
    elseif (isequal(MatrixType, rigidtform2d))
        type = 'rigid';
    elseif (isequal(MatrixType, simtform2d))
        type = 'similarity';
    elseif (isequal(MatrixType, transltform2d))
        type = 'translation';
    else 
        fprintf("Weird Transformation type provided\n");
        return;
    end

    if length(indexPairs(:,1)) >= 2
        try
            [TransfMatrix, inlierIndex] = estgeotform2d(matchedPointsRight, matchedPointsLeft,...
                    type , 'Confidence', 99.9, 'MaxNumTrials', 10000, 'MaxDistance',2);
        catch
            % fprintf('Failed to produce transf matrix (not enough inliers)\n')
            TransfMatrix = MatrixType;
            inlierIndex = false(length(matchedPointsLeft), 1);
        end
        % fprintf(['Tform matrix returned error:', status])
    else
        TransfMatrix = MatrixType;
        inlierIndex = false(length(matchedPointsLeft), 1);
        % fprintf("Not enough matched points found\n");
    end
end

