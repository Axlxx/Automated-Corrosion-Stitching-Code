% SaveStripsAsDepthMaps('C:\Users\an4223\OneDrive - Imperial College London\PhD MATLAB\Datasets\Realistic depths\22July_2.mat', strips);
function depth_strips = SaveStripsAsDepthMaps(folder_file, strips)
%SAVESTRIPSASDEPTHMAPS Summary of this function goes here
%   Detailed explanation goes here
    depth_strips = double(size(strips));
    for i = 1:size(strips,1)
        for j = 1:size(strips,2)
            for k = 1:size(strips,3)
                depth_strips(i, j, k) = strips(i, j, k).depth;
            end
        end
    end

    if ~isempty(folder_file)
        save(folder_file, 'depth_strips');
    end
end

