function strips = AddNoiseToStrips(strips_noisefree,noise_intensity, type)
    % Adds Noise to strips 
    strips = strips_noisefree;
    % if noise_intensity == 0
    %     return;
    % end
    noStrips = numel(strips_noisefree(:,1,1));
    for i = 1:noStrips
        if isequal(type, 'GAUSSIAN')
            strips(i,:,:) = imnoise(squeeze(strips(i,:,:)), 'gaussian', 0, noise_intensity/10000);
        elseif isequal(type, 'WHITE')
            noise = (rand(size(squeeze(strips_noisefree(i,:,:)))) - 0.5) / 100 * noise_intensity;
            strips(i,:,:) = squeeze(strips_noisefree(i,:,:)) + noise;
        else 
            printf("No noise Added")
        end
    end

    strips(strips > 1) = 1;
    strips(strips < 0) = 0;
end

