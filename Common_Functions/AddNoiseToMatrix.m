function strips = AddNoiseToMatrix(strips_noisefree,noise_intensity, type)
    % Adds Noise to strips 
    strips = strips_noisefree;

    noStrips = size(strips_noisefree, 1);
    for i = 1:noStrips
        % maxim = max(strips_noisefree(i,:,:), [], 'all');
        % minim = min(strips_noisefree(i,:,:), [], 'all');

        if isequal(type, 'GAUSSIAN')
            strips(i,:,:) = imnoise(squeeze(strips(i,:,:)), 'gaussian', 0, noise_intensity/100);
        elseif isequal(type, 'WHITE')
            noise = (rand(size(squeeze(strips_noisefree(i,:,:)))) - 0.5) * noise_intensity;
            strips(i,:,:) = squeeze(strips_noisefree(i,:,:)) + noise;
        else 
            printf("No noise Added")
        end

        % strips(strips > maxim) = maxim;
        % strips(strips < minim) = minim;
    end
end

