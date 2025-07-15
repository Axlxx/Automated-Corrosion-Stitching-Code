function strips = AddMismatchToStrips(strips, mismatch_values, mismatch_intensity)
    % strips: strips to be modified
    % mismatch_values: pecentage between 0 and 1 of values that are changed
    % mismatch_intensity: how much the change occurs, in percentage 
    %   e.g. 0.25 means change of +- 25% of original value
    noStrips = numel(strips(:,1,1));
    totalValues = size(strips, 2) * size(strips, 3);
    replacedValues = int32(totalValues * mismatch_values);

    for i = 1:noStrips
        replaced_indeces = randperm(totalValues,replacedValues);
        curr_strip = squeeze(strips(i,:,:));
        for j = 1:replacedValues
            val = (1-mismatch_intensity) + (rand(1, 1)*mismatch_intensity*2);
            % before = curr_strip(replaced_indeces(j))
            % after = before*val
            curr_strip(replaced_indeces(j)) = curr_strip(replaced_indeces(j))*val;
        end

        curr_strip(curr_strip<0) = 0;
        curr_strip(curr_strip>1) = 1;

        strips(i,:,:) = curr_strip;
    end

end

