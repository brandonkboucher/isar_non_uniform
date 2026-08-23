function [err,pairs,missed,false_alarms,d] = calculate_reconstruction_error(...
    true_locations, ...
    estimated_locations ...
    )

    % number of scatterers
    K = size(true_locations, 1);

    % make sure the array is in the expected dimensions
    if size(estimated_locations, 1) == 2 ...
            && size(estimated_locations, 2) ~= 2

        estimated_locations = estimated_locations.';

    elseif size(estimated_locations, 1) ~= 2 ...
            && size(estimated_locations, 2) ~= 2

        error('estimate location dimension is not (x,y)')

    end

    % euclidean distance matrix, cost matrix. computed directly rather than
    % with pdist2 so that scoring does not depend on the Statistics toolbox
    % being licensed.
    dx = true_locations(:,1) - estimated_locations(:,1).';
    dy = true_locations(:,2) - estimated_locations(:,2).';
    C = sqrt(dx.^2 + dy.^2);

    % find the matches using MATLAB's matchpairs function
    [pairs, missed, false_alarms] = matchpairs(C, 1e10);

    % calculate the total error 
    d = C(sub2ind(size(C), pairs(:,1), pairs(:,2)));
    err = sqrt(mean(d.^2));
end

