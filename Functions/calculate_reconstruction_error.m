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

    % euclidean distance matrix, cost matrix
    C = pdist2(true_locations, estimated_locations, "euclidean");

    % find the matches using MATLAB's matchpairs function
    [pairs, missed, false_alarms] = matchpairs(C, 1e10);

    % calculate the total error 
    d = C(sub2ind(size(C), pairs(:,1), pairs(:,2)));
    err = sqrt(mean(d.^2));
end

