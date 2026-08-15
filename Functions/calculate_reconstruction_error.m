function err = calculate_reconstruction_error(...
    true_locations, ...
    estimated_locations ...
    )

    fprintf('Calculating error.\n')

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

    % squared l2-norm
    C = pdist2(true_locations, estimated_locations, "squaredeuclidean");

    % find the matches using MATLAB's matchpairs function
    M = matchpairs(C, 1e10);

    % calculate the total error 
    err = sum(C(sub2ind(size(C), M(:,1), M(:,2))));

end

