function pos = extract_target_positions(...
    image, ...
    x_array, ...
    y_array, ...
    K, ... % number of scatterers
    interpolation_type ... % 'none' or 'linear'
    )

    if ~strcmp(interpolation_type, 'none') ...
            && ~strcmp(interpolation_type, 'linear')
        error('interpolation type not found')
    end

    if strcmp(interpolation_type, 'none') % assumes omp

        % If there isn't interpolation simply find the
        % maximum on-grid scatterer values
        [sorted_vals, scatterer_idx] = sort(abs(image(:)),"descend");
        
        % extract only the non-zero values
        scatterer_idx = scatterer_idx(sorted_vals > 0);

        % if there are more non-zero values than scatterers,
        % extract the top K scatterers
        if size(scatterer_idx,1) > K
            scatterer_idx = scatterer_idx(1:K);
        end

        % convert to row/col to find grid point
        [y_idx, x_idx] = ind2sub(size(image), scatterer_idx);

        % form the position (x,y) pairs
        pos = [x_array(x_idx).', y_array(y_idx).'];


    elseif strcmp(interpolation_type, 'linear') % assumes bp and/or dense

        % find the maximum on-grid scatterer values
        [sorted_vals, scatterer_idx] = sort(abs(image(:)),"descend");
        
        % if there are more non-zero values than scatterers,
        % extract the top K scatterers
        if size(scatterer_idx,1) > K
            scatterer_idx = scatterer_idx(1:K);
        end

        % convert to row/col to find grid point
        [rows, cols] = ind2sub(size(image), scatterer_idx);

        pos = zeros(K,2);
        for iscatterer = 1:K

            % extract the candidate scatterer on-grid index
            row = rows(iscatterer);
            col = cols(iscatterer);

            % extract the corresponding position
            xpos = x_array(col);
            ypos = y_array(row);

            % ensure the local patch is entirely within the
            % grid
            row_min = max(1,row-1);
            row_max = min(size(image,1),row+1);
            col_min = max(1,col-1);
            col_max = min(size(image,2),col+1);

            % extract the local neighborhood of the
            % scatterer
            local_patch = abs(image(row_min:row_max, col_min:col_max));

            % find the weighted average position of the
            % target scatterer along each dimension
            local_patch_col = mean(local_patch,1) ...
                ./ sum(mean(local_patch,1)); 
            local_patch_row = mean(local_patch,2) ...
                ./ sum(mean(local_patch,2));

            % weighted average
            xpos = sum(x_array(col_min:col_max) .* local_patch_col);
            ypos = sum(y_array(row_min:row_max) .* local_patch_row.');
            pos(iscatterer, :) = [xpos, ypos];

        end
        

    end
end

