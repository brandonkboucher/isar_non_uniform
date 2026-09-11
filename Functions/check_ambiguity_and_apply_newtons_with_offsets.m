
% From the initial atom selected in the image former, check
% whether any of the atom's ambiguous doppelgangers have a 
% higher correlation

% NOTE: We should probably add cyclical refinement in future

function [selected_a, selected_ambiguity,...
    selected_x, selected_y] = check_ambiguity_and_apply_newtons_with_offsets(...
    r, ...          % [ML x 1] residual
    x0, ...         % (m) crossrange position of selected atom
    y0, ...         % (m) range position of selected atom
    Wx, ...         % (m) unambiguous crossrange extent
    num_of_amb, ... % the number of ambiguities to check
    u0, ...         % (m) distance from radar to target rotation axis
    theta_m, ...    % [M x 1] yawing angle as function of slow-time
    f_hat_l, ...    % [L x 1]
    fc, ...         % (Hz) center frequency
    Rs, ...         % number of refinement steps
    cross_range_resolution, ... % (m) crossrange pixel size
    num_offsets_pixels, ... % number of pixels to traverse using Newton's method
    options ...     % additional scenario options
    )

    % define the ambiguities
    ii        = 0:(num_of_amb-1);
    amb_index = ceil(ii/2) .* (-1).^ii;
    amb_index = sort(amb_index);

    % as a result of the acceleration of the target, Wx
    % isn't constant, meaning the aliased scatterer is
    % precisely located at x +/- ambiguity * Wx. To find,
    % the ambiguous target scatterer we need to search in
    % its local neighborhood

    x_offsets = ...
        (-num_offsets_pixels*cross_range_resolution)...
        :(cross_range_resolution)...
        :(num_offsets_pixels*cross_range_resolution);

    
    best_c             = -inf;
    selected_a         = [];
    selected_ambiguity = 0;
    selected_x         = x0;
    selected_y         = y0;

    % iterate over each ambiguity and each neighborhood
    for iamb = 1:num_of_amb
        for ioffset = 1:numel(x_offsets)

            % define the cross-range using the ambiguity
            xk = x0 + amb_index(iamb) * Wx + x_offsets(ioffset);
    
            % for each ambiguity and each offset perform
            % Newton's method as an additional refinement
            p_hat = newton_method(...
                r, ... % [ML x  1] residual
                xk, ... % (m) cross-range position
                y0, ... % (m) range position
                u0, ... % (m) radar-to-target center distance
                theta_m, ... % [M x 1] (rad) target yaw
                f_hat_l, ... % [L x 1] (Hz) range-frequency
                fc, ... % (Hz) center frequency
                Rs, ... % number of refinement steps
                options.use_range_approx... % range approx boolean
                );
    
            % compute optimized atom
            ak = compute_atom(...
                p_hat(1), ...
                p_hat(2), ...
                u0, ...
                theta_m, ...
                f_hat_l, ...
                fc, ...
                options.use_range_approx);
    
            % inner product of the residual and this candidate
            ck = abs(ak' * r);
    
            % save the atom that has the highest correlation
            if ck > best_c
                best_c             = ck;
                selected_a         = ak;
                selected_ambiguity = amb_index(iamb);

                % both coordinates come from the winning refinement. Assigning
                % the range outside the loops would report whichever seed
                % happened to be evaluated last, not the one selected.
                selected_x         = p_hat(1);
                selected_y         = p_hat(2);
            end
        end
    end
end

