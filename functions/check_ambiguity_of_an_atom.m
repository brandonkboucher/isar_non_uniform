
% From the initial atom selected in the image former, check
% whether any of the atom's ambiguous doppelgangers have a 
% higher correlation

% NOTE: We should probably add cyclical refinement in future

function [selected_a, selected_ambiguity, selected_x, selected_y] = check_ambiguity_of_an_atom(...
    r, ...          % [ML x 1] residual
    x, ...          % (m) crossrange position of selected atom
    y, ...          % (m) range position of selected atom
    Wx, ...         % (m) unambiguous crossrange extent
    px, ...         % (m) crossrange pixel size; 0 disables local refinement
    num_of_amb, ... % the number of ambiguities to check
    u0, ...         % (m) distance from radar to target rotation axis
    theta_m, ...    % [M x 1] yawing angle as function of slow-time
    f_hat_l, ...    % [L x 1]
    fc, ...         % (Hz) center frequency
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
    if isfield(options, 'amb_refine_pixels')
        n_px = options.amb_refine_pixels;
    else
        n_px = 5;
    end

    % keep in mind, as a consequence we are moving from
    % on-grid to off-grid
    if nargin < 5 || isempty(px) || px <= 0 || n_px <= 0
        x_offsets = 0;              % single-point test (original behaviour)
    else
        x_offsets = (-n_px*px):(px/10):(n_px*px);
    end

    best_c             = -inf;
    selected_a         = [];
    selected_ambiguity = 0;
    selected_x         = x;

    % iterate over each ambiguity and each neighborhood
    for iamb = 1:num_of_amb
        for io = 1:numel(x_offsets)

            xk = x + amb_index(iamb)*Wx + x_offsets(io);

            ak = compute_atom(...
                xk, ...
                y, ...
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
                selected_x         = xk;
            end
        end
    end
    selected_y = y;
end

