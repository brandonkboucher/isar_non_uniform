
% From the initial atom selected in the image former, check
% whether any of the atom's ambiguous doppelgangers have a 
% higher correlation

% NOTE: We should probably add cyclical refinement in future

function [selected_a, selected_ambiguity, selected_x] = check_ambiguity_of_an_atom(...
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

    ii        = 0:(num_of_amb-1);
    amb_index = ceil(ii/2) .* (-1).^ii;
    amb_index = sort(amb_index);

    % Wx is not an integer number of pixels, and a scatterer's ghost does not
    % sit exactly at x + Wx either, so the atom the image former selected is
    % already offset from the peak it was fitting. Shifting that atom by Wx
    % carries the offset into every other ambiguity, where it is pure coherence
    % loss -- which biases the test toward the ambiguity the grid search
    % already optimised. Measured on a two-scatterer scenario, the true
    % ambiguity lost a single-point test by 0.04% and won a locally refined
    % one by 0.68%. So give every ambiguity the same freedom: maximise over a
    % scan of one pixel either side rather than testing one point.
    % The scan half-width must exceed the ghost's displacement from x + Wx.
    % That displacement grows with yaw acceleration, because acceleration
    % changes the effective fold distance rather than destroying the alias:
    % measured 0.2 pixels at w1 = 0 and 4.25 pixels at w1 = 170 rad/s/s.
    if isfield(options, 'amb_refine_pixels')
        n_px = options.amb_refine_pixels;
    else
        n_px = 5;
    end
    if nargin < 5 || isempty(px) || px <= 0 || n_px <= 0
        x_offsets = 0;              % single-point test (original behaviour)
    else
        x_offsets = (-n_px*px):(px/10):(n_px*px);
    end

    best_c             = -inf;
    selected_a         = [];
    selected_ambiguity = 0;
    selected_x         = x;

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

            if ck > best_c
                best_c             = ck;
                selected_a         = ak;
                selected_ambiguity = amb_index(iamb);
                selected_x         = xk;
            end
        end
    end

end

