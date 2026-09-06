function [range_resolution,cross_range_resolution] = ...
    calculate_resolution(...
        theta, ...      % [M x 1] yawing angle as function of slow-time
        f_hat_l, ...    % [L x 1] range-frequency
        fc...           % (Hz) center frequency
        )
    
    const = Constants();
    c = const.c;

    % determine the angular span which determines the cross
    % range resolution
    sin_span = max(sin(theta)) - min(sin(theta));
    cross_range_resolution = const.c ...
        / (2 * (fc + max(f_hat_l)) * sin_span);
    range_resolution = ...
        const.c / (2 * (max(f_hat_l) - min(f_hat_l))); % [m]
end

