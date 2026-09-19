
% determine the derivatives of the atom:
% da/dx = -1j*\gamma*theta_m*a(x,y)
% compute_atom fills an [M x L] array and reshapes it column-major, so
% entry i of the vectorized atom is pulse m = mod(i-1,M)+1 and range
% frequency l = floor((i-1)/M)+1. gamma must therefore be held constant
% across each block of M entries (repelem), while theta cycles with
% period M (repmat). Using repmat for both silently transposes the
% range-frequency axis whenever M == L.

function [dadx, dady, d2adx2, d2ady2, d2adxdy] = compute_atom_derivatives(...
    a, ...                  % atom
    x, ...                  % (m) x position
    y, ...                  % (m) y position
    u0, ...                 % (m) target center range
    theta_m,...             % [M x 1] target yaw as a function of slow-time
    f_hat_l, ...            % [L x 1] (Hz) range-frequency
    fc, ...                 % (Hz) center frequency
    use_range_approx ...    % boolean to use r = u0 + y approximation
    )   

    % instantiate constants
    const = Constants;
    c = const.c;

    % define the dimensions of the sensing matrix
    M = size(theta_m, 1);
    L = size(f_hat_l, 2);
    
    % calculate atom derivatives for newton's method and Gauss-Newton
    gamma = -1j * 4 * pi * (fc + f_hat_l) / c;
    gamma = gamma(:);
    gamma = repelem(gamma, M, 1); % [ML x 1]

    % determine the derivative of r wrt position
    theta = repmat(theta_m, L, 1); % [ML x 1]
    if use_range_approx

        % r = u0 + x*sin(theta) + y*cos(theta) is linear in x and y, so it
        % has no curvature. the (dr/dx)^2 contribution is already carried by
        % the gamma^2 term in d2adx2 below -- putting it here as well would
        % count it twice, at the wrong power of gamma
        drdx    = sin(theta);
        drdy    = cos(theta);
        d2rdx2  = zeros(size(theta));
        d2rdy2  = zeros(size(theta));
        d2rdxdy = zeros(size(theta));

    else
        
        r = sqrt((x.*cos(theta) - y.*sin(theta)).^2 ...
            + (u0 + (x.*sin(theta) + y.*cos(theta))).^2);
        drdx    = (x + u0.*sin(theta)) ./ r;
        drdy    = (y + u0.*cos(theta)) ./ r;
        d2rdx2  = (1 ./ r) ...
            - (x + u0.*sin(theta)).^2 ./ r.^3;
        d2rdy2  = (1 ./ r) ...
            - (y + u0.*cos(theta)).^2 ./ r.^3;
        d2rdxdy = -(x + u0.*sin(theta)) ...
            .* (y + u0.*cos(theta)) ./ r.^3;

    end

    dadx =      gamma .* drdx .* a; % [ML x 1]
    dady =      gamma .* drdy .* a; % [ML x 1]
    d2adx2 =    gamma.*(gamma .* (drdx).^2 + d2rdx2) .* a; % [ML x 1]
    d2ady2 =    gamma.*(gamma .* (drdy).^2 + d2rdy2) .* a; % [ML x 1]
    d2adxdy =   gamma.*(gamma .* (drdx.*drdy) + d2rdxdy) .* a; % [ML x 1]

end

