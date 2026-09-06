function [p_hat, p_hat_hist] = newton_method(...
    r, ...
    x, ...
    y, ...
    u0, ...
    theta_m, ...
    f_hat_l, ...
    fc, ...
    Ns, ...
    use_range_approx)
    
    const = Constants;
    c = const.c;

    % compute initial atom and G
    a = compute_atom(x, y, u0, theta_m, f_hat_l, fc, use_range_approx);
    G_old = abs((a' * r))^2 / (a' * a);
    p_hat = [x;y];
    
    % save estimated position. Ns is the upper bound; the loop below can stop
    % early, so nfilled tracks how many columns are real and the array is
    % trimmed to that before returning -- an unfilled column would otherwise
    % read as an estimate at the origin
    p_hat_hist = zeros(2, Ns);
    nfilled = 0;

    % extract dimensions
    L = size(f_hat_l,2);
    M = size(theta_m,1);

    for istep = 1:Ns

        % extract estimated target position
        x = p_hat(1); y = p_hat(2);

        % calculate the complex reflectivity
        alpha_k = a' * r / (a' * a);

        % determine the derivative of r wrt position
        [dadx, dady, d2adx2, d2ady2, d2adxdy] = compute_atom_derivatives(...
            a, ...                  % atom
            x, ...                  % (m) x position
            y, ...                  % (m) y position
            u0, ...                 % (m) target center range
            theta_m,...             % [M x 1] target yaw as a function of slow-time
            f_hat_l, ...            % [L x 1] (Hz) range-frequency
            fc, ...                 % (Hz) center frequency
            use_range_approx ...    % boolean to use r = u0 + y approximation
            );

        % newton's method for optimization:
        % p = p - H^{-1}F, F (eq 7) and H (eq 8), p =
        % [x, y]^{T}
        F(1,1) = real((r - a * alpha_k)' * alpha_k * dadx);
        F(2,1) = real((r - a * alpha_k)' * alpha_k * dady);

        % define the Hessian
        H(1,1) = real((r - a * alpha_k)'* alpha_k * d2adx2) ...
            - abs(alpha_k)^2 * (dadx' * dadx);
        H(2,2) = real((r - a * alpha_k)'* alpha_k * d2ady2) ...
            - abs(alpha_k)^2 * (dady' * dady);
        H(1,2) = real((r - a * alpha_k)'* alpha_k * d2adxdy) ...
            - abs(alpha_k)^2 * real(dadx' * dady);
        H(2,1) = H(1,2);

        % update positional approximation
        p = [x; y];
        if ~(trace(H) < 0 && det(H) > 0), break; end
        p = p - inv(H)*F;
        x = p(1,1); y = p(2,1);
        p_hat_hist(:,istep) = [x;y];
        nfilled = istep;
        
         %-------- Refinement Acceptance Condition -------------
        a = compute_atom(...
                x, y, u0, theta_m, f_hat_l, fc, use_range_approx);
        G_new = abs((a' * r))^2 / (a' * a);

        if  G_new > G_old

            p_hat = [x; y];
            G_old = G_new;
            
        else
            break
        end

    end

    % hand back only the steps that were actually taken
    p_hat_hist = p_hat_hist(:, 1:nfilled);

end
