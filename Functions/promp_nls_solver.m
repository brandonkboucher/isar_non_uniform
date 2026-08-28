
function [p_hat, alpha_hat] = promp_nls_solver(...
    measurement, ...  % measurement
    p, ...          % previous estimates (i-1-th) [2 x i-1]
    x0, ...         % new coarse estimate in x-direction
    y0, ...         % new coarse estimate in y-direction
    alpha, ...      % previous reflectivity estimate
    i, ...          % estimate index 
    h_max, ...      % max number of Gauss-Newton iterations
    u0, ...         % [m] distance from radar to scene center
    theta_m, ...    % [rad] angle as a function of slow-time
    f_hat_l, ...    % [Hz] range-frequency
    fc, ...         % [Hz] radar center frequency
    zeta, ...       % [2 x 1] constraint bound, [crossrange; range]
    p_anchor ...    % [2 x i] grid nodes the constraint is measured from
    )

    const = Constants;
    c = const.c;

    M = size(theta_m,1);
    L = size(f_hat_l,2);

    % the goal of this script is to perform a non-linear least squares
    % estimate of the parameters, complex reflectivity and position. i 
    % am going to try to use the notation of the promp paper where 

    % xi^{i} = [xi_{1}^{i}, ..., xi_{k}^{i} ..., xi_{i}^{i}] (eq 17a)
    % xi_k^{i} = [\alpha^{R}_{k}, \alpha^{I}_{k}, x_{k}, y_{k}] (eq 17b)

    % where \alpha^{R} is the real-component of the reflectivity and 
    % \alpha^{I} is the imaginary-component

    % --------------------------------------
    % build z_tilde

    z_tilde = [real(measurement); imag(measurement)]; % (eq 15)

    % --------------------------------------
    % seed xi^{i} once, from the converged solution of iteration i-1 stacked
    % with the newly selected atom (eq 23, eq 24). every Gauss-Newton step
    % below then reads its parameters back out of xi.
    xi = zeros(4*i,1);

    for k = 1:i

        if k ~= i
            x = p(1,k);
            y = p(2,k);
        else
            x = x0;
            y = y0;
        end

        xi(4 * (k-1) + 1) = real(alpha(k));
        xi(4 * (k-1) + 2) = imag(alpha(k));
        xi(4 * (k-1) + 3) = x;
        xi(4 * (k-1) + 4) = y;

    end

    % indices of the position entries of xi, used for the stopping rule
    pos_idx = reshape([3:4:4*i; 4:4:4*i], [], 1);

    Gamma = zeros(2*M*L, 4*i);

    % --------------------------------------
    % perform Gauss-Newton steps
    for h = 1:h_max

        % z and Gamma are both evaluated at the current iterate, so they are
        % rebuilt over every atom on every step
        z = zeros(2*M*L, 1);
        A = zeros(M*L, i);

        for k = 1:i

            alpha_real = xi(4 * (k-1) + 1);
            alpha_imag = xi(4 * (k-1) + 2);

            % determine the phase from the estimated position
            [a, phase] = ...
                compute_atom(...
                    xi(4 * (k-1) + 3), ... % x
                    xi(4 * (k-1) + 4), ... % y
                    u0, theta_m, f_hat_l, fc, false);
            
            % divide into real and imaginary components
            zi = [alpha_real * cos(phase) - alpha_imag * sin(phase); ... % (eq30a)
                alpha_real * sin(phase) + alpha_imag * cos(phase)]; % (eq 30b)

            % sum across all atoms
            z = z + zi; % (eq 16)

            Gamma(:,4 * (k-1) + 1:4 * (k-1) + 4 ) = compute_jacobian(...
                alpha_real + 1j * alpha_imag, ...
                phase, ...
                theta_m, ...
                fc, ...
                f_hat_l);
            
            % add to sensing matrix formulation
            A(:,k) = a;

        end

        % Gauss-Newton update, solved by QR rather than by forming the
        % normal equations explicitly
        delta = Gamma \ (z_tilde - z); % (eq 19)
        xi = xi + delta;

        % --------------------------------------
        % check the new estimates against the constraint, measured from the
        % grid node that was selected for each atom (eq 20a, eq 20b)
        clamped = false;

        for k = 1:i

            xG = p_anchor(1,k);
            yG = p_anchor(2,k);

            x_hat = xi(4 * (k-1) + 3);
            y_hat = xi(4 * (k-1) + 4);

            if abs(x_hat - xG) > zeta(1)
                xi(4 * (k-1) + 3) = xG + sign(x_hat - xG) * zeta(1);
                clamped = true;
            end
            if abs(y_hat - yG) > zeta(2)
                xi(4 * (k-1) + 4) = yG + sign(y_hat - yG) * zeta(2);
                clamped = true;
            end

        end

        % re-estimate the reflectivities by least squares, but only when the
        % constraint is actually in effect (eq 21). rebuilding A first keeps
        % the reflectivities consistent with the clamped positions
        if clamped

            for k = 1:i
                A(:,k) = compute_atom(...
                    xi(4 * (k-1) + 3), ...
                    xi(4 * (k-1) + 4), ...
                    u0, theta_m, f_hat_l, fc, false);
            end

            alpha = A \ measurement;

            for k = 1:i
                xi(4 * (k-1) + 1) = real(alpha(k));
                xi(4 * (k-1) + 2) = imag(alpha(k));
            end

        end

        % terminate once the position update is small relative to the
        % positions themselves
        if norm(delta(pos_idx)) / norm(xi(pos_idx)) < 1e-3
            break
        end

    end
    
    % extract the final estimates (eq 22)
    p_hat = zeros(2, i);
    alpha_hat = zeros(i, 1);
    
    for k = 1:i
        alpha_hat(k) = xi(4 * (k-1) + 1) + 1j * xi(4 * (k-1) + 2);
        p_hat(:,k) = [xi(4 * (k-1) + 3); xi(4 * (k-1) + 4)];
    end

    
end
