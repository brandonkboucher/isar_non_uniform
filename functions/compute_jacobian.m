function Gamma = compute_jacobian(...
    alpha, ...
    phase, ...
    theta_m, ...
    fc, ...
    f_hat_l, ...
    x, ...
    y, ...
    u0, ...
    use_range_approx ...
    )


    % compute the range
    if use_range_approx
    
        drdx = sin(theta_m);                % [ML x 1], theta cycles
        drdy = cos(theta_m);                % [ML x 1], theta cycles

    else

        r = sqrt((x .* cos(theta_m) - y .* sin(theta_m)).^2 ...
            + (u0 + (x .* sin(theta_m) + y .* cos(theta_m))).^2);
        drdx    = (x + u0 .* sin(theta_m)) ./ r;
        drdy    = (y + u0 .* cos(theta_m)) ./ r;
        
    end

    const = Constants;
    c = const.c;
    M = size(theta_m,1);
    L = size(f_hat_l,2);

    alpha_real = real(alpha);
    alpha_imag = imag(alpha);

    % Jacobian used for NLS steps (eq 32b)
    dzdar = [cos(phase); sin(phase)];
    dzdai = [-sin(phase); cos(phase)];
    drdx = repmat(drdx, L, 1);
    drdy = repmat(drdy, L, 1);
    
    % find derivative of range wrt x,y
    kvec = repelem(4*pi*(fc + f_hat_l(:))/c, M, 1);   % [ML x 1], k_l blockwise
    
    % find the derivative of the phase wrt x,y
    d_phase_dx = -kvec .* drdx; 
    d_phase_dy = -kvec .* drdy;

    dzrdx = -d_phase_dx .* (alpha_real * sin(phase) + alpha_imag * cos(phase)); % (eq 32j)
    dzidx = d_phase_dx .* (-alpha_imag * sin(phase) + alpha_real * cos(phase)); % (eq 32k)

    dzdx = [dzrdx; dzidx]; % (eq 32i)

    dzrdy = -d_phase_dy .* (alpha_real * sin(phase) + alpha_imag * cos(phase)); % (eq 32m)
    dzidy = d_phase_dy .* (-alpha_imag * sin(phase) + alpha_real * cos(phase)); % (eq 32n)

    dzdy = [dzrdy; dzidy]; % (eq 32l)

    Gamma = [dzdar, dzdai, dzdx, dzdy]; % (eq 32b)

end