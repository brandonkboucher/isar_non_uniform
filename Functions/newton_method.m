function p_hat = newton_method(...
    r, x, y, u0, theta_m, f_hat_l, fc, Ns)
    
    const = Constants;
    c = const.c;

    % compute initial atom and G
    a = compute_atom(x, y, u0, theta_m, f_hat_l, fc, false);
    G_old = abs((a' * r))^2 / (a' * a);
    p_hat = [x;y];

    % extract dimensions
    L = size(f_hat_l,2);
    M = size(theta_m,1);

    for istep = 1:Ns

        % extract estimated target position
        x = p_hat(1); y = p_hat(2);

        % calculate the complex reflectivity
        alpha_k = a' * r / (a' * a);

        % determine the derivatives of the atom:
        % da/dx = -1j*\gamma*theta_m*a(x,y)
        gamma = -1j * 4 * pi * (fc + f_hat_l) / c;
        gamma = gamma(:);
        gamma = repmat(gamma, M, 1); % [ML x 1]

        % determine the derivative of r wrt position
        theta = repmat(theta_m, L, 1); % [ML x 1]
        drdx = sin(theta);
        drdy = cos(theta);
        d2rdx2 = drdx .* drdx;
        d2rdy2 = drdy .* drdy;
        d2rdxdy = sin(theta) .* cos(theta);

        dadx =      gamma       .* drdx     .* a; % [ML x 1]
        dady =      gamma       .* drdy     .* a; % [ML x 1]
        d2adx2 =    gamma.^2    .* d2rdx2   .* a; % [ML x 1]
        d2ady2 =    gamma.^2    .* d2rdy2   .* a; % [ML x 1]
        d2adxdy =   gamma.^2    .* d2rdxdy  .* a; % [ML x 1]

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

         %-------- Refinement Acceptance Condition -------------
        a = compute_atom(...
                x, y, u0, theta_m, f_hat_l, fc, false);
        G_new = abs((a' * r))^2 / (a' * a);

        if  G_new > G_old

            p_hat = [x; y];
            G_old = G_new;
            fprintf('new atom accepted.\n')
        else
            fprintf('new atom not accepted.\n')
            break
        end

    end

end

