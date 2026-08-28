function [a, phase, dadx, dady, d2adx2, d2ady2, d2adxdy] = compute_atom(...
    x, ...
    y, ...
    u0, ...
    theta_m,...
    f_hat_l, ...
    fc, ...
    use_range_approx)   % optional; false (default) = exact range

    % USE_RANGE_APPROX selects the range model:
    %   false (default) : exact  rk = sqrt((u0 + yk_m)^2 + xk_m^2)
    %   true            : the linearized range of Cheng et al. (2019),
    %                     IEEE Access 7:157019, eq (1):
    %                     R(t) ~ R0 + y_k*cos(dtheta) + x_k*sin(dtheta),
    %                     i.e. rk = u0 + yk_m (drop the sqrt / the
    %                     r^2/(2*R0) cross-range term). This is the
    %                     approximation under which the (f,beta) chirp
    %                     model and the off-grid OMP of that paper hold.
    if nargin < 7 || isempty(use_range_approx)
        use_range_approx = false;
    end

    % instantiate constants
    const = Constants;
    c = const.c;

    % define the dimensions of the sensing matrix
    M = size(theta_m, 1);
    L = size(f_hat_l, 2);
    a = zeros(M,L);
    phase = zeros(M,L);

    for m = 1:M

        % define the rotation matrix
        theta = theta_m(m);
        R = [cos(theta), -sin(theta); ...
                 sin(theta), cos(theta)];
    
        % apply the rotation matrix and extract the
        % positional components
        uk_m = R * [x; y];
        xk_m = uk_m(1);
        yk_m = uk_m(2);
    
        % calculate the range
        if use_range_approx
            rk = u0 + yk_m;                       % Cheng et al. eq (1)
        else
            rk = sqrt((u0 + yk_m)^2 + xk_m^2);    % exact range
        end
        tau_k = 2 * rk / c;
    
        for l = 1:L
            
            % extract the range frequency
            f_hat = f_hat_l(l);
            
            % form the sensing matrix
            phase(m,l) = -2 * pi  * (fc + f_hat) * tau_k;
            a(m,l) = exp(1j * phase(m,l));
        end
    end

    % reshape the sensing matrix to a vector
    a = reshape(a, [], 1); % M*L x 1
    phase = reshape(phase, [], 1); % M*L x 1

    % calculate atom derivatives for newton's method and Gauss-Newton
    gamma = -1j * 4 * pi * (fc + f_hat_l) / c;
    gamma = gamma(:);
    gamma = repelem(gamma, M, 1); % [ML x 1]

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

end

