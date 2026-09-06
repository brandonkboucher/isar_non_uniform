function A = compute_atoms_batch( ...
    x, ...                  % [K x 1] crossrange of each atom [m]
    y, ...                  % [K x 1] range of each atom [m]
    u0, ...                 % [m] distance to the center of rotation
    thetas, ...             % [M x 1] yaw angle as a function of slow time
    f_hat_l, ...            % [1 x L] range frequencies
    fc, ...                 % [Hz] center frequency
    use_range_approx, ...   % optional; false (default) = exact range
    progress_fcn ...        % optional; @(frac) called with progress in [0,1]
    )
% COMPUTE_ATOMS_BATCH  Vectorized equivalent of COMPUTE_ATOM over many atoms.
%
%   A = COMPUTE_ATOMS_BATCH(x, y, ...) returns [M*L x K], where column k is
%   identical to COMPUTE_ATOM(x(k), y(k), ...). The two share the same
%   column-major (m fastest, then l) ordering. Written as a batch so the
%   dictionary build is a handful of array operations rather than K nested
%   loops -- the difference between seconds and minutes for a dense grid.

    if nargin < 7 || isempty(use_range_approx)
        use_range_approx = false;
    end
    if nargin < 8
        progress_fcn = [];
    end

    const = Constants;
    c = const.c;

    x = x(:).';             % [1 x K]
    y = y(:).';             % [1 x K]
    thetas = thetas(:);     % [M x 1]
    f_hat_l = f_hat_l(:).'; % [1 x L]

    M = numel(thetas);
    L = numel(f_hat_l);
    K = numel(x);

    % rotate every atom position by every slow-time yaw angle -> [M x K]
    ct = cos(thetas);
    st = sin(thetas);
    xk_m = ct * x - st * y;
    yk_m = st * x + ct * y;

    % calculate the range
    if use_range_approx
        rk = u0 + yk_m;                         % Cheng et al. eq (1)
    else
        rk = sqrt((u0 + yk_m).^2 + xk_m.^2);    % exact range
    end
    tau_k = 2 * rk / c;                         % [M x K]

    % a(m,l) = exp(-1j*2*pi*(fc + f_hat_l(l))*tau_k(m)), stacked column-major
    A = zeros(M*L, K);
    for l = 1:L
        A((l-1)*M + (1:M), :) = exp(-1j * 2 * pi * (fc + f_hat_l(l)) * tau_k);

        if ~isempty(progress_fcn)
            progress_fcn(l / L);
        end
    end
end
