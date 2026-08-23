function out = isar_run_scenario(cfg, opts)
% ISAR_RUN_SCENARIO  Run one ISAR reconstruction scenario.
%
%   out = ISAR_RUN_SCENARIO(cfg, opts) synthesizes the phase history for a
%   single row of the scenario matrix, forms the image with each requested
%   algorithm, and scores the recovered scatterer positions against truth.
%   It is the body of isar_testing_v2.m lifted into a function so that both
%   the batch script and isar_gui can drive the same simulation.
%
%   cfg -- the scenario columns of documentation/scenarios_new.xlsx:
%     .NumberOfAmbiguitiesHavingScatterers  1 or 3
%     .NumberOfAmbiguitiesInImageFormer     1 or 3
%     .ScattererLocations                   'On-Grid' | 'Off-Grid'
%     .ImageFormerGridDensity               'Critical' | 'Oversampled'
%     .AngleRate                            'Constant' | 'Accelerating'
%     .Noise                                'None' or an SNR in dB (numeric
%                                           or a string like '20 dB')
%     .TargetSpacing                        'normal' | 'close'
%
%   opts -- everything that is not a scenario variable (radar parameters,
%   which algorithms to run, solver settings, seed). ISAR_DEFAULT_OPTIONS
%   returns the defaults used by isar_testing_v2.m; any field left out here
%   is filled in from there.
%
%   out contains the formed images, the estimated and true positions, the
%   per-algorithm error metrics, the plotting grids, and a sim_config struct
%   matching the one written to the results spreadsheet.

    if nargin < 2 || isempty(opts)
        opts = struct();
    end
    opts = isar_default_options(opts);

    log_fcn      = opts.log_fcn;
    progress_fcn = opts.progress_fcn;

    if ~isempty(opts.seed)
        rng(opts.seed)
    end

    %% radar parameters
    const = Constants;
    c = const.c;

    Nd  = opts.Nd;
    fc  = opts.fc;
    prf = opts.prf;
    fs  = opts.fs;

    Tp = (1/fs) * Nd;   % [s] pulse width

    % slow-time sampling and the target's yaw come from one place so the
    % GUI's trajectory preview cannot drift from the simulated motion
    traj  = isar_target_trajectory(cfg, opts);
    t_m   = traj.t_m;                    % [s] slow-time
    t_hat = (0:(1/fs):(Tp - 1/fs)).';    % [s] fast-time

    M = size(t_m,1);    % number of pulses
    L = size(t_hat,1);  % number of fast-time samples

    % define the range-frequency
    df_l = (fs/L);
    f_hat_l = (-L/2)*df_l:df_l:(L/2 - 1)*df_l;

    range_resolution = c / (2 * (max(f_hat_l) - min(f_hat_l))); % [m]

    %% target motion
    theta = traj.theta;
    w0 = traj.w0;
    w1 = traj.w1;
    w2 = traj.w2;

    % the angular span sets the cross range resolution
    sin_span = max(sin(theta)) - min(sin(theta));
    cross_range_resolution = c / (2 * (fc + max(f_hat_l)) * sin_span);

    % crossrange unambiguous extent: a scatterer beyond this extent is
    % ambiguous
    Wx = c * prf / (2 * fc * w0);

    % number of pixels within the crossrange unambiguous extent
    n_pix_per_amb = round(Wx / cross_range_resolution);

    n_amb_scat = cfg.NumberOfAmbiguitiesHavingScatterers;
    n_amb_if   = cfg.NumberOfAmbiguitiesInImageFormer;

    %% latent image grid
    switch lower(cfg.ImageFormerGridDensity)
        case 'critical'
            Kr = 1;
            Kc = 1;
        case 'oversampled'
            Kr = opts.oversampling_factor;
            Kc = opts.oversampling_factor;
        otherwise
            error('isar_run_scenario:gridDensity', ...
                'unknown ImageFormerGridDensity ''%s''', ...
                char(string(cfg.ImageFormerGridDensity)));
    end

    % redefine the grid in order to handle the number of ambiguous scatterers
    Ny = Kr * opts.N_critical;
    Nx = Kc * n_amb_if * n_pix_per_amb;

    % ensure the number of pixels is odd in both directions
    Nx = Nx + (mod(Nx,2)==0);
    Ny = Ny + (mod(Ny,2)==0);

    K = Nx * Ny; % number of vectorized samples

    range_pixel_res       = range_resolution / Kr;
    cross_range_pixel_res = cross_range_resolution / Kc;

    u0 = opts.u0; % [m] distance from the origin to the target center

    % range and crossrange grid of the latent image, centered on 0 so the
    % critical grid is a subset of the oversampled grid
    y_array = ((0:Ny-1) - (Ny-1)/2) * range_pixel_res;
    x_array = ((0:Nx-1) - (Nx-1)/2) * cross_range_pixel_res;
    [X,Y] = meshgrid(x_array,y_array);

    yk = reshape(Y, [], 1);
    xk = reshape(X, [], 1);

    %% scatterer placement
    Ks = opts.Ks; % number of point scatterers

    % assign scatterers depending on the number of ambiguities. if the number
    % of ambiguities is 1 -> [0], if 3 -> [-1,0,1], etc
    ii        = 0:(n_amb_scat-1);
    amb_index = ceil(ii/2) .* (-1).^ii;
    amb_index = sort(amb_index);

    is_off_grid = strcmpi(cfg.ScattererLocations, 'Off-grid');
    is_close    = strcmpi(cfg.TargetSpacing, 'close');

    if is_close && is_off_grid
        [target_locations, amb_of_k] = create_closely_spaced_target_scatterers( ...
            Ks, ...
            amb_index, ...
            Wx, ...
            max(y_array) - min(y_array), ...
            cross_range_pixel_res);
    else
        % round-robin so every requested ambiguity holds at least one scatterer
        amb_of_k = amb_index(mod(0:Ks-1, n_amb_scat) + 1).';

        % Columns are [crossrange, range] to match compute_atom's (x, y)
        % argument order. Crossrange is uniform within the scatterer's own
        % ambiguity; range is uniform over the grid extent (range has no
        % ambiguity structure here -- it is set by bandwidth, not the PRF).
        target_locations = [ ...
            amb_of_k * Wx + (rand(Ks,1) - 0.5) * Wx, ...
            (rand(Ks,1) - 0.5) * (y_array(end) - y_array(1)) ];
    end

    % On-Grid rows additionally snap onto the critically-sampled grid. The
    % ambiguity assignment above already happened, so on-grid scatterers are
    % distributed across ambiguities exactly like the off-grid ones.
    if ~is_off_grid
        target_locations(:,1) = round(target_locations(:,1)/cross_range_resolution) * cross_range_resolution;
        target_locations(:,2) = round(target_locations(:,2)/range_resolution)       * range_resolution;
    end

    % The latent image is supported only on the grid, which spans the first
    % n_amb_if ambiguities (same ordering as the scatterer assignment above).
    % Scatterers outside those ambiguities have no atom that represents them:
    % their true atoms correlate with the dictionary only at its sidelobe
    % floor and are not sparsely representable, so they act as unmodeled
    % interference rather than as ghosts folded into ambiguity 0.
    jj               = 0:(n_amb_if-1);
    amb_if           = ceil(jj/2) .* (-1).^jj;
    is_latent        = ismember(amb_of_k, amb_if);
    Ks_latent        = Ks;
    latent_locations = target_locations(is_latent, :);

    % determine if Doppler aliasing will occur
    theta_dot = w0 + w1*t_m + w2*t_m.^2;
    fd_max = max((2*fc/c) * abs(target_locations(:,1)) * max(abs(theta_dot)));
    doppler_aliasing = fd_max > prf/2;

    log_fcn(sprintf('W_x = %.3f m (%d critical pixels); scatterers in ambiguities %s (requested %d)', ...
        Wx, n_pix_per_amb, mat2str(unique(amb_of_k).'), n_amb_scat));
    log_fcn(sprintf('image former spans %.2f ambiguities (requested %d), covering %s', ...
        Nx*cross_range_pixel_res/Wx, n_amb_if, mat2str(sort(amb_if))));
    log_fcn(sprintf('%d of %d scatterers are representable in the latent image; %d act as interference', ...
        sum(is_latent), Ks, Ks - sum(is_latent)));
    if doppler_aliasing
        log_fcn('Doppler aliasing will occur');
    else
        log_fcn('No Doppler aliasing');
    end
    log_fcn(sprintf('grid is %d x %d (%d atoms), measurement is %d samples', ...
        Ny, Nx, K, M*L));

    %% sensing matrix and measurement
    progress_fcn(0.02, 'building sensing matrix');
    as = compute_atoms_batch(target_locations(:,1), target_locations(:,2), ...
        u0, theta, f_hat_l, fc, opts.use_range_approx);

    progress_fcn(0.05, 'building reconstruction dictionary');
    A = compute_atoms_batch(xk, yk, u0, theta, f_hat_l, fc, ...
        opts.use_range_approx, @(f) progress_fcn(0.05 + 0.35*f, 'building reconstruction dictionary'));

    log_fcn(sprintf('A has dimension %d by %d', size(A,1), size(A,2)));
    if opts.compute_rank
        progress_fcn(0.42, 'computing rank(A)');
        log_fcn(sprintf('A has a rank of %d', rank(A)));
    end

    % the measurement superposes the exact phase histories of the (off-grid)
    % scatterers; each column of `as` is one scatterer's response
    alpha_s = ones(Ks,1);   % complex scattering amplitudes (unit for now)
    y = as * alpha_s;

    % additive white complex Gaussian noise at the requested SNR
    snr_db = parse_noise(cfg.Noise);
    if ~isnan(snr_db)
        sigma2 = mean(abs(y).^2) / 10^(snr_db/10);
        y = y + sqrt(sigma2/2) * (randn(size(y)) + 1j*randn(size(y)));
        log_fcn(sprintf('added white Gaussian noise at %.1f dB SNR', snr_db));
    end

    %% image formation
    x_hat = struct();

    if opts.execute_omp
        progress_fcn(0.45, 'running OMP');
        x_hat_omp = omp_vec(y, A, Ks_latent);
        x_hat.omp.image = reshape(x_hat_omp, Ny, Nx);

        x_hat.omp.positions = extract_target_positions( ...
            x_hat.omp.image, x_array, y_array, Ks_latent, 'none');

        [x_hat.omp.error, x_hat.omp.pairs, x_hat.omp.missed, ...
            x_hat.omp.false_alarms, x_hat.omp.d] = ...
            calculate_reconstruction_error(latent_locations, x_hat.omp.positions);
    end

    if opts.execute_bp
        progress_fcn(0.70, 'running backprojection');
        x_hat_bp = A' * y;
        x_hat_bp = reshape(x_hat_bp, Ny, Nx);
        x_hat.bp.image = x_hat_bp / norm(x_hat_bp, 'fro');

        if is_off_grid
            interpolation_type = 'linear';
        else
            interpolation_type = 'none';
        end

        x_hat.bp.positions = extract_target_positions( ...
            x_hat.bp.image, x_array, y_array, Ks_latent, interpolation_type);

        [x_hat.bp.error, x_hat.bp.pairs, x_hat.bp.missed, ...
            x_hat.bp.false_alarms, x_hat.bp.d] = ...
            calculate_reconstruction_error(latent_locations, x_hat.bp.positions);
    end

    if opts.execute_nomp
        progress_fcn(0.78, 'running NOMP');
        [alpha_hat, p_hat] = nomp_vec(y, A, opts.Rs, opts.Rc, Ks_latent, ...
            xk, yk, u0, theta, f_hat_l, fc);

        x_hat.nomp.positions = p_hat.';
        x_hat.nomp.alpha = alpha_hat;

        [x_hat.nomp.error, x_hat.nomp.pairs, x_hat.nomp.missed, ...
            x_hat.nomp.false_alarms, x_hat.nomp.d] = ...
            calculate_reconstruction_error(latent_locations, x_hat.nomp.positions);
    end

    progress_fcn(0.98, 'collecting results');

    % summarize the reconstruction error for each algorithm that ran
    log_fcn(' ');
    log_fcn('  algorithm   RMS position error [m]');
    for alg = ["omp" "nomp" "bp"]
        if isfield(x_hat, alg) && isfield(x_hat.(alg), 'error')
            log_fcn(sprintf('  %-10s  %.4f', alg, x_hat.(alg).error));
        end
    end

    %% output
    out.x_hat                 = x_hat;
    out.x_array               = x_array;
    out.y_array               = y_array;
    out.u0                    = u0;
    out.theta                 = theta;
    out.t_m                   = t_m;
    out.traj                  = traj;
    out.target_locations      = target_locations;
    out.latent_locations      = latent_locations;
    out.amb_of_k              = amb_of_k;
    out.is_latent             = is_latent;
    out.cfg                   = cfg;
    out.opts                  = opts;

    out.sim_config = struct( ...
        'W_x_m',                Wx, ...
        'PixelsPerAmbiguity',   n_pix_per_amb, ...
        'Nx',                   Nx, ...
        'Ny',                   Ny, ...
        'CrossRangePixelRes_m', cross_range_pixel_res, ...
        'RangePixelRes_m',      range_pixel_res, ...
        'NumScatterers',        Ks, ...
        'NumLatentScatterers',  Ks_latent, ...
        'OversamplingFactor',   Kc, ...
        'NewtonSteps_Rs',       opts.Rs, ...
        'CyclicRefinements_Rc', opts.Rc, ...
        'DopplerAliasing',      string(doppler_aliasing));

    progress_fcn(1, 'done');
end

function snr_db = parse_noise(noise)
% PARSE_NOISE  Map the spreadsheet Noise column onto an SNR in dB.
%   Returns NaN for 'None' (or an empty/missing entry), meaning noiseless.

    snr_db = NaN;

    if isempty(noise) || (isnumeric(noise) && isnan(noise))
        return
    end

    if isnumeric(noise)
        snr_db = double(noise);
        return
    end

    s = strtrim(char(string(noise)));
    if isempty(s) || strcmpi(s, 'none')
        return
    end

    % accept '20', '20 dB', '20dB'
    v = sscanf(s, '%f');
    if isempty(v)
        error('isar_run_scenario:noise', 'unrecognized Noise value ''%s''', s);
    end
    snr_db = v(1);
end
