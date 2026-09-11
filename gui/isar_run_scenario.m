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
%   per-algorithm error metrics, the plotting grids, the measurement vector y,
%   and a sim_config struct matching the one written to the results
%   spreadsheet. It also carries the reconstruction dictionary A and the
%   measurement sensing matrix as, so a caller can save or reuse them without
%   having to have asked for them before the run.

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
    M     = size(t_m,1);                 % number of pulses

    range_cell_mode = opts.range_cell_mode;

    if range_cell_mode
        % One range cell: a single fast-time sample and a single range
        % frequency. There is no bandwidth left to resolve range with, so the
        % latent image collapses from a picture to a line in crossrange, and
        % every scatterer sits in the one cell being imaged. Doppler is the
        % only thing separating them, which is the point -- it isolates the
        % ambiguity question from the range dimension entirely.
        t_hat   = 0;
        L       = 1;
        f_hat_l = 0;

        % no range resolution exists; the single cell is the whole extent
        range_resolution = Inf;
    else
        t_hat = (0:(1/fs):(Tp - 1/fs)).';    % [s] fast-time
        L     = size(t_hat,1);               % number of fast-time samples

        df_l    = (fs/L);
        f_hat_l = (-L/2)*df_l:df_l:(L/2 - 1)*df_l;

        range_resolution = c / (2 * (max(f_hat_l) - min(f_hat_l))); % [m]
    end

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

    if range_cell_mode
        Ny = 1;                     % one cell, so one row
        Kr = 1;
        range_pixel_res = 1;        % nothing to scale; kept finite for the
                                    % reports and the plotting limits
    else
        Ny = Kr * opts.N_critical;
        Ny = Ny + (mod(Ny,2)==0);
        range_pixel_res = range_resolution / Kr;
    end

    cross_range_pixel_res = cross_range_resolution / Kc;

    u0 = opts.u0; % [m] distance from the origin to the target center

    % range grid, centered on 0 so the critical grid is a subset of the
    % oversampled grid. Only crossrange depends on the ambiguity span.
    if range_cell_mode
        y_array = opts.range_cell;      % the one cell being imaged
    else
        y_array = ((0:Ny-1) - (Ny-1)/2) * range_pixel_res;
    end

    %% ambiguity span per algorithm
    % Each algorithm can be given its own number of ambiguities in the image
    % former, so that mod-OMP working from a one-block dictionary can be
    % compared against an OMP handed the full three-block dictionary -- the
    % comparison the computational argument rests on. Algorithms that share a
    % span share one dictionary; only a distinct span costs a second build.
    alg_names  = {'omp', 'mod_omp', 'nomp', 'promp', 'bp'};
    alg_enabled = [opts.execute_omp, opts.execute_mod_omp, opts.execute_nomp, ...
        opts.execute_promp, opts.execute_bp];

    amb_if_alg = struct();
    for ia = 1:numel(alg_names)
        amb_if_alg.(alg_names{ia}) = n_amb_if;
    end
    if isfield(opts, 'amb_in_image_former') && isstruct(opts.amb_in_image_former)
        for fn = fieldnames(opts.amb_in_image_former).'
            a = fn{1};
            if ~ismember(a, alg_names)
                error('isar_run_scenario:ambOverride', ...
                    'unknown algorithm ''%s'' in amb_in_image_former', a);
            end
            v = opts.amb_in_image_former.(a);
            if ~isempty(v) && isfinite(v) && v >= 1
                amb_if_alg.(a) = v;
            end
        end
    end

    % the spans actually needed: the scenario's own (it is what the reports
    % and out.x_array describe) plus every enabled algorithm's
    spans = n_amb_if;
    for ia = 1:numel(alg_names)
        if alg_enabled(ia)
            spans(end+1) = amb_if_alg.(alg_names{ia}); %#ok<AGROW>
        end
    end
    spans = unique(spans);

    grids = struct();
    for v = spans(:).'
        Nxv = Kc * v * n_pix_per_amb;
        Nxv = Nxv + (mod(Nxv,2)==0);        % keep the number of pixels odd

        x_arr    = ((0:Nxv-1) - (Nxv-1)/2) * cross_range_pixel_res;
        [Xv, Yv] = meshgrid(x_arr, y_array);

        grids.(span_key(v)) = struct( ...
            'n_amb_if', v, ...
            'Nx',       Nxv, ...
            'K',        Nxv * Ny, ...
            'x_array',  x_arr, ...
            'xk',       reshape(Xv, [], 1), ...
            'yk',       reshape(Yv, [], 1), ...
            'A',        []);
    end

    % the scenario's own span is what the reports, out.x_array and out.A
    % describe, so keep those names bound to it
    g_default = grids.(span_key(n_amb_if));
    Nx = g_default.Nx;  K = g_default.K;
    x_array = g_default.x_array;
    xk = g_default.xk;  yk = g_default.yk;

    % A missed scatterer is charged the crossrange width of the imaged scene.
    % It has to be the same charge for every algorithm or the comparison tilts
    % toward whoever was given the narrower grid, so it comes from the widest.
    g_widest = grids.(span_key(max(spans)));
    miss_penalty = max(g_widest.x_array) - min(g_widest.x_array);

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

    % with a single range cell every scatterer is in it by construction, so
    % the range coordinate is the cell itself rather than a random draw
    if range_cell_mode
        target_locations(:,2) = opts.range_cell;
    end

    % On-Grid rows additionally snap onto the critically-sampled grid. The
    % ambiguity assignment above already happened, so on-grid scatterers are
    % distributed across ambiguities exactly like the off-grid ones.
    if ~is_off_grid
        target_locations(:,1) = round(target_locations(:,1)/cross_range_resolution) * cross_range_resolution;
        if ~range_cell_mode
            target_locations(:,2) = round(target_locations(:,2)/range_resolution) * range_resolution;
        end
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

    % Algorithms can be given different ambiguity spans, so "which scatterers
    % are representable" is a per-span question, not a per-scenario one. The
    % plots need it per algorithm: a panel spanning three ambiguities has to
    % draw the truth in all three, not only the band the scenario named.
    latent_by_span = struct();
    for v = spans(:).'
        jv = 0:(v-1);
        latent_by_span.(span_key(v)) = ...
            ismember(amb_of_k, ceil(jv/2) .* (-1).^jv);
    end

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
    for ia = 1:numel(alg_names)
        if alg_enabled(ia) && amb_if_alg.(alg_names{ia}) ~= n_amb_if
            log_fcn(sprintf('  %s uses %d ambiguities in the image former', ...
                alg_names{ia}, amb_if_alg.(alg_names{ia})));
        end
    end

    %% sensing matrix and measurement
    progress_fcn(0.02, 'building sensing matrix');
    as = compute_atoms_batch(target_locations(:,1), target_locations(:,2), ...
        u0, theta, f_hat_l, fc, opts.use_range_approx);

    n_spans = numel(spans);
    for isp = 1:n_spans
        v = spans(isp);
        key = span_key(v);
        g   = grids.(key);

        if v == 1
            msg = 'building dictionary (1 ambiguity)';
        else
            msg = sprintf('building dictionary (%d ambiguities)', v);
        end
        base = 0.05 + 0.35*(isp-1)/n_spans;
        progress_fcn(base, msg);

        g.A = compute_atoms_batch(g.xk, g.yk, u0, theta, f_hat_l, fc, ...
            opts.use_range_approx, ...
            @(f) progress_fcn(base + 0.35*f/n_spans, msg));
        grids.(key) = g;

        log_fcn(sprintf('dictionary for %d-ambiguity grid: %d by %d', ...
            v, size(g.A,1), size(g.A,2)));
        if opts.compute_rank
            progress_fcn(0.42, 'computing rank(A)');
            log_fcn(sprintf('  rank %d', rank(g.A)));
        end
    end

    A = grids.(span_key(n_amb_if)).A;   % the scenario's own dictionary

    % the measurement superposes the exact phase histories of the (off-grid)
    % scatterers; each column of `as` is one scatterer's response
    % complex scattering amplitudes, identical across scatterers. The noise
    % power below is set from mean(|y|^2), so the requested SNR is preserved
    % whatever magnitude is chosen -- this scales the signal, not the contrast
    % against the noise.
    alpha_s = opts.target_magnitude * ones(Ks,1);
    y = as * alpha_s;

    % additive white complex Gaussian noise at the requested SNR
    snr_db = parse_noise(cfg.Noise);
    opts.residual_threshold = [];
    if ~isnan(snr_db)
        sigma2 = mean(abs(y).^2) / 10^(snr_db/10);
        y = y + sqrt(sigma2/2) * (randn(size(y)) + 1j*randn(size(y)));
        log_fcn(sprintf('added white Gaussian noise at %.1f dB SNR', snr_db));

        % Nguyen et al. stop every algorithm when the signal residual reaches
        % the noise level (Sec IV-A) rather than at a fixed atom count. For
        % per-sample variance sigma2 over N samples, E{||e||^2} = N*sigma2.
        if opts.residual_stopping
            opts.residual_threshold = sqrt(numel(y) * sigma2);
            log_fcn(sprintf(['algorithms stop when ||r|| reaches %.4g ' ...
                '(residual stopping)'], opts.residual_threshold));
        end
    elseif opts.residual_stopping
        log_fcn(['residual stopping requested but the measurement is ' ...
            'noiseless, so the sparsity cap governs']);
    end

    %% image formation
    x_hat = struct();

    if opts.execute_omp
        progress_fcn(0.45, 'running OMP');
        g = grids.(span_key(amb_if_alg.omp));

        x_hat_omp = omp_vec(y, g.A, Ks_latent, opts);
        x_hat.omp.image = reshape(x_hat_omp, Ny, g.Nx);

        % each algorithm carries the crossrange axis of the grid it was
        % actually given, since those axes now differ between algorithms
        x_hat.omp.x_array  = g.x_array;
        x_hat.omp.n_amb_if = g.n_amb_if;
        x_hat.omp.is_latent = latent_by_span.(span_key(g.n_amb_if));

        x_hat.omp.positions = extract_target_positions( ...
            x_hat.omp.image, g.x_array, y_array, Ks_latent, 'none');

        [x_hat.omp.error, x_hat.omp.pairs, x_hat.omp.missed, ...
            x_hat.omp.false_alarms, x_hat.omp.d] = ...
            calculate_reconstruction_error(target_locations, x_hat.omp.positions, miss_penalty);
    end

    if opts.execute_mod_omp
        progress_fcn(0.60, 'running modified OMP');

        % mod-OMP keeps the one-block dictionary and recovers each atom's
        % ambiguity index by testing its doppelgangers, refining the winner
        % off the grid. It therefore reports continuous positions in the same
        % [alpha_hat, p_hat] form as PROMP, not a gridded latent image.
        n_amb = n_amb_scat;
        g = grids.(span_key(amb_if_alg.mod_omp));

        grid_s = struct( ...
            'xk',                    g.xk, ...
            'yk',                    g.yk, ...
            'Wx',                    Wx, ...
            'Nx',                    g.Nx, ...
            'Ny',                    Ny, ...
            'x_array',               g.x_array, ...
            'y_array',               y_array, ...
            'cross_range_pixel_res', cross_range_pixel_res);

        % mod_omp_vec takes the scenario struct separately from the options:
        % it chooses between refining the doppelganger by a local offset
        % search and refining it by Newton's method, and the latter needs a
        % step count. Rs is that step count -- the same Newton/Gauss-Newton
        % budget NOMP and PROMP are given, so the comparison stays even.
        sc_mod = struct( ...
            'optimization_method',     opts.optimization_method, ...
            'num_optimization_steps',  opts.Rs, ...
            'num_offsets_pixels',      opts.num_offsets_pixels);

        [alpha_hat, p_hat] = mod_omp_vec(y, g.A, Ks_latent, grid_s, u0, ...
            theta, f_hat_l, fc, n_amb, sc_mod, opts);

        x_hat.mod_omp.positions = p_hat.';
        x_hat.mod_omp.alpha     = alpha_hat;
        x_hat.mod_omp.x_array   = g.x_array;
        x_hat.mod_omp.n_amb_if  = g.n_amb_if;
        x_hat.mod_omp.is_latent = latent_by_span.(span_key(g.n_amb_if));

        % carried for the plots: the unambiguous crossrange band
        x_hat.mod_omp.Wx        = Wx;

        [x_hat.mod_omp.error, x_hat.mod_omp.pairs, x_hat.mod_omp.missed, ...
            x_hat.mod_omp.false_alarms, x_hat.mod_omp.d] = ...
            calculate_reconstruction_error(target_locations, ...
            x_hat.mod_omp.positions, miss_penalty);
    end

    if opts.execute_bp
        progress_fcn(0.70, 'running backprojection');
        g = grids.(span_key(amb_if_alg.bp));

        x_hat_bp = g.A' * y;
        x_hat_bp = reshape(x_hat_bp, Ny, g.Nx);
        x_hat.bp.image = x_hat_bp / norm(x_hat_bp, 'fro');

        x_hat.bp.x_array  = g.x_array;
        x_hat.bp.n_amb_if = g.n_amb_if;
        x_hat.bp.is_latent = latent_by_span.(span_key(g.n_amb_if));

        if is_off_grid
            interpolation_type = 'linear';
        else
            interpolation_type = 'none';
        end

        x_hat.bp.positions = extract_target_positions( ...
            x_hat.bp.image, g.x_array, y_array, Ks_latent, interpolation_type);

        [x_hat.bp.error, x_hat.bp.pairs, x_hat.bp.missed, ...
            x_hat.bp.false_alarms, x_hat.bp.d] = ...
            calculate_reconstruction_error(target_locations, x_hat.bp.positions, miss_penalty);
    end

    if opts.execute_nomp
        progress_fcn(0.78, 'running NOMP');
        g = grids.(span_key(amb_if_alg.nomp));

        [alpha_hat, p_hat, p_hat_hist] = nomp_vec(y, g.A, opts.Rs, opts.Rc, ...
            Ks_latent, g.xk, g.yk, u0, theta, f_hat_l, fc, opts);

        x_hat.nomp.positions = p_hat.';
        x_hat.nomp.alpha = alpha_hat;
        x_hat.nomp.x_array  = g.x_array;
        x_hat.nomp.n_amb_if = g.n_amb_if;
        x_hat.nomp.is_latent = latent_by_span.(span_key(g.n_amb_if));

        % per-step position estimates for the traced atom, [nsteps x 2].
        % empty unless opts.save_histories is set
        x_hat.nomp.p_hat_hist = p_hat_hist.';

        [x_hat.nomp.error, x_hat.nomp.pairs, x_hat.nomp.missed, ...
            x_hat.nomp.false_alarms, x_hat.nomp.d] = ...
            calculate_reconstruction_error(target_locations, x_hat.nomp.positions, miss_penalty);
    end

    if opts.execute_promp
        progress_fcn(0.88, 'running PROMP');
        g = grids.(span_key(amb_if_alg.promp));

        [alpha_hat, p_hat, p_hat_hist] = promp_vec(y, g.A, opts.Rs, Ks_latent, ...
            g.xk, g.yk, u0, theta, f_hat_l, fc, opts);

        x_hat.promp.positions = p_hat.';
        x_hat.promp.alpha = alpha_hat;
        x_hat.promp.x_array  = g.x_array;
        x_hat.promp.n_amb_if = g.n_amb_if;
        x_hat.promp.is_latent = latent_by_span.(span_key(g.n_amb_if));

        % per-step position estimates for the traced atom, [nsteps x 2].
        % empty unless opts.save_histories is set
        x_hat.promp.p_hat_hist = p_hat_hist.';

        [x_hat.promp.error, x_hat.promp.pairs, x_hat.promp.missed, ...
            x_hat.promp.false_alarms, x_hat.promp.d] = ...
            calculate_reconstruction_error(target_locations, x_hat.promp.positions, miss_penalty);
    end

    progress_fcn(0.98, 'collecting results');

    % summarize the reconstruction error for each algorithm that ran
    log_fcn(' ');
    log_fcn('  algorithm   RMS position error [m]');
    for alg = ["omp" "mod_omp" "nomp" "promp" "bp"]
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
    out.y                     = y;
    out.A                     = A;     % reconstruction dictionary, [ML x Nx*Ny]
    out.as                    = as;    % measurement sensing matrix, [ML x Ks]

    out.amb_if_alg            = amb_if_alg;
    out.grids                 = rmfield_all_A(grids);

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
        'DopplerAliasing',      string(doppler_aliasing), ...
        'RangeCellMode',        string(range_cell_mode));

    progress_fcn(1, 'done');
end

function key = span_key(v)
% SPAN_KEY  Struct field name for the grid of a given ambiguity span.

    key = sprintf('span_%d', v);
end

function g = rmfield_all_A(grids)
% RMFIELD_ALL_A  The per-span grid geometry without the dictionaries.
%
%   The dictionaries are large and the caller already gets the scenario's own
%   as out.A, so the reported grids carry geometry only.

    g = grids;
    for fn = fieldnames(g).'
        g.(fn{1}) = rmfield(g.(fn{1}), 'A');
    end
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
