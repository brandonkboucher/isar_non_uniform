function opts = isar_default_options(user_opts)
% ISAR_DEFAULT_OPTIONS  Non-scenario settings for ISAR_RUN_SCENARIO.
%
%   opts = ISAR_DEFAULT_OPTIONS() returns the radar, solver and bookkeeping
%   settings used by isar_testing_v2.m. ISAR_DEFAULT_OPTIONS(user_opts)
%   overlays the fields of user_opts on top of those defaults, so a caller
%   only has to name what it wants to change.

    const = Constants;

    opts = struct( ...
        ... % --- which image formation algorithms to run ---
        'execute_omp',          true, ...
        'execute_nomp',         true, ...
        'execute_promp',        true, ...
        'execute_bp',           true, ...
        ... % --- radar parameters ---
        'Nd',                   16, ...                     % phase-history dimension
        'fc',                   1 * const.GHz2Hz, ...       % [Hz] center frequency
        'prf',                  200, ...                    % [Hz] pulse repetition frequency
        'fs',                   300 * const.MHz2Hz, ...     % [Hz] sampling frequency
        'u0',                   1000, ...                   % [m] range to the center of rotation
        ... % --- target motion ---
        'w0',                   pi, ...                     % [rad/s] yawing rate
        'w1',                   1e3, ...                    % [rad/s/s] yawing acceleration
        'w2',                   1e3, ...                    % [rad/s/s/s] yawing jerk
        'jerk_mag',             5*pi, ...                    % jerk magnitude of the maneuvering trajectory
        'complex_maneuver',     true, ...                   % use create_complex_target_trajectory when Accelerating
        ... % --- geometry / grid ---
        'Ks',                   10, ...                     % number of point scatterers
        'N_critical',           41, ...                     % critically sampled range pixels
        'oversampling_factor',  4, ...                      % grid oversampling when 'Oversampled'
        ... % --- solver settings ---
        'Rs',                   4, ...                      % Newton steps per atom selection (NOMP)
        'Rc',                   2, ...                      % cyclic refinements (NOMP)
        'use_range_approx',     false, ...                  % linearized range of Cheng et al. eq (1)
        ... % --- bookkeeping ---
        'compute_rank',         false, ...                  % rank(A) is an SVD; off by default
        'seed',                 0, ...                      % [] to leave the RNG alone
        'log_fcn',              @(s) fprintf('%s\n', s), ...
        'progress_fcn',         @(frac, msg) []);

    if nargin >= 1 && ~isempty(user_opts)
        if ~isstruct(user_opts)
            error('isar_default_options:type', 'user_opts must be a struct');
        end

        unknown = setdiff(fieldnames(user_opts), fieldnames(opts));
        if ~isempty(unknown)
            error('isar_default_options:unknownOption', ...
                'unknown option(s): %s', strjoin(unknown.', ', '));
        end

        for f = fieldnames(user_opts).'
            opts.(f{1}) = user_opts.(f{1});
        end
    end
end
