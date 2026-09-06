function opts = isar_test_options()
% ISAR_TEST_OPTIONS  Frozen simulation settings for the regression sweep.
%
%   opts = ISAR_TEST_OPTIONS() returns the complete set of options that
%   TEST_SCENARIO_REGRESSION runs with. Every field ISAR_RUN_SCENARIO reads is
%   named here, so nothing is inherited from ISAR_DEFAULT_OPTIONS.
%
%   That separation is the point. ISAR_DEFAULT_OPTIONS follows whatever
%   scenario the study is currently interested in, and it moved from X-band to
%   Ka-band (fc 1 -> 30 GHz, prf 200 -> 6000 Hz, w1 1e3 -> 170) when the GUI
%   defaults were aligned with isar_testing_v3.m. Because the test inherited
%   those fields, that change moved Wx and the crossrange resolution, took Nx
%   from 17 to 15, and failed all five scenarios -- reporting a numeric
%   difference rather than the settings change that caused it.
%
%   The values below are the ones in force when the baseline was recorded, so
%   the sweep answers "does the code still behave as it did", not "have the
%   defaults moved". Change them only when you intend to retire the old
%   baseline, and regenerate in the same commit.
%
%   Note that these are deliberately NOT physically representative: at
%   fc = 1 GHz with fs = 300 MHz the fractional bandwidth is 0.3, which smears
%   the Doppler ambiguity rather than aliasing it, and the range grid spans
%   more than the unambiguous range window. That is fine here -- a regression
%   test detects change, it does not have to be a valid experiment. Use
%   isar_testing_v3.m for physics.
%
%   See also TEST_SCENARIO_REGRESSION, ISAR_DEFAULT_OPTIONS, ISAR_RUN_SCENARIO.

    opts = struct( ...
        ... % --- which image formation algorithms to run ---
        'execute_omp',          true, ...
        'execute_nomp',         true, ...
        'execute_promp',        true, ...
        ... % mod-OMP postdates the baseline. Turning it on adds mod_omp_*
        ... % fields to every fingerprint, so enable it and regenerate
        ... % together, in a commit that says so.
        'execute_mod_omp',      false, ...
        'execute_bp',           true, ...
        ... % no per-algorithm override: every algorithm uses the
        ... % scenario's own ambiguity span, as when the baseline was recorded
        'amb_in_image_former',  struct(), ...
        ... % --- radar parameters ---
        'Nd',                   16, ...             % phase-history dimension
        'fc',                   1e9, ...            % [Hz] center frequency
        'prf',                  200, ...            % [Hz] pulse repetition frequency
        'fs',                   300e6, ...          % [Hz] sampling frequency
        'u0',                   1000, ...           % [m] range to the center of rotation
        ... % --- target motion ---
        'w0',                   pi, ...             % [rad/s] yawing rate
        'w1',                   1e3, ...            % [rad/s/s] yawing acceleration
        'w2',                   1e3, ...            % [rad/s/s/s] yawing jerk
        'jerk_mag',             5*pi, ...           % maneuvering trajectory jerk
        'complex_maneuver',     true, ...           % complex trajectory when Accelerating
        ... % --- geometry / grid ---
        'Ks',                   4, ...              % number of point scatterers
        ... % the baseline was recorded with unit amplitudes, before
        ... % target_magnitude existed; pinning 1 keeps those numbers
        'target_magnitude',     1, ...              % reflectivity magnitude of every scatterer
        'N_critical',           15, ...             % critically sampled range pixels
        'oversampling_factor',  4, ...              % grid oversampling when 'Oversampled'
        ... % --- solver settings ---
        'Rs',                   4, ...              % Newton steps per atom selection
        'Rc',                   2, ...              % cyclic refinements (NOMP)
        'use_range_approx',     false, ...          % exact range geometry
        ... % inert while execute_mod_omp is false; 'offset' is the path
        ... % that existed when the baseline was recorded
        'optimization_method',  'offset', ...       % mod-OMP doppelganger refinement
        'num_offsets_pixels',   10, ...             % 'newtons_and_offset' only
        'amb_refine_pixels',    5, ...              % mod-OMP offset-search half-width
        'residual_stopping',    false, ...          % fixed atom count, as when the baseline was recorded
        'residual_threshold',   [], ...             % unused unless residual_stopping
        'range_cell_mode',      false, ...          % full 2-D image, not one range cell
        'range_cell',           0, ...              % [m] unused unless range_cell_mode
        'save_histories',       true, ...           % record per-step position estimates
        'save_atom_idx',        1, ...              % which atom to trace
        ... % --- bookkeeping ---
        'compute_rank',         false, ...          % rank(A) is an SVD; off by default
        'seed',                 0);                 % fixed, so the sweep is deterministic
end
