function [sc, options, radar] = paper_sweep_settings()
% PAPER_SWEEP_SETTINGS  Frozen scenario for the paper's Monte Carlo sweep.
%
%   [sc, options, radar] = PAPER_SWEEP_SETTINGS() returns the complete
%   scenario, options and radar structs that RUN_PAPER_SWEEP runs with. Every
%   field the sweep reads is named here, so nothing is inherited from
%   CREATE_SCENARIO.
%
%   That separation is the point, and it is the same one ISAR_TEST_OPTIONS
%   makes: create_scenario.m follows whatever the study is currently
%   interested in, and a field moving there (a flag, an acceleration, an
%   offset count) would change the published numbers without anyone saying so.
%   The values below are the ones in force when the baseline was recorded.
%
%   Change them only when you intend to retire the baseline, and regenerate it
%   in the same commit:  test_paper_sweep('regenerate')
%
%   See also RUN_PAPER_SWEEP, TEST_PAPER_SWEEP, CREATE_SCENARIO.

    const = Constants;
    c = const.c;

    %% options
    options.calculate_mutual_coherence  = false;
    options.log_scale_plotting          = false;
    options.execute_bp                  = false;
    options.execute_itsa                = false;
    options.execute_sbl                 = false;
    options.execute_omp                 = true;
    options.execute_nomp                = true;
    options.execute_promp               = false;
    options.execute_mod_omp             = true;
    options.execute_mod_omp_newton      = true; % mod-OMP with newton_method_exact + cyclic refinement
    options.execute_nomp_newton         = true;
    options.save_results                = false;
    options.save_plots                  = false;   % never delete plots/ from a test
    options.save_histories              = false;
    options.save_atom_idx               = 1;
    options.seed                        = 0;       % iteration i uses seed + i - 1
    options.use_range_approx            = false;   % exact range geometry
    options.manual_close_spacing        = false;
    options.manual_maneuvering_target   = false;
    options.debug_printing              = false;
    options.residual_threshold          = [];      % no noise, so no early stop

    %% scenario
    sc.num_of_scatterers            = 7;
    sc.is_target_accelerating       = true;
    sc.is_target_maneuvering        = false;
    sc.is_grid_oversampled          = true;
    sc.is_closely_spaced            = false;
    sc.is_off_grid                  = true;
    sc.num_amb_having_scatterers    = 7;   % scatterers spread over three bands
    sc.num_amb_in_image_former      = 1;   % mod-OMP's image former
    sc.num_amb_in_image_former_baselines = 7;  % OMP's and NOMP's dictionary
    sc.yaw_acceleration             = 170; % [rad/s/s]
    sc.yaw_jerk                     = 0;   % [rad/s/s/s]
    sc.target_magnitude             = 5;
    sc.snr_db                       = [];  % noiseless
    sc.num_optimization_steps       = 4;   % Rs
    sc.num_optimization_cycles      = 2;   % Rc
    sc.oversampling_factor          = 4;
    sc.optimization_method          = 'newtons_and_offset';
    sc.num_offsets_pixels           = 10;  % randomized variant's candidate count
    sc.N_critical                   = 15;  % range cells

    %% radar
    radar.Nd     = 16;
    radar.fc     = 30 * const.GHz2Hz;      % [Hz] Ka-band
    radar.B      = 149.9 * const.MHz2Hz;   % [Hz] not used
    radar.prf    = 6000;                   % [Hz]
    radar.fs     = 300 * const.MHz2Hz;     % [Hz]
    radar.lambda = c / radar.fc;
    radar.Tp     = (1/radar.fs) * radar.Nd;
    radar.T      = (1/radar.prf) * radar.Nd;
end
