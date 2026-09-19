function [sc,options,radar] = create_scenario_folded_collision()
% CREATE_SCENARIO_FOLDED_COLLISION  A draw where mod-OMP assigns the wrong
% ambiguity and OMP/NOMP do not.
%
%   A copy of CREATE_SCENARIO with one field changed: options.seed is pinned
%   to seed 188 of a 500-draw Monte Carlo sweep.
%
%   Truth (crossrange, range, band, folded crossrange):
%
%       -8.831  -1.625   -1   +0.107
%       +2.589  +2.470    0   +2.589
%       +9.480  -1.427   +1   +0.542
%
%   Why. The two outer scatterers are in different ambiguities but fold to
%   points 0.435 m apart in crossrange and 0.198 m in range, closer than one
%   resolution cell (0.59 m x 0.50 m). In the one-ambiguity image former the
%   two ghosts merge into a single peak. mod-OMP places one band-0 atom at
%   crossrange +0.34 m, between the two folded positions and matching
%   neither scatterer, and the leftover energy is then fitted wrongly. OMP and
%   NOMP have both ambiguities in their dictionary, so they see two separate
%   atoms rather than one merged peak.
%
%   Measured with isar_monte_carlo.m (mod-OMP on one ambiguity, baselines on
%   three), RMS position error:
%
%       mod-OMP  5.28 m     OMP  0.122 m     NOMP  0.081 m
%
%   (with the projected-Hessian Newton fix: mod-OMP 5.29 m, NOMP 0.053 m)
%
%   Causal check: moving the band -1 scatterer 2 m in range (breaking the
%   folded collision) and changing nothing else drops mod-OMP to ~0.01-0.07 m.
%   The failure persists with the projected-Hessian Newton fix, so it is not
%   the range-refinement bug that caused most of the earlier gross errors.
%
%   Scope. This is the mod-OMP-specific failure that is left once Newton
%   refines range correctly: about 3 draws in 500. A folded collision at a
%   *single* range is hard for NOMP too, so most collisions fail for both.
%   It is not a Doppler-width effect: raising the yaw acceleration widens the
%   Wx_m spread and lowers ghost coherence, and collisions then fail less.
%
%   Run it through isar_monte_carlo.m or isar_testing_v3.m; both give mod-OMP
%   one ambiguity and OMP/NOMP three.
%
%   See also CREATE_SCENARIO, ISAR_MONTE_CARLO.


    % the radar block below is written in terms of the physical constants and
    % of radar's own fields, so both have to exist before it
    const = Constants;
    c = const.c;

    %% options
    options.calculate_mutual_coherence  = false;
    options.log_scale_plotting          = true;
    
    % backprojection, sbl and itsa
    options.execute_bp                  = true;
    options.execute_itsa                = false;
    options.execute_sbl                 = false;
    
    % orthogonal matching pursuit and related algorithms
    options.execute_omp                 = true;
    options.execute_nomp                = true;
    options.execute_nomp_newton         = true; % NOMP with newton_method_exact
    options.execute_promp               = false;
    options.execute_mod_omp             = true;
    options.execute_mod_omp_newton      = true; % mod-OMP with newton_method_exact + cyclic refinement

    % candidate plot for mod_omp_newton: one figure per iteration that
    % detected an aliased scatterer's ghost, or every iteration if the second
    % flag is set
    options.plot_mod_omp_candidates         = true;
    options.plot_candidates_all_iterations  = true;
    options.plot_best_candidate_paths       = true; % highlight each ambiguity's best Newton trajectory

    options.save_results                = true;
    options.save_plots                  = true;
    options.save_histories              = false; % only for one scatterer
    options.save_atom_idx               = 1; % only if save_histories and must be less than sc.num_of_scatterers
    
    % THE ONLY DEPARTURE FROM create_scenario. Seed 188 of the Monte Carlo
    % sweep puts two scatterers in bands -1 and +1 whose folded crossranges
    % (+0.107 and +0.542 m) sit 0.48 m apart, inside one resolution cell.
    options.seed                        = 188;

    if options.save_plots
        delete('plots/*')
    end
    
    % use the linearized range model of Cheng et al. (2019) eq (1) throughout
    % (both the measurement atoms `as` and the reconstruction dictionary `a`);
    % false = exact range geometry. See compute_atom.m.
    options.use_range_approx             = false;
    
    % ---------------- testing ------------------
    options.manual_close_spacing        = false;
    options.manual_maneuvering_target   = false;
    options.debug_printing              = false;
    
    %% scenario parameters
    
    % define target parameters
    sc.num_of_scatterers            = 3;
    
    sc.is_target_accelerating       = true;
    sc.is_target_maneuvering        = false;
    sc.is_grid_oversampled          = true;
    sc.is_closely_spaced            = false;
    sc.is_off_grid                  = true;
    sc.num_amb_having_scatterers    = 3;
    sc.num_amb_in_image_former      = 1;
    sc.num_amb_in_image_former_baselines = 3; % OMP/NOMP dictionary span in isar_testing_v3
    
    sc.yaw_acceleration             = 170;  % [rad/s/s]
    sc.yaw_jerk                     = 0;    % [rad/s/s/s]
    sc.target_magnitude             = 5;
    
    % SNR in dB for the additive white complex Gaussian noise, or [] for a
    % noiseless measurement. Nguyen et al. stop every algorithm when the signal
    % residual reaches the noise level (Sec IV-A), so that criterion needs a noise
    % level to exist: leave this empty and the pursuits run to the sparsity cap
    % instead, which is what they have always done here.
    sc.snr_db                   = [];
    
    sc.num_optimization_steps   = 4;
    sc.num_optimization_cycles  = 2; % only NOMP
    sc.oversampling_factor      = 4;
    sc.optimization_method      = 'newtons_and_offset'; % only mod-OMP, 'newtons' or 'offset'
    sc.num_offsets_pixels       = 10; % only mod-OMP 'newtons_and_offset'

    % define the latent image grid dimensions
    % Range is bounded by its own ambiguity, exactly as crossrange is bounded by
    % Wx. The range-frequency samples are spaced df = fs/L, which makes the
    % unambiguous range window c/(2*df) -- and that window holds L-1 range
    % resolution cells regardless of fs. A grid taller than this folds in range:
    % at N_critical = 41 the grid spanned 2.7 windows, so scatterers aliased in
    % range (coherence 0.998) as well as in crossrange.
    sc.N_critical = 15; % range cells, must be <= L = size(t_hat,1)


    % define the dimensionality of the phase-history
    radar.Nd          = 16; % dimensionality of latent image
    radar.fc          = 30 * const.GHz2Hz; % [Hz] center frequency - Ka-band
    radar.B           = 149.9 * const.MHz2Hz; % [Hz] bandwidth, not used
    radar.prf         = 6000; % [Hz] pulse repetition frequency
    radar.fs          = 300 * const.MHz2Hz; % [Hz] sampling frequency
    radar.lambda      = c / radar.fc; % [m] wavelength
    radar.Tp          = (1/radar.fs) * radar.Nd; % [s] pulse width
    radar.T           = (1/radar.prf) * radar.Nd; % [s] simulation duration

end

