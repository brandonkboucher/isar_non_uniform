function [sc,options,radar] = create_scenario_baselines_win()
% CREATE_SCENARIO_BASELINES_WIN  A draw where exact-Newton mod-OMP picks the
% wrong ambiguity because Newton cannot reach the true peak.
%
%   A copy of CREATE_SCENARIO with one field changed: options.seed = 4.
%
%   Measured with isar_testing_v3.m (mod-OMP on one ambiguity, OMP/NOMP on
%   three), RMS position error:
%
%       mod-OMP, exact Newton + cyclic (mod_omp_newton)   5.146 m
%       mod-OMP (repo, mod_omp_vec)                       0.079 m
%       OMP                                               0.031 m
%       NOMP                                              0.014 m
%
%   Why. Scatterer 3 sits at (+13.16, -2.27) m in band +1. At iteration 3 the
%   one-ambiguity image former detects its band-0 ghost at (4.30, -2.13). The
%   residual correlates 0.999 with the true atom and 0.876 with the ghost,
%   but the band +1 candidates never reach the true peak:
%
%     - every candidate is seeded at the detection's range, which the ghost
%       has shifted 0.14 m (0.28 of a range cell) from the truth;
%     - the true peak is narrow, so that seed sits near its inflection point,
%       where Newton's quadratic model is poor: the first exact-Newton step
%       overshoots about 2x, lands lower on the far flank, and is rejected;
%     - with no backtracking the search stops, 19 of 21 band +1 candidates
%       never leave their seeds, and they are scored there (best 0.857).
%
%   The smeared band-0 ghost is broad, so Newton does climb it (0.876), and it
%   wins. The candidate plot in isar_testing_v3.m (plot_mod_omp_candidates)
%   shows all of this for iteration 3.
%
%   Representative. In a 500-draw sweep, 43 draws fail grossly for exact
%   mod-OMP + cyclic but not for NOMP; 35 of them are this same search
%   failure (a band +1 scatterer assigned to band 0 although the true atom
%   correlates better). See claude_scratch/classify_exact_failures.m.
%
%   History. Seed 4 was first pinned for a different failure (a ghost folding
%   outside band 0 because Wx used the initial rotation rate); computing Wx
%   from the mean rotation rate removed that.
%
%   See also CREATE_SCENARIO, CREATE_SCENARIO_FOLDED_COLLISION, ISAR_TESTING_V3.


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
    
    % THE ONLY DEPARTURE FROM create_scenario. Seed 4 puts a scatterer at
    % (+13.16, -2.27) m in band +1 whose true peak exact-Newton mod-OMP cannot
    % reach from the ghost's seeds.
    options.seed                        = 4;

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
    
    % DEPRECIATED, we use Wx for the crossrange offset
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

