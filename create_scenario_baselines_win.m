function [sc,options,radar] = create_scenario_baselines_win()
% CREATE_SCENARIO_BASELINES_WIN  A draw where OMP and NOMP beat mod-OMP.
%
%   A copy of CREATE_SCENARIO with one field changed: options.seed is pinned
%   to a draw observed during the Monte Carlo sweep in which both baselines
%   outperform mod-OMP by roughly two orders of magnitude.
%
%   Measured, with mod-OMP on one ambiguity and OMP/NOMP on three (the split
%   isar_monte_carlo.m makes):
%
%       mod-OMP  5.1495 m
%       OMP      0.0358 m
%       NOMP     0.0138 m
%
%   Why. The three scatterers sit at crossrange [-5.09, +0.45, +14.05] m with
%   Wx = 9.543 m, so band 0 covers [-4.746, +4.746]. Acceleration makes the
%   fold distance that actually applies the one set by the mean rotation rate
%   over the aperture,
%
%       Wx_eff = Wx * w0 / (w0 + w1*T/2) = 8.938 m,
%
%   which is 0.605 m shorter than the nominal Wx. The scatterer at +14.054
%   therefore folds to +14.054 - 8.938 = +5.115, which is OUTSIDE band 0. Its
%   ghost is not in the image former at all, so mod-OMP's detection step
%   cannot anchor on it: it picks a sidelobe at -3.856 instead, and every
%   doppelganger candidate is then offset from that wrong anchor. It chose
%   ambiguity +1 and reported +5.136 against a true +14.054.
%
%   OMP and NOMP are handed all three ambiguities, so they address that
%   scatterer's atom directly and never face the question.
%
%   Nothing here is a tuning artefact: no neighbourhood size fixes it, because
%   the anchor is wrong rather than the offset from it. The scatterers in the
%   outer (Wx - Wx_eff)/Wx = 6.3% of each band have this property, so a draw
%   like this appears whenever one lands there.
%
%   Run it through isar_monte_carlo.m, which gives mod-OMP one ambiguity and
%   the baselines three. isar_testing_v3.m puts every algorithm on the same
%   span, so it will not show the comparison this scenario is built for.
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
    options.execute_promp               = true;
    options.execute_mod_omp             = true;

    options.save_results                = true;
    options.save_plots                  = true;
    options.save_histories              = false; % only for one scatterer
    options.save_atom_idx               = 1; % only if save_histories and must be less than sc.num_of_scatterers

    % THE ONLY DEPARTURE FROM create_scenario. Seed 4 of the Monte Carlo sweep
    % puts a scatterer at +14.05 m, whose ghost falls outside the single
    % ambiguity mod-OMP is given. Seeds 6, 8 and 9 are the same story with
    % scatterers at +12.61, +13.07 and -14.22 m.
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
