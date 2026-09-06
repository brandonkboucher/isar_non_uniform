
% The goal of this script is to demonstrate reconstruction 
% with Doppler ambiguity with and without non-uniform 
% sampling. Reconstruction using compressed sensing image 
% formation techniques should fail for uniform sampling as 
% the dictionary atoms will be coherent. By contrast, the 
% dictionary atoms of the sensing matrix with non-uniform 
% sampling should be linearly independent, reducing mutual 
% coherence and improving the probability of meeting the 
% restrictive isometry property and the probability of a 
% successful reconstruction.

clear
clc
rng(0)

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
options.execute_mod_omp_with_promp  = false;

options.save_results                = true;
options.save_plots                  = true;
options.save_histories              = true; % only for one scatterer
options.save_atom_idx               = 1; % only if save_histories and must be less than sc.num_of_scatterers

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

%% scenario parameters

% define target parameters
sc.num_of_scatterers            = 10;

sc.is_target_accelerating       = true;
sc.is_target_maneuvering        = false;
sc.is_grid_oversampled          = true;
sc.is_closely_spaced            = false;
sc.is_off_grid                  = false;
sc.num_amb_having_scatterers    = 1;
sc.num_amb_in_image_former      = 1;

sc.yaw_acceleration             = 170;  % [rad/s/s]
sc.yaw_jerk                     = 0;    % [rad/s/s/s]
sc.target_magnitude             = 5;

% SNR in dB for the additive white complex Gaussian noise, or [] for a
% noiseless measurement. Nguyen et al. stop every algorithm when the signal
% residual reaches the noise level (Sec IV-A), so that criterion needs a noise
% level to exist: leave this empty and the pursuits run to the sparsity cap
% instead, which is what they have always done here.
sc.snr_db                       = [];

sc.num_optimization_steps   = 4;
sc.num_optimization_cycles  = 2; % only NOMP
sc.oversampling_factor      = 4;
sc.optimization_method      = 'newtons_and_offset'; % only mod-OMP, 'newtons' or 'offset'
sc.num_offsets_pixels       = 10; % only mod-OMP 'newtons_and_offset'

%% radar parameters

% instantiate constants
const = Constants;
c = const.c;

% define the dimensionality of the phase-history
Nd = 16;

fc      = 30 * const.GHz2Hz; % [Hz] center frequency - Ka-band
B       = 149.9 * const.MHz2Hz; % [Hz] bandwidth, not used
prf     = 6000; % [Hz] pulse repetition frequency
fs      = 300 * const.MHz2Hz; % [Hz] sampling frequency

lambda  = c / fc; % [m] wavelength

Tp = (1/fs) * Nd; % [s] pulse width
T = (1/prf) * Nd; % [s] simulation duration

t_m = (0:(1/prf):(T - 1/prf)).'; % [s] slow-time
t_hat = (0:(1/fs):(Tp - 1/fs)).'; % [s] fast-time

M = size(t_m,1); % number of pulses
L = size(t_hat,1); % number of fast-time samples

% define the range-frequency using the 
df_l = (fs/L);
f_hat_l = (-L/2)*df_l:df_l:(L/2 - 1)*df_l;

range_array = t_hat .* const.c / 2;

% define the latent image grid dimensions
% Range is bounded by its own ambiguity, exactly as crossrange is bounded by
% Wx. The range-frequency samples are spaced df = fs/L, which makes the
% unambiguous range window c/(2*df) -- and that window holds L-1 range
% resolution cells regardless of fs. A grid taller than this folds in range:
% at N_critical = 41 the grid spanned 2.7 windows, so scatterers aliased in
% range (coherence 0.998) as well as in crossrange.
N_critical = 15; % range cells, must be <= L = size(t_hat,1)

% create target scatterers and grid
[sc,target_locations, grid, theta_m, ...
    u0, is_doppler_aliasing] ...
    = create_target_and_grid(...
        sc, ...         % scenario parameters
        t_m, ...        % [s] (M x 1) slow-time
        fc, ...         % [Hz] center frequency
        prf, ...        % [Hz] pulse repetition frequency
        f_hat_l, ...    % [Hz] (L x 1) range-frequency
        N_critical ...  % dimension in x and y (for critically sampled)
        );
K = size(grid.xk,1);

% A scatterer an algorithm never reports is charged the crossrange width of
% the imaged scene, so declining to report a hard target cannot flatter the
% RMS the way it did for PROMP (one estimate for two scatterers scored
% 0.0034 m). See calculate_reconstruction_error.
miss_penalty = max(grid.x_array) - min(grid.x_array);

% compute sensing dictionary
as = zeros(M*L,sc.num_of_scatterers);
for k = 1:sc.num_of_scatterers

    ak = compute_atom(...
        target_locations(k,1), ...
        target_locations(k,2), ...
        u0, ...
        theta_m, ...
        f_hat_l, ...
        fc, ...
        options.use_range_approx);

    as(:,k) = ak;

    progress_bar('Sensing matrix', k, sc.num_of_scatterers)
end
fprintf('\n')

% compute dictionary for reconstruction
A = zeros(M*L,K);
for k = 1:K
    
    ak = compute_atom(...
        grid.xk(k), ...
        grid.yk(k), ...
        u0, ...
        theta_m, ...
        f_hat_l, ...
        fc, ...
        options.use_range_approx);

    % the reconstruction dictionary (grid atoms) is A; the off-grid scatterer
    % atoms `as` are used only to synthesize the measurement below.
    A(:,k) = ak;

    if mod(k,round(K/20)) == 0
        progress_bar('A matrix', k, K)
    end
end

% print matrix dimensions and conditioning
fprintf('\n')
fprintf(['A has dimension ', num2str(size(A,1)), ' by ', num2str(size(A,2)), '\n'])
fprintf(['A has a rank of ', num2str(rank(A)), '\n'])

% calculate the measurement: superpose the exact phase histories of the
% off-grid scatterers (each column of `as` is one scatterer's response).
alpha_s = sc.target_magnitude * ones(sc.num_of_scatterers,1);       % complex scattering amplitudes (unit for now)
y = as * alpha_s;           % = sum_k alpha_s(k) * as(:,k)

% additive white complex Gaussian noise, and the residual level the pursuits
% stop at. For per-sample variance sigma2 over N samples, E{||e||^2} =
% N*sigma2, so the noise level the paper refers to is sqrt(N*sigma2).
options.residual_threshold = [];
if isfield(sc, 'snr_db') && ~isempty(sc.snr_db) && isfinite(sc.snr_db)

    sigma2 = mean(abs(y).^2) / 10^(sc.snr_db/10);
    y = y + sqrt(sigma2/2) * (randn(size(y)) + 1j*randn(size(y)));

    options.residual_threshold = sqrt(numel(y) * sigma2);
    fprintf(['noise at %.1f dB SNR; algorithms stop when ||r|| reaches ' ...
        '%.4g (Nguyen et al. Sec IV-A)\n'], sc.snr_db, options.residual_threshold);
end
Y = reshape(y, M, L);

%% output
if options.calculate_mutual_coherence
    mu_mat = calculate_mutual_coherence(...
        A,...
        grid.x_array,...
        grid.y_array,...
        target_locations);
end

% modified OMP: each selected atom is compared against its ambiguous
% doppelgangers and refined off the grid, so it reports continuous scatterer
% positions in the same [alpha_hat, p_hat] form as PROMP
if isfield(options, 'execute_mod_omp') ...
    && options.execute_mod_omp

    % number of ambiguities the doppelganger test searches over
    n_amb = sc.num_amb_having_scatterers;

    % mod-OMP refines each selected atom off the grid while it tests the
    % ambiguous doppelgangers, so it returns continuous positions like PROMP
    % rather than a gridded latent image
    [alpha_hat, p_hat] = mod_omp_vec(...
        y,... % measurement
        A,... % gridded sensing matrix
        sc.num_of_latent_scatterers, ... % sparsity
        grid, ...
        u0, ...
        theta_m, ...
        f_hat_l, ...
        fc, ...
        n_amb, ...
        sc, ...
        options ...
        );

    x_hat.mod_omp.positions = p_hat.';
    x_hat.mod_omp.alpha     = alpha_hat;

    % carried for the plots: the unambiguous crossrange band, so a recovered
    % ambiguity can be read off the position panel
    x_hat.mod_omp.Wx        = grid.Wx;

    [x_hat.mod_omp.error, x_hat.mod_omp.pairs, ...
        x_hat.mod_omp.missed, x_hat.mod_omp.false_alarms, ...
        x_hat.mod_omp.d] = ...
        calculate_reconstruction_error(...
        target_locations, ...
        x_hat.mod_omp.positions, ...
        miss_penalty);

end

if isfield(options, 'execute_omp') ...
    && options.execute_omp

    x_hat_omp = omp_vec(...
        y,...
        A,...
        sc.num_of_latent_scatterers, ...
        options);

    x_hat.omp.image = reshape(x_hat_omp,grid.Ny,grid.Nx);

    x_hat.omp.positions =...
        extract_target_positions(...
        x_hat.omp.image, ...
        grid.x_array, ...
        grid.y_array, ...
        sc.num_of_scatterers, ...
        'none');

    [x_hat.omp.error, x_hat.omp.pairs, ...
        x_hat.omp.missed, x_hat.omp.false_alarms, ...
        x_hat.omp.d] = ...
        calculate_reconstruction_error(...
        target_locations, ...
        x_hat.omp.positions, ...
        miss_penalty);

end

% filtered backprojection using pseudo-inverse
if isfield(options, 'execute_bp') ...
    && options.execute_bp

    x_hat_bp = A'*y;
    x_hat_bp = reshape(x_hat_bp,grid.Ny,grid.Nx);
    % x_hat.bp.image = x_hat_bp / norm(x_hat_bp, "fro");
    x_hat.bp.image = x_hat_bp;

    interpolation_type = 'linear';
    if ~sc.is_off_grid
        interpolation_type = 'none';
    end

    x_hat.bp.positions = extract_target_positions(...
        x_hat.bp.image, ...
        grid.x_array, ...
        grid.y_array, ...
        sc.num_of_scatterers, ...
        interpolation_type);

    [x_hat.bp.error, x_hat.bp.pairs, ...
        x_hat.bp.missed, x_hat.bp.false_alarms, ...
        x_hat.bp.d] = ...
        calculate_reconstruction_error(...
        target_locations, ...
        x_hat.bp.positions, ...
        miss_penalty);
end

if isfield(options, 'execute_promp') ...
    && options.execute_promp

    [alpha_hat, p_hat] = promp_vec(...
        y, ...                          % measurement [ML x 1]
        A, ...                          % sensing matrix [ML x K]
        sc.num_optimization_steps, ...  % number of Gauss-Newton steps per atom selection
        sc.num_of_scatterers, ...       % sparsity
        grid.xk, ...                    % x position for each k-index [K]
        grid.yk, ...                    % y position for each k-index [K]
        u0, ...                         % center of rotation
        theta_m, ...                    % yaw angle as a function of time [M]
        f_hat_l, ...                    % range-frequencies [L]
        fc, ...                          % center frequency
        options);

    x_hat.promp.positions = p_hat.';
    x_hat.promp.alpha = alpha_hat;

    [x_hat.promp.error, x_hat.promp.pairs, ...
        x_hat.promp.missed, x_hat.promp.false_alarms, ...
        x_hat.promp.d] = ...
        calculate_reconstruction_error(...
        target_locations, ...
        x_hat.promp.positions, ...
        miss_penalty);

end

if isfield(options, 'execute_nomp') ...
        && options.execute_nomp

    [alpha_hat, p_hat, p_hat_hist] = nomp_vec(...
        y, ...
        A, ...
        sc.num_optimization_steps, ...
        sc.num_optimization_cycles, ...
        sc.num_of_scatterers, ...
        grid.xk, ...
        grid.yk, ...
        u0, ...
        theta_m, ...
        f_hat_l, ...
        fc, ...
        options);

    x_hat.nomp.positions = p_hat.';
    x_hat.nomp.alpha = alpha_hat;
    x_hat.nomp.p_hat_hist = p_hat_hist.';

    [x_hat.nomp.error, x_hat.nomp.pairs, ...
        x_hat.nomp.missed, x_hat.nomp.false_alarms, ...
        x_hat.nomp.d] = ...
        calculate_reconstruction_error(...
        target_locations, ...
        x_hat.nomp.positions, ...
        miss_penalty);

end

% summarize the reconstruction error for each algorithm that ran
fprintf('\n  algorithm   RMS position error [m]\n');
for alg = ["omp" "mod_omp" "nomp" "promp" "bp"]
    if isfield(x_hat, alg) && isfield(x_hat.(alg), 'error')
        fprintf('  %-10s  %.4f\n', alg, x_hat.(alg).error);
    end
end
fprintf('\n');

if options.save_plots

    % plotting
    plot_reconstruction(...
        1, ...
        target_locations,...
        u0,...
        grid.x_array, ...
        grid.y_array, ...
        x_hat,...
        options, ...
        'Scenario #1');

    plot_all_pairings(...
        1, ...
        target_locations,...
        u0,...
        grid.x_array, ...
        grid.y_array, ...
        x_hat,...
        options, ...
        'Scenario #1')
end
