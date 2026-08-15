
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

% extract scenarios
scenarios = readtable('documentation/scenarios.xlsx');
delete('plots/*')

options.calculate_mutual_coherence  = false;
options.execute_itsa                = false;
options.execute_omp                 = true;
options.execute_nomp                = true;
options.execute_sbl                 = false;
options.execute_bp                  = true;
options.log_scale_plotting          = false;

Rc = 1;
Rs = 1;

% use the linearized range model of Cheng et al. (2019) eq (1) throughout
% (both the measurement atoms `as` and the reconstruction dictionary `a`);
% false = exact range geometry. See compute_atom.m.
use_range_approx            = false;

%% radar parameters

% instantiate constants
const = Constants;
c = const.c;

% define the dimensionality of the phase-history
Nd = 16;

fc      = 1 * const.GHz2Hz; % [Hz] center frequency - X-band
B       = 149.9 * const.MHz2Hz; % [Hz] bandwidth, not used
prf     = 200; % [Hz] pulse repetition frequency, previously: 10 * const.kHz2Hz
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

range_resolution = const.c / (2 * (max(f_hat_l) - min(f_hat_l))); % [m]
range_array = t_hat .* const.c / 2;

% define the latent image grid dimensions
N_critical = 41; % dimension in x and y (for critically sampled)

% iterate through scenarios
num_scenarios = size(scenarios,1);

% Skip rows that differ only by image formation algorithm. The options struct
% above runs every algorithm on every scenario, so those rows would repeat the
% same simulation. Case #, algorithm and prediction are descriptive rather than
% parameters of the collection, so they are excluded from the comparison.
param_vars = setdiff(scenarios.Properties.VariableNames, ...
    {'Case_', 'ImageFormationAlgorithm', 'Prediction'}, 'stable');

% one key string per row; lower() so 'Off-Grid' and 'Off-grid' compare equal
scenario_keys = strings(num_scenarios,1);
for i = 1:num_scenarios
    fields = table2cell(scenarios(i, param_vars));
    scenario_keys(i) = lower(strjoin( ...
        cellfun(@(v) char(string(v)), fields, 'UniformOutput', false), '|'));
end

% keep the first row of each group of identical parameter sets
[~, isc_unique, scenario_group] = unique(scenario_keys, 'stable');
fprintf('running %d of %d scenarios (%d skipped as algorithm-only duplicates)\n\n', ...
    numel(isc_unique), num_scenarios, num_scenarios - numel(isc_unique));

for isc = 1:num_scenarios

    if ~any(isc == isc_unique)
        continue
    end

    covers = scenarios.Case_(scenario_group == scenario_group(isc));
    fprintf('--- case %s ---\n', strjoin(string(covers), ', '));

    %% target definition
    % Define the target's motion
    w0 = pi; % [rad/s] yawing rate
    if strcmp(scenarios.AngleRate(isc),'Accelerating')
        w1 = 1e3; % [rad/s/s] yawing acceleration
        w2 = 1e3;    % [rad/s/s] yawing jerk
    else
        w1 = 0; % [rad/s/s] yawing acceleration
        w2 = 0;    % [rad/s/s] yawing jerk
    end

    % determine the targets rotation as a function of slow time
    theta = w0 * t_m ...
        + (1/2) * w1 * t_m .* t_m ...
        + (1/3) * w2 * t_m .* t_m .* t_m;

    % determine the angular span which determines the cross
    % range resolution
    sin_span = max(sin(theta)) - min(sin(theta));
    cross_range_resolution = const.c ...
        / (2 * (fc + max(f_hat_l)) * sin_span);

    % define the crossrange unambiguous extent, this is the extent at which the 
    % crossrange is unambiguous. if a scatterer's crossrange exceeds this extent
    % then it is ambiguous
    W_x = c * prf / (2 * fc * w0);

    % determine the number of pixels within the crossrange unambiguous extent
    n_pix_per_amb = round(W_x / cross_range_resolution);

    % extract the number of ambiguous scatterers 
    n_amb_scat = scenarios.NumberOfAmbiguitiesHavingScatterers(isc);
    n_amb_if   = scenarios.NumberOfAmbiguitiesInImageFormer(isc);

    if strcmp(scenarios.ImageFormerGridDensity(isc), 'Critical')
        Kr = 1;
        Kc = 1;
    elseif strcmp(scenarios.ImageFormerGridDensity(isc), 'Oversampled')
        Kr = 2;
        Kc = 2;
    end

    % redefine the grid in order to handle the number of ambiguous scatterers
    Ny = Kr * N_critical;
    Nx = Kc * n_amb_if * n_pix_per_amb;

    % ensure the number of pixels is odd in both directions
    Nx = Nx + (mod(Nx,2)==0);
    Ny = Ny + (mod(Ny,2)==0);

    K = Nx * Ny; % number of vectorized samples

    % define the resolution of the grid
    range_pixel_res = range_resolution / Kr;
    cross_range_pixel_res = cross_range_resolution / Kc;
    
    % distance from the origin to the target center
    u0 = 1000; % [m]
    
    % define the range and crossrange grid of the latent image
    % centered on 0 so the critical grid is a subset of the oversampled grid
    y_array = ((0:Ny-1) - (Ny-1)/2) * range_pixel_res;
    x_array = ((0:Nx-1) - (Nx-1)/2) * cross_range_pixel_res;
    [X,Y] = meshgrid(x_array,y_array);
    
    % define the x and y component for each grid point
    yk = reshape(Y, [], 1);
    xk = reshape(X, [], 1);
   
    %% sensing matrix and latent image
    
    % number of point scatterers
    Ks = 10;
    
    % assign scatterers depending on the number of ambiguities. if the number
    % of ambiguities is 1 -> [0], if 3 -> [-1,0,1], etc
    ii        = 0:(n_amb_scat-1);
    amb_index = ceil(ii/2) .* (-1).^ii;

    % round-robin so every requested ambiguity holds at least one scatterer
    amb_of_k = amb_index(mod(0:Ks-1, n_amb_scat) + 1).';

    % define the target scatterer locations, in meters. Columns are
    % [crossrange, range] to match compute_atom's (x, y) argument order (as used
    % for xk, yk below). Crossrange is uniform within the scatterer's own
    % ambiguity; range is uniform over the grid extent (range has no ambiguity
    % structure here -- it is set by bandwidth, not by the PRF).
    target_locations = [ ...
        amb_of_k * W_x + (rand(Ks,1) - 0.5) * W_x, ...
        (rand(Ks,1) - 0.5) * (y_array(end) - y_array(1)) ];

    % On-Grid rows additionally snap onto the critically-sampled grid. The
    % ambiguity assignment above already happened, so on-grid scatterers are
    % distributed across ambiguities exactly like the off-grid ones.
    if ~strcmpi(scenarios.ScattererLocations(isc), 'Off-grid')
        target_locations(:,1) = round(target_locations(:,1)/cross_range_resolution) * cross_range_resolution;
        target_locations(:,2) = round(target_locations(:,2)/range_resolution)       * range_resolution;
    end
    
    % determine if Doppler aliasing will occur
    theta_dot = w0 + w1*t_m + w2*t_m.^2;
    fd_max = max((2*fc/c) * abs(target_locations(:,1)) * max(abs(theta_dot)));
    fprintf('W_x = %.3f m (%d critical pixels); scatterers in ambiguities %s (requested %d)\n', ...
        W_x, n_pix_per_amb, mat2str(unique(amb_of_k).'), n_amb_scat);
    fprintf('image former spans %.2f ambiguities (requested %d)\n', Nx*cross_range_pixel_res/W_x, n_amb_if);
    if fd_max > prf/2
        disp('Doppler aliasing will occur')
    else
        disp('No Doppler aliasing')
    end

    % compute sensing dictionary
    as = zeros(M*L,Ks);
    tt = tic;
    for k = 1:Ks
        
        if mod(k,round(k/20)) == 0
            ttt = toc(tt);
            fprintf(['k: ', num2str(k), ', ', num2str(ttt), ' seconds\n'])
            tt = tic;
        end
    
        ak = compute_atom(...
            target_locations(k,1), ...
            target_locations(k,2), ...
            u0, ...
            theta, ...
            f_hat_l, ...
            fc, ...
            use_range_approx);
    
        as(:,k) = ak;
    
    end
    
    % compute dictionary for reconstruction
    A = zeros(M*L,K);
    tt = tic;
    for k = 1:K
        
        if mod(k,round(k/20)) == 0
            ttt = toc(tt);
            fprintf(['k: ', num2str(k), ', ', num2str(ttt), ' seconds\n'])
            tt = tic;
        end
    
        ak = compute_atom(...
            xk(k), ...
            yk(k), ...
            u0, ...
            theta, ...
            f_hat_l, ...
            fc, ...
            use_range_approx);
    
        % the reconstruction dictionary (grid atoms) is A; the off-grid scatterer
        % atoms `as` are used only to synthesize the measurement below.
        A(:,k) = ak;
    end
    
    % print matrix dimensions and conditioning
    fprintf(['A has dimension ', num2str(size(A,1)), ' by ', num2str(size(A,2)), '\n'])
    fprintf(['A has a rank of ', num2str(rank(A)), '\n'])
    
    % calculate the measurement: superpose the exact phase histories of the
    % off-grid scatterers (each column of `as` is one scatterer's response).
    alpha_s = ones(Ks,1);       % complex scattering amplitudes (unit for now)
    y = as * alpha_s;           % = sum_k alpha_s(k) * as(:,k)
    Y = reshape(y, M, L);
    
    %% output
    if options.calculate_mutual_coherence
        mu_mat = calculate_mutual_coherence(...
            A,...
            x_array,...
            y_array,...
            target_locations);
    end

     % filtered backprojection using pseudo-inverse
    if isfield(options, 'execute_omp') ...
        && options.execute_omp

        x_hat_omp = omp_vec(...
            y,...
            A,...
            Ks+2);
        x_hat.omp.image = reshape(x_hat_omp,Ny,Nx);
        
        x_hat.omp.positions = extract_target_positions(...
            x_hat.omp.image, ...
            x_array, ...
            y_array, ...
            Ks, ...
            'none');

        x_hat.omp.error = calculate_reconstruction_error(...
            target_locations, ...
            x_hat.omp.positions);

    end


    % filtered backprojection using pseudo-inverse
    if isfield(options, 'execute_bp') ...
        && options.execute_bp

        x_hat_bp = A'*y;
        x_hat_bp = reshape(x_hat_bp,Ny,Nx);
        x_hat.bp.image = x_hat_bp / norm(x_hat_bp, "fro");

        if ~strcmpi(scenarios.ScattererLocations(isc), 'Off-grid')
            interpolation_type = 'none';
        else
            interpolation_type = 'linear';
        end

        x_hat.bp.positions = extract_target_positions(...
            x_hat.bp.image, ...
            x_array, ...
            y_array, ...
            Ks, ...
            interpolation_type);

        x_hat.bp.error = calculate_reconstruction_error(...
            target_locations, ...
            x_hat.bp.positions);
    end

    if isfield(options, 'execute_nomp') ...
            && options.execute_nomp

        [alpha_hat, p_hat] = nomp_vec(...
            y, ...
            A, ...
            Rs, ...
            Rc, ...
            Ks, ...
            xk, ...
            yk, ...
            u0, ...
            theta, ...
            f_hat_l, ...
            fc);

        x_hat.nomp.positions = p_hat;
        x_hat.nomp.alpha = alpha_hat;

        x_hat.nomp.error = calculate_reconstruction_error(...
            target_locations, ...
            x_hat.nomp.positions);

    end

    % plotting
    plot_reconstruction(...
        scenarios.Case_(isc), ...
        target_locations,...
        u0,...
        x_array, ...
        y_array, ...
        x_hat,...
        options, ...
        ['Scenario #', num2str(scenarios.Case_(isc))]);
    
end