
% Monte Carlo characterisation of the ambiguity-refined OMP.
%
% Every iteration draws a new set of scatterers spread over
% sc.num_amb_having_scatterers ambiguities and forms one measurement that all
% four algorithms see. What differs is the dictionary each is handed:
%
%   mod-OMP             one ambiguity in the image former
%   OMP, NOMP, PROMP    all sc.num_amb_having_scatterers ambiguities
%
% That is the comparison the computational argument rests on: mod-OMP recovers
% a scatterer lying outside the imaged band by testing its ambiguous
% doppelgangers, so it never needs the wider dictionary the others are given.
%
% Both grids, and both dictionaries, are fixed by the radar parameters rather
% than by the random draw, so they are built once before the loop. Only the
% scatterer positions and the measurement change from iteration to iteration.
%
% This script drives the existing solvers unchanged; nothing here alters how
% any of them behave.

clear
clc

%% monte carlo settings

n_iterations = 1;

results_file = fullfile('results', 'monte_carlo.mat');
figure_file  = fullfile('plots',   'monte_carlo_cdf.png');

%% scenario, options and radar parameters
% Taken from create_scenario so this sweep and isar_testing_v3.m cannot drift
% apart. Only what the sweep itself needs is overridden below; the physics
% (scatterer count, on/off grid, motion, band) belongs in create_scenario.

[sc, options, radar] = create_scenario_baselines_win();
base_seed    = options.seed;        % iteration i uses seed base_seed + i - 1

% --- monte carlo overrides ---------------------------------------------

% the four algorithms compared here. Backprojection forms an image rather than
% a scatterer list, so it is not part of this comparison.
options.execute_bp      = false;

% nothing is plotted or written per iteration; the sweep writes one summary,
% one .mat and one figure at the end
options.save_plots      = false;
options.save_results    = false;
options.save_histories  = false;   % per-step traces suit one run, not 500

% create_target_and_grid gates its reporting, and the ghost-coherence scan
% that goes with it, on this. Off for the sweep: that scan costs an atom build
% per scan point per scatterer on every iteration.
options.debug_printing  = false;

% save a position figure every this many iterations, 0 for none. The sweep
% reports distributions; these are for looking at what an individual run
% actually did, which is the only way to see an ambiguity assigned wrongly.
options.plot_every              = 1;
positions_dir                   = fullfile('plots', 'monte_carlo');

% the two image former spans being compared. num_amb_in_image_former is set
% per arm below rather than once for the scenario.
amb_if_mod_omp                  = 1;
amb_if_baselines                = 3;
sc.num_amb_in_image_former      = amb_if_mod_omp;

%% derived radar quantities

t_m   = (0:(1/radar.prf):(radar.T  - 1/radar.prf)).'; % [s] slow-time
t_hat = (0:(1/radar.fs):(radar.Tp - 1/radar.fs)).';   % [s] fast-time

M = size(t_m,1);   % number of pulses
L = size(t_hat,1); % number of fast-time samples

df_l    = (radar.fs/L);
f_hat_l = (-L/2)*df_l:df_l:(L/2 - 1)*df_l;

%% grids and dictionaries
% Fixed by the radar parameters, so built once. The seed used here only
% affects the throwaway target draw that comes back with them.

fprintf('building the two image formers\n\n');

sc_mod  = sc; sc_mod.num_amb_in_image_former  = amb_if_mod_omp;
sc_base = sc; sc_base.num_amb_in_image_former = amb_if_baselines;

rng(base_seed)
[sc_mod, target_check_mod, grid_mod, theta_m, u0] = ...
    create_target_and_grid(sc_mod, t_m, radar, f_hat_l, options);

rng(base_seed)
[sc_base, target_check_base, grid_base] = ...
    create_target_and_grid(sc_base, t_m, radar, f_hat_l, options);

% The comparison is only meaningful if both arms see the same scatterers. The
% draw uses Wx and the range extent, neither of which depends on the image
% former span, so the same seed has to reproduce it -- check rather than trust.
if ~isequal(target_check_mod, target_check_base)
    error('isar_monte_carlo:targetMismatch', ...
        ['the two image former spans produced different scatterers from ' ...
         'the same seed, so the arms would not be comparable']);
end

A_mod  = build_dictionary(grid_mod,  u0, theta_m, f_hat_l, radar.fc, ...
    options.use_range_approx, 'mod-OMP dictionary');
A_base = build_dictionary(grid_base, u0, theta_m, f_hat_l, radar.fc, ...
    options.use_range_approx, 'baseline dictionary');

% A missed scatterer is charged the crossrange width of the imaged scene. That
% charge has to be identical for every algorithm or the comparison tilts
% toward whichever was given the narrower grid, so it comes from the wider one.
miss_penalty = max(grid_base.x_array) - min(grid_base.x_array);

fprintf('\n  mod-OMP  : %d ambiguity,  %3d x %3d grid (%5d atoms)\n', ...
    amb_if_mod_omp, grid_mod.Ny, grid_mod.Nx, size(A_mod,2));
fprintf('  baselines: %d ambiguities, %3d x %3d grid (%5d atoms)\n', ...
    amb_if_baselines, grid_base.Ny, grid_base.Nx, size(A_base,2));
fprintf('  Wx %.3f m, miss penalty %.3f m, %d scatterers over %d ambiguities\n\n', ...
    grid_mod.Wx, miss_penalty, sc.num_of_scatterers, sc.num_amb_having_scatterers);

%% monte carlo

algs      = ["mod_omp" "omp" "nomp" "promp"];
alg_label = ["mod-OMP" "OMP" "NOMP" "PROMP"];
enabled   = [options.execute_mod_omp, options.execute_omp, ...
             options.execute_nomp,    options.execute_promp];
n_alg     = numel(algs);

mc.error        = nan(n_iterations, n_alg);
mc.missed       = nan(n_iterations, n_alg);
mc.false_alarms = nan(n_iterations, n_alg);
mc.n_atoms      = nan(n_iterations, n_alg);
mc.seed         = nan(n_iterations, 1);

t_start = tic;
for iter = 1:n_iterations

    seed = base_seed + iter - 1;
    mc.seed(iter) = seed;
    rng(seed)

    % a fresh scatterer draw; the grid that comes back is the one already
    % built above, so only target_locations is taken from here
    [sc_iter, target_locations] = ...
        create_target_and_grid(sc_mod, t_m, radar, f_hat_l, options);

    % the measurement superposes the exact phase histories of the (off-grid)
    % scatterers, each column of `as` being one scatterer's response
    as = zeros(M*L, sc.num_of_scatterers);
    for k = 1:sc.num_of_scatterers
        as(:,k) = compute_atom(target_locations(k,1), target_locations(k,2), ...
            u0, theta_m, f_hat_l, radar.fc, options.use_range_approx);
    end
    y = as * (sc.target_magnitude * ones(sc.num_of_scatterers,1));

    % additive white complex Gaussian noise, and the residual level the
    % pursuits stop at: for per-sample variance sigma2 over N samples,
    % E{||e||^2} = N*sigma2
    options.residual_threshold = [];
    if ~isempty(sc.snr_db) && isfinite(sc.snr_db)
        sigma2 = mean(abs(y).^2) / 10^(sc.snr_db/10);
        y = y + sqrt(sigma2/2) * (randn(size(y)) + 1j*randn(size(y)));
        options.residual_threshold = sqrt(numel(y) * sigma2);
    end

    positions = cell(1, n_alg);
    pairings  = cell(1, n_alg);

    % ---- mod-OMP, on the one-ambiguity dictionary --------------------
    if options.execute_mod_omp
        [~, p_hat] = quietly(@mod_omp_vec, 2, ...
            y, A_mod, sc_iter.num_of_latent_scatterers, grid_mod, u0, ...
            theta_m, f_hat_l, radar.fc, sc.num_amb_having_scatterers, ...
            sc_iter, options);
        positions{1} = p_hat.';
    end

    % ---- OMP, on the full dictionary ---------------------------------
    if options.execute_omp
        x_hat_omp = quietly(@omp_vec, 1, ...
            y, A_base, sc_iter.num_of_latent_scatterers, options);
        positions{2} = extract_target_positions( ...
            reshape(x_hat_omp, grid_base.Ny, grid_base.Nx), ...
            grid_base.x_array, grid_base.y_array, ...
            sc.num_of_scatterers, 'none');
    end

    % ---- NOMP, on the full dictionary --------------------------------
    if options.execute_nomp
        [~, p_hat] = quietly(@nomp_vec, 2, ...
            y, A_base, sc.num_optimization_steps, sc.num_optimization_cycles, ...
            sc.num_of_scatterers, grid_base.xk, grid_base.yk, u0, ...
            theta_m, f_hat_l, radar.fc, options);
        positions{3} = p_hat.';
    end

    % ---- PROMP, on the full dictionary -------------------------------
    if options.execute_promp
        [~, p_hat] = quietly(@promp_vec, 2, ...
            y, A_base, sc.num_optimization_steps, sc.num_of_scatterers, ...
            grid_base.xk, grid_base.yk, u0, theta_m, f_hat_l, radar.fc, ...
            options);
        positions{4} = p_hat.';
    end

    % ---- score every algorithm against the same truth ----------------
    for a = 1:n_alg
        if isempty(positions{a})
            continue
        end
        [e, pairs, missed, false_alarms] = calculate_reconstruction_error( ...
            target_locations, positions{a}, miss_penalty);
        pairings{a} = pairs;

        mc.error(iter, a)        = e;
        mc.missed(iter, a)       = numel(missed);
        mc.false_alarms(iter, a) = numel(false_alarms);
        mc.n_atoms(iter, a)      = size(positions{a}, 1);
    end

    if options.plot_every > 0 && mod(iter, options.plot_every) == 0
        save_positions_figure(iter, target_locations, positions, ...
            mc.error(iter,:), pairings, alg_label, enabled, ...
            grid_mod.Wx, u0, positions_dir);
    end

    progress_bar('Monte Carlo', iter, n_iterations);
end
fprintf('\n');
elapsed = toc(t_start);

%% summary

fprintf('\n%d iterations in %.1f s (%.2f s each)\n\n', ...
    n_iterations, elapsed, elapsed/n_iterations);

fprintf('  %-9s %9s %9s %9s %9s %8s %8s %7s\n', ...
    'algorithm', 'mean', 'median', 'p90', 'worst', 'atoms', 'missed', 'best');
fprintf('  %-9s %9s %9s %9s %9s %8s %8s %7s\n', ...
    '', '[m]', '[m]', '[m]', '[m]', 'mean', 'mean', '%');
fprintf('  %s\n', repmat('-', 1, 76));

[~, winner] = min(mc.error, [], 2);
for a = 1:n_alg
    if ~enabled(a)
        continue
    end
    e = mc.error(~isnan(mc.error(:,a)), a);
    if isempty(e)
        continue
    end
    fprintf('  %-9s %9.4f %9.4f %9.4f %9.4f %8.2f %8.2f %6.1f%%\n', ...
        alg_label(a), mean(e), median(e), prctile(e,90), max(e), ...
        mean(mc.n_atoms(:,a),'omitnan'), mean(mc.missed(:,a),'omitnan'), ...
        100*mean(winner == a));
end
fprintf('\n');

mc.algs         = algs;
mc.alg_label    = alg_label;
mc.enabled      = enabled;
mc.sc           = sc;
mc.radar        = radar;
mc.options      = options;
mc.amb_if       = struct('mod_omp', amb_if_mod_omp, ...
                         'baselines', amb_if_baselines);
mc.plot_every   = options.plot_every;
mc.miss_penalty = miss_penalty;
mc.elapsed_s    = elapsed;
mc.grid_atoms   = struct('mod_omp', size(A_mod,2), ...
                         'baselines', size(A_base,2));

if ~isfolder(fileparts(results_file))
    mkdir(fileparts(results_file));
end
save(results_file, 'mc');
fprintf('  wrote %s\n', results_file);

%% empirical CDF of the position error

f = figure('Visible', 'off');
hold on
shown = {};
for a = 1:n_alg
    e = sort(mc.error(~isnan(mc.error(:,a)), a));
    if isempty(e)
        continue
    end
    plot(e, (1:numel(e))/numel(e), 'LineWidth', 2);
    shown{end+1} = char(alg_label(a)); %#ok<SAGROW>
end
hold off
set(gca, 'XScale', 'log')
grid on
xlabel('RMS position error [m]')
ylabel('empirical CDF')
title(sprintf(['mod-OMP (%d ambiguity) vs baselines (%d ambiguities), ' ...
    '%d runs'], amb_if_mod_omp, amb_if_baselines, n_iterations))
legend(shown, 'Location', 'southeast')
set(gca, 'FontSize', 12)

if ~isfolder(fileparts(figure_file))
    mkdir(fileparts(figure_file));
end
saveas(f, figure_file);
close(f)
fprintf('  wrote %s\n\n', figure_file);

%% ------------------------------------------------------------------ local

function A = build_dictionary(grid, u0, theta_m, f_hat_l, fc, ...
    use_range_approx, label)
% BUILD_DICTIONARY  One atom per grid node, as the batch scripts build it.

    K  = numel(grid.xk);
    ML = size(theta_m,1) * size(f_hat_l,2);
    A  = zeros(ML, K);

    for k = 1:K
        A(:,k) = compute_atom(grid.xk(k), grid.yk(k), u0, theta_m, ...
            f_hat_l, fc, use_range_approx);

        if mod(k, max(1,round(K/20))) == 0
            progress_bar(label, k, K);
        end
    end
    fprintf('\n');
end

function save_positions_figure(iter, truth, positions, err, pairings, ...
    labels, enabled, Wx, u0, out_dir)
% SAVE_POSITIONS_FIGURE  One row of position panels for a single iteration.
%
%   Every panel shows the same truth against one algorithm's estimates, on
%   shared axes, with the unambiguous crossrange bands marked. A scatterer
%   assigned to the wrong ambiguity is then obvious by eye: its pairing line
%   runs roughly Wx across the plot.

    n = numel(labels);

    % shared limits, so the panels can be read against one another
    xs = truth(:,1);
    ys = truth(:,2);
    for a = 1:n
        if enabled(a) && ~isempty(positions{a})
            xs = [xs; positions{a}(:,1)]; %#ok<AGROW>
            ys = [ys; positions{a}(:,2)]; %#ok<AGROW>
        end
    end
    padx = max(0.05*(max(xs) - min(xs)), 0.5);
    pady = max(0.05*(max(ys) - min(ys)), 0.5);
    xl = [min(xs) - padx, max(xs) + padx];
    yl = [min(ys) - pady, max(ys) + pady] + u0;

    f = figure('Visible', 'off', 'Position', [0 0 1900 480]);

    for a = 1:n
        subplot(1, n, a)
        hold on

        % the ambiguity band edges sit at half-integer multiples of Wx
        for edge = (-3.5:1:3.5) * Wx
            if edge > xl(1) && edge < xl(2)
                plot([edge edge], yl, '--', 'Color', [0.82 0.82 0.82], ...
                    'HandleVisibility', 'off');
            end
        end

        P = positions{a};

        % pairing lines first, so the markers draw on top of them
        if enabled(a) && ~isempty(P) && ~isempty(pairings{a})
            for j = 1:size(pairings{a}, 1)
                it = pairings{a}(j,1);
                ie = pairings{a}(j,2);
                plot([truth(it,1) P(ie,1)], [truth(it,2) P(ie,2)] + u0, '-', ...
                    'Color', [0.65 0.65 0.65], 'LineWidth', 1, ...
                    'HandleVisibility', 'off');
            end
        end

        h = plot(truth(:,1), truth(:,2) + u0, 'o', ...
            'MarkerEdgeColor', [0.85 0.33 0.10], 'MarkerSize', 9, ...
            'LineWidth', 1.4, 'DisplayName', 'true');

        if enabled(a) && ~isempty(P)
            h(end+1) = plot(P(:,1), P(:,2) + u0, 'x', ...
                'MarkerEdgeColor', [0 0.45 0.74], 'MarkerSize', 11, ...
                'LineWidth', 1.6, 'DisplayName', 'estimated'); %#ok<AGROW>
        end
        hold off

        xlim(xl); ylim(yl); grid on; box on
        set(gca, 'YDir', 'reverse', 'FontSize', 11)
        xlabel('crossrange [m]')
        if a == 1
            ylabel('range [m]')
            legend(h, 'Location', 'southoutside', 'Orientation', 'horizontal', ...
                'FontSize', 9, 'Box', 'off');
        end
        title(sprintf('%s  (%.3f m)', labels(a), err(a)), 'FontSize', 12)
    end

    sgtitle(sprintf('iteration %d', iter), 'FontSize', 14);

    if ~isfolder(out_dir)
        mkdir(out_dir);
    end
    saveas(f, fullfile(out_dir, sprintf('positions_iter%04d.png', iter)));
    close(f)
end

function varargout = quietly(fcn, n_out, varargin)
% QUIETLY  Call FCN, discarding whatever it prints.
%
%   The pursuits report their progress on every call, which is useful for a
%   single run and unreadable over hundreds. Their behaviour is untouched;
%   this only swallows the text.

    varargout = cell(1, n_out);
    txt = evalc('[varargout{1:n_out}] = fcn(varargin{:});'); %#ok<NASGU>
end
