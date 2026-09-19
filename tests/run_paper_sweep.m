function res = run_paper_sweep(n_iter, progress)
% RUN_PAPER_SWEEP  The Monte Carlo comparison behind the paper's numbers.
%
%   res = RUN_PAPER_SWEEP(n_iter) draws n_iter scatterer sets (default 500)
%   from PAPER_SWEEP_SETTINGS and runs every algorithm arm on each. mod-OMP
%   and its variants see an image former spanning one ambiguity; OMP and NOMP
%   are handed a dictionary spanning sc.num_amb_in_image_former_baselines, the
%   split isar_monte_carlo.m makes. That is the comparison: mod-OMP recovers a
%   scatterer outside the imaged band by testing its ambiguous doppelgangers,
%   so it never needs the wider dictionary.
%
%   The arms, in order:
%
%     1  NOMP (paper)              nomp_vec, newton_method as published
%     2  mod-OMP (repo)            mod_omp_vec, the repository's own version
%     3  mod-OMP orig, new offsets mod_omp_newton + newton_method
%     4  mod-OMP proj              projected Hessian, no backtracking
%     5  mod-OMP proj + cyclic     the same with Rc cyclic refinements
%     6  OMP                       on-grid reference
%     7  mod-OMP exact             exact Hessian + backtracking
%     8  mod-OMP exact + cyclic    the same with Rc cyclic refinements
%     9  mod-OMP 2D exact          exact, candidates also seeded in range
%    10  mod-OMP 2D exact + cyclic the same with Rc cyclic refinements
%    11  NOMP orig, wrapper        nomp_newton + newton_method; must equal 1
%    12  NOMP exact                nomp_newton + newton_method_exact
%    13  mod-OMP 2D rand + cyclic  random 2-D candidates, screened, exact
%
%   RES fields:
%     E          [n_iter x n_arm] RMS position error of each draw
%     labels     arm names
%     T          [1 x n_arm] total seconds spent in each arm
%     arm_dict   1 = searched A_mod (one ambiguity), 2 = A_base
%     time_A_mod, time_A_base, K_mod, K_base   dictionary build times and sizes
%     Wx, range_offset, settings, n_iter, base_seed
%
%   Timings are wall clock and machine dependent, so TEST_PAPER_SWEEP checks E
%   and reports times without asserting on them.
%
%   See also TEST_PAPER_SWEEP, PAPER_SWEEP_SETTINGS, ISAR_MONTE_CARLO.

    if nargin < 1 || isempty(n_iter),   n_iter = 500;    end
    if nargin < 2 || isempty(progress), progress = true; end

    [sc, options, radar] = paper_sweep_settings();
    base_seed = options.seed;

    t_m   = (0:(1/radar.prf):(radar.T  - 1/radar.prf)).';
    t_hat = (0:(1/radar.fs):(radar.Tp - 1/radar.fs)).';
    M = numel(t_m); L = numel(t_hat);
    df_l = radar.fs/L; f_hat_l = (-L/2)*df_l:df_l:(L/2-1)*df_l;
    fc = radar.fc; ura = options.use_range_approx;

    % the two image formers, from the same draw
    sc_mod  = sc; sc_mod.num_amb_in_image_former  = sc.num_amb_in_image_former;
    sc_base = sc; sc_base.num_amb_in_image_former = sc.num_amb_in_image_former_baselines;
    rng(base_seed); [sc_mod, tgt_mod, grid_mod, theta_m, u0] = ...
        create_target_and_grid(sc_mod, t_m, radar, f_hat_l, options);
    rng(base_seed); [~, tgt_base, grid_base] = ...
        create_target_and_grid(sc_base, t_m, radar, f_hat_l, options);

    % the draw uses Wx and the range extent, neither of which depends on the
    % image former span, so the same seed has to reproduce it
    if ~isequal(tgt_mod, tgt_base)
        error('run_paper_sweep:targetMismatch', ...
            ['the two image former spans produced different scatterers from ' ...
             'the same seed, so the arms would not be comparable']);
    end

    t_dict = tic;
    [A_mod, walk_mod] = build_dict(grid_mod, u0, theta_m, f_hat_l, fc, ura);
    time_A_mod = toc(t_dict);
    t_dict = tic;
    A_base = build_dict(grid_base, u0, theta_m, f_hat_l, fc, ura);
    time_A_base = toc(t_dict);

    % the 2-D variants seed candidates over the largest range walk across the
    % aperture, scaled by the ambiguity span, as isar_testing_v3.m sets it
    grid_mod.range_offset = sc.num_amb_having_scatterers * walk_mod;

    % a missed scatterer is charged the crossrange width of the imaged scene;
    % the charge must be identical for every arm, so it comes from the wider
    % grid, or the comparison tilts toward the narrower one
    miss_penalty = max(grid_base.x_array) - min(grid_base.x_array);

    K    = sc.num_of_scatterers;
    Namb = sc.num_amb_having_scatterers;
    Rs   = sc.num_optimization_steps;
    Rc   = sc.num_optimization_cycles;

    labels = ["NOMP (paper)", "mod-OMP (repo)", "mod-OMP orig, new offsets", ...
              "mod-OMP proj", "mod-OMP proj + cyclic", "OMP", ...
              "mod-OMP exact", "mod-OMP exact + cyclic", ...
              "mod-OMP 2D exact", "mod-OMP 2D exact + cyclic", ...
              "NOMP orig, wrapper", "NOMP exact", "mod-OMP 2D rand + cyclic"];
    arm_dict = [2 1 1 1 1 2 1 1 1 1 2 2 1];

    E = nan(n_iter, numel(labels));
    T = zeros(1, numel(labels));

    for it = 1:n_iter
        rng(base_seed + it - 1)
        [sc_it, tl] = create_target_and_grid(sc_mod, t_m, radar, f_hat_l, options);

        as = zeros(M*L, K);
        for k = 1:K
            as(:,k) = compute_atom(tl(k,1), tl(k,2), u0, theta_m, f_hat_l, fc, ura);
        end
        y  = as * (sc.target_magnitude * ones(K,1));
        ns = sc_it.num_of_latent_scatterers;
        P  = cell(1, numel(labels));

        mod_args  = {y, A_mod, ns, grid_mod, u0, theta_m, f_hat_l, fc, Namb, sc_it, options};
        nomp_args = {y, A_base, Rs, Rc, K, grid_base.xk, grid_base.yk, u0, theta_m, f_hat_l, fc, options};

        tic; [~,p] = quiet(@nomp_vec, nomp_args{:});                                   P{1}  = p.'; T(1)  = T(1)  + toc;
        tic; [~,p] = quiet(@mod_omp_vec, mod_args{:});                                 P{2}  = p.'; T(2)  = T(2)  + toc;
        tic; [~,p] = quiet(@mod_omp_newton, mod_args{:}, 0,  @newton_method);          P{3}  = p.'; T(3)  = T(3)  + toc;
        tic; [~,p] = quiet(@mod_omp_newton, mod_args{:}, 0,  @newton_method_projected); P{4} = p.'; T(4)  = T(4)  + toc;
        tic; [~,p] = quiet(@mod_omp_newton, mod_args{:}, Rc, @newton_method_projected); P{5} = p.'; T(5)  = T(5)  + toc;
        tic; xh = omp_vec(y, A_base, ns, options);
        P{6} = extract_target_positions(reshape(xh, grid_base.Ny, grid_base.Nx), ...
            grid_base.x_array, grid_base.y_array, K, 'none');                          T(6)  = T(6)  + toc;
        tic; [~,p] = quiet(@mod_omp_newton, mod_args{:}, 0,  @newton_method_exact);    P{7}  = p.'; T(7)  = T(7)  + toc;
        tic; [~,p] = quiet(@mod_omp_newton, mod_args{:}, Rc, @newton_method_exact);    P{8}  = p.'; T(8)  = T(8)  + toc;
        tic; [~,p] = quiet(@mod_omp_newton_2d, mod_args{:}, 0,  @newton_method_exact); P{9}  = p.'; T(9)  = T(9)  + toc;
        tic; [~,p] = quiet(@mod_omp_newton_2d, mod_args{:}, Rc, @newton_method_exact); P{10} = p.'; T(10) = T(10) + toc;
        tic; [~,p] = quiet(@nomp_newton, nomp_args{:}, @newton_method);                P{11} = p.'; T(11) = T(11) + toc;
        tic; [~,p] = quiet(@nomp_newton, nomp_args{:}, @newton_method_exact);          P{12} = p.'; T(12) = T(12) + toc;
        tic; [~,p] = quiet(@mod_omp_newton_2d_rand, mod_args{:}, Rc, @newton_method_exact); P{13} = p.'; T(13) = T(13) + toc;

        for a = 1:numel(labels)
            E(it,a) = calculate_reconstruction_error(tl, P{a}, miss_penalty);
        end

        if progress && mod(it, 50) == 0
            fprintf('    %d/%d\n', it, n_iter);
        end
    end

    res = struct('E', E, 'labels', labels, 'T', T, 'arm_dict', arm_dict, ...
        'time_A_mod', time_A_mod, 'time_A_base', time_A_base, ...
        'K_mod', size(A_mod,2), 'K_base', size(A_base,2), ...
        'Wx', grid_mod.Wx, 'range_offset', grid_mod.range_offset, ...
        'miss_penalty', miss_penalty, 'n_iter', n_iter, 'base_seed', base_seed, ...
        'settings', struct('sc', sc, 'options', options, 'radar', radar));
end

function [A, walk] = build_dict(g, u0, th, f, fc, ura)
% one atom per grid node; WALK is the largest range change over the aperture
% (first sample to last) of any atom
    A = zeros(numel(th)*numel(f), numel(g.xk));
    walk = 0;
    for k = 1:numel(g.xk)
        [A(:,k), ~, rk] = compute_atom(g.xk(k), g.yk(k), u0, th, f, fc, ura);
        walk = max(walk, abs(rk(1) - rk(end)));
    end
end

function varargout = quiet(fcn, varargin)
% call FCN, discarding whatever it prints
    varargout = cell(1,2);
    evalc('[varargout{1:2}] = fcn(varargin{:});');
end
