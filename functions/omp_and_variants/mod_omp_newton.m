function [alpha_hat, p_hat, cand] = mod_omp_newton(...
    s,...
    A,...
    sparsity,...
    grid,...
    u0, ...
    theta_m,...
    f_hat_l,...
    fc,...
    num_of_amb,...
    sc,...
    options,...
    Rc,...
    newton_fn)
% MOD_OMP_NEWTON  mod_omp_vec ('newtons_and_offset') with the Newton routine
% passed in, plus optional NOMP-style cyclic refinement.
%
%   Identical to mod_omp_vec + check_ambiguity_and_apply_newtons_with_offsets
%   except that every Newton call goes through NEWTON_FN (e.g.
%   @newton_method, @newton_method_projected or @newton_method_exact), and,
%   when Rc > 0, every atom selected so far is re-refined Rc times against the
%   measurement minus the other atoms after each new selection (the loop in
%   nomp_vec). Rc = 0 and newton_fn = @newton_method reproduces mod_omp_vec.
%
%   CAND (optional third output) records each iteration's doppelganger search,
%   for plotting where the candidate atoms land. It is only filled when
%   requested. Each element has
%
%       r         [ML x 1] residual the search was run against
%       l         index of the atom detected on the one-ambiguity grid
%       x0, y0    (m) that grid atom's position
%       amb       [1 x Nc] ambiguity index of each candidate
%       seed      [2 x Nc] Newton starting point of each candidate
%       p         [2 x Nc] where Newton ended
%       path      {1 x Nc} accepted Newton iterates, [2 x n] each
%       corr      [1 x Nc] |a(p)' r| / (||a|| ||r||) at the end point
%       best      index of the candidate that was selected
%
%   (Candidate positions are before cyclic refinement; p_hat is after.)

    record = nargout >= 3;

    ML  = size(A,1);
    Rs  = sc.num_optimization_steps;
    ura = options.use_range_approx;

    ii        = 0:(num_of_amb-1);
    amb_index = sort(ceil(ii/2) .* (-1).^ii);
    res       = grid.cross_range_pixel_res;

    % the crossrange offset defining the neighborhood is
    % dependent on the range of Doppler ambiguities
    Wx_offset_range = abs(grid.Wx_m(end) - grid.Wx_m(1));
    num_offsets = floor(Wx_offset_range / res);
    x_offsets = -floor(num_offsets/2)*res:res:floor(num_offsets/2)*res;

    % x_offsets = (-sc.num_offsets_pixels*res):res:(sc.num_offsets_pixels*res);
    [AMB, OFF] = ndgrid(amb_index, x_offsets);
    AMB = AMB.'; OFF = OFF.';                 % ambiguity-major ordering
    Nc  = numel(AMB);

    selected_idx = [];
    As        = zeros(ML, sparsity);
    p_hat     = zeros(2, sparsity);
    alpha_hat = zeros(sparsity, 1);
    cand      = struct('r', {}, 'l', {}, 'x0', {}, 'y0', {}, 'amb', {}, ...
                       'seed', {}, 'p', {}, 'path', {}, 'corr', {}, 'best', {});

    % Section IV-A stops every algorithm when the signal residual reaches the
    % noise level; Inf when no threshold is set, leaving the sparsity cap
    tau    = residual_stop_threshold(options);
    n_kept = sparsity;

    r = s;
    for i = 1:sparsity

        %------------- detection on the one-ambiguity grid ----------------
        c = A' * r;
        c(selected_idx) = 0;
        [~, l] = max(abs(c), [], "all");
        selected_idx = [selected_idx; l]; %#ok<AGROW>
        x0 = grid.xk(l); y0 = grid.yk(l);

        %------------- doppelganger search with Newton --------------------
        seed = [x0 + AMB(:).'*grid.Wx + OFF(:).'; repmat(y0, 1, Nc)];
        pend = zeros(2, Nc); corr = zeros(1, Nc); path = cell(1, Nc);
        best = 0; best_c = -inf;
        for ic = 1:Nc
            if record
                [p, path{ic}] = newton_fn(r, seed(1,ic), seed(2,ic), u0, theta_m, f_hat_l, fc, Rs, ura);
            else
                p = newton_fn(r, seed(1,ic), seed(2,ic), u0, theta_m, f_hat_l, fc, Rs, ura);
            end
            ak = compute_atom(p(1), p(2), u0, theta_m, f_hat_l, fc, ura);
            ck = abs(ak' * r);
            pend(:,ic) = p;
            corr(ic)   = ck / (norm(ak) * norm(r));
            if ck > best_c
                best_c = ck; best = ic;
                As(:,i) = ak; p_hat(:,i) = p;
            end
        end

        if record
            cand(i) = struct('r', r, 'l', l, 'x0', x0, 'y0', y0, ...
                'amb', AMB(:).', 'seed', seed, 'p', pend, 'path', {path}, ...
                'corr', corr, 'best', best);
        end

        alpha_hat(1:i) = As(:,1:i) \ s;

        %------------- cyclic refinement (as in nomp_vec) -----------------
        if sparsity ~= 1
            for icycle = 1:Rc
                for j = 1:i
                    idx = 1:i; idx(j) = [];
                    rj = s - As(:,idx) * alpha_hat(idx);
                    p_hat(:,j) = newton_fn(rj, p_hat(1,j), p_hat(2,j), ...
                        u0, theta_m, f_hat_l, fc, Rs, ura);
                    As(:,j) = compute_atom(p_hat(1,j), p_hat(2,j), u0, ...
                        theta_m, f_hat_l, fc, ura);
                end
                alpha_hat(1:i) = As(:,1:i) \ s;
            end
        end
        r = s - As(:,1:i) * alpha_hat(1:i);

        % the residual has reached the noise level
        if norm(r) <= tau
            n_kept = i;
            break
        end
    end

    % an early stop leaves the tail of the preallocated arrays untouched
    p_hat     = p_hat(:, 1:n_kept);
    alpha_hat = alpha_hat(1:n_kept);
end
