function [alpha_hat, p_hat] = nomp_newton(s, A, Rs, Rc, sparsity, xk, yk, ...
    u0, theta_m, f_hat_l, fc, options, newton_fn)
% NOMP_NEWTON  nomp_vec with the Newton routine passed in.
%
%   Identical to nomp_vec (detection on the grid, Newton refinement of the new
%   atom, Rc cyclic refinements of every atom so far, least squares) except
%   that every Newton call goes through NEWTON_FN, e.g. @newton_method (which
%   reproduces nomp_vec) or @newton_method_exact. The step-history output of
%   nomp_vec is not kept.

    ura       = options.use_range_approx;
    alpha_hat = zeros(sparsity, 1);
    p_hat     = zeros(2, sparsity);

    % Section IV-A stops every algorithm when the signal residual reaches the
    % noise level; Inf when no threshold is set, leaving the sparsity cap
    tau    = residual_stop_threshold(options);
    n_kept = sparsity;

    r = s;
    for iatom = 1:sparsity

        %---------------- detection ----------------
        [~, k] = max(abs(A' * r), [], "all");

        %---------------- refinement ----------------
        p_hat(:,iatom) = newton_fn(r, xk(k), yk(k), u0, theta_m, f_hat_l, fc, Rs, ura);

        AS = build_AS(p_hat(:,1:iatom), u0, theta_m, f_hat_l, fc, ura);
        alpha_hat(1:iatom) = AS \ s;

        %------------- cyclic refinement -------------
        if sparsity ~= 1
            for icycle = 1:Rc
                for jatom = 1:iatom
                    idx = 1:iatom; idx(jatom) = [];
                    AS = build_AS(p_hat(:,idx), u0, theta_m, f_hat_l, fc, ura);
                    rj = s - AS * alpha_hat(idx);
                    p_hat(:,jatom) = newton_fn(rj, p_hat(1,jatom), p_hat(2,jatom), ...
                        u0, theta_m, f_hat_l, fc, Rs, ura);
                end
                AS = build_AS(p_hat(:,1:iatom), u0, theta_m, f_hat_l, fc, ura);
                alpha_hat(1:iatom) = AS \ s;
            end
        end
        r = s - AS * alpha_hat(1:iatom);

        % the residual has reached the noise level
        if norm(r) <= tau
            n_kept = iatom;
            break
        end
    end

    p_hat     = p_hat(:, 1:n_kept);
    alpha_hat = alpha_hat(1:n_kept);
end
